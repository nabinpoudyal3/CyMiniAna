"""
Created:         6 March 2018
Last Updated:    6 March 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Make data/mc plots

To run:
 python python/runDataMC_ttbarAC.py --hists <histograms.txt> -o <output_path>


/store/group/lpctop/ttbarAC/flatNtuples/ttbarAC_skim-v0.2/
"""
import sys
import ROOT
from time import strftime
from array import array
from argparse import ArgumentParser

from Analysis.CyMiniAna.hepPlotter.hepPlotterDataMC import HepPlotterDataMC
from Analysis.CyMiniAna.hepPlotter.hepPlotterSystematics import HepPlotterSystematics
import Analysis.CyMiniAna.hepPlotter.hepPlotterTools as hpt
import Analysis.CyMiniAna.hepPlotter.hepPlotterLabels as hpl
import Analysis.CyMiniAna.util as util


contain = [
'NONE',
'QONLY',
'BONLY',
'BQ',
'W',
'FULL',
]
metadata = hpt.getMetadata()  # use default metadatafile


def getHistograms(files,histograms,sel):
    """Aggregate histograms from many files"""
    pd    = ''
    name  = ''
    hists = {}

    for fi,file in enumerate(files):
        f = ROOT.TFile.Open(file)

        if not fi:
            pd   = util.getPrimaryDataset(f)
            name = metadata[pd].sampleType     # compare primary dataset with metadatafile

        for h in histograms:
            if name == 'ttbar' and any(ic in h for ic in ['ljet_BEST','ljet_pt','ljet_SDmass']):
                for c in contain:
                    hc = h.replace('_'+sel,'_'+c+'_'+sel)
                    h_temp = getattr(f,hc)
                    h_temp.SetDirectory(0)
                    try:
                        hists[hc].Add( h_temp )
                    except KeyError:
                        hists[hc] = h_temp
            else:
                h_temp = getattr(f,h)
                h_temp.SetDirectory(0)
                try:
                    hists[h].Add( h_temp )
                except KeyError:
                    hists[h] = h_temp

    return {"hists":hists,"primaryDataset":pd,"name":name}


x_labels      = hpl.variable_labels()
sample_labels = hpl.sample_labels()
extralabel    = {"ejets":sample_labels['ejets'].label,"mujets":sample_labels['mujets'].label,'afb':r'A$_\text{FB}$'}
extralabel["cwola"] = "CWoLa"
extralabel["cwolaejets"]  = "CWoLa "+extralabel["ejets"]
extralabel["cwolamujets"] = "CWoLa "+extralabel["mujets"]

today  = strftime("%d%b%Y")
parser = ArgumentParser(description="DataMC Plotter")

parser.add_argument('--hists', action='store',dest='listOfHists',
                    default='config/listOfHists.txt',
                    help='Name of file that contains histograms to plot')
parser.add_argument('-o','--outpath', action='store',default='plots/'+today+"/",
                    dest='outpath',
                    help='Directory for storing output plots')
parser.add_argument('--selection',action='store',default='ejets',
                    dest='selection',
                    help='Plot results from a specific selection')
results = parser.parse_args()


outpath       = results.outpath
selection     = results.selection
histograms    = util.file2list(results.listOfHists)
histograms    = [i.format(selection) for i in histograms]

samples = [
'diboson_ww',
'diboson_wz',
'diboson_zz',
'zjets',
'singletop-s',
'singletop-t',
'singletop-tW',
'singletop-tbar',
'singletop-tbarW',
#'ttbar-ext',
'ttbar',
'wjets1',
'wjets2',
'wjets3',
'wjets4',
'data'
]
dataname   = 'SingleElectron' if 'ejets' in selection else 'SingleMuon'
sample_dir = 'config/cwola_samples'

filelists = dict( (k,util.file2list('{0}/{1}.txt'.format(sample_dir,k))) for k in samples if k!='data')

data_files  = []
for d in ['B','C','D','E','F','G','H','Hv3']:
    data_files += util.file2list("{0}/{1}{2}.txt".format(sample_dir,dataname,d))
filelists['data'] = data_files

# Load histograms
h_hists_dict = {}
for sample in samples:
    print " Load histograms from ",sample
    h_hists_dict[sample] = getHistograms(filelists[sample],histograms,selection)


## Make the plots
for histogram in histograms:

    histogramName = histogram.replace("_"+selection,"")
    if histogramName.startswith("h_"): histogramName = histogramName[2:]

    if histogram.startswith("h_mu_") and "ejets" in selection: continue
    if histogram.startswith("h_el_") and ("mujets" in selection or selection=="afb"): continue

    print "  :: Plotting "+histogram

    ## setup histogram
    hist = HepPlotterDataMC()

    xlimits   = None
    legendLoc = 1

    hist.drawStatUncertainty = True      # draw stat uncertainty separately
    hist.drawSystUncertainty = False     # draw syst uncertainty separately
    hist.drawStatSystUncertainty = False # draw stat+syst uncertainty
    hist.legendLoc   = legendLoc
    hist.stackSignal = False
    hist.blind       = False
    hist.xlim        = xlimits
    hist.ratio_plot  = True        # plot a ratio of things [Data/MC]
    hist.ratio_type  = "ratio"     # "ratio"
    hist.stacked     = True        # stack plots
    hist.rebin       = x_labels[histogramName].binning  # rebin per histogram
    hist.logplot     = {"y":False,"x":False,"data":False}
    hist.x_label     = x_labels[histogramName].label
    hist.y_label     = "Events"
    hist.y_ratio_label = "Data/MC"
    hist.lumi          = '35.9'
    hist.CMSlabel    = 'top left'  # 'top left', 'top right'; hack code for something else
    hist.numLegendColumns = 1
    hist.CMSlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
    hist.format           = 'pdf'       # file format for saving image
    hist.saveAs           = outpath+"/datamc_ttbarAC_splitTtbar_{0}_{1}".format(selection,histogramName) # save figure with name

    hist.extra_text.Add(extralabel[selection],coords=[0.03,0.90])

    hist.initialize()

    ## -- Add the results from each sample -- group some together
    # diboson
    h_hist_diboson = h_hists_dict['diboson_ww']['hists'][histogram]
    h_hist_diboson.Add( h_hists_dict['diboson_wz']['hists'][histogram] )
    h_hist_diboson.Add( h_hists_dict['diboson_zz']['hists'][histogram] )
    hist.Add(h_hist_diboson, name="diboson", sampleType="background", systematics=None)

    # zjets
    hist.Add(h_hists_dict['zjets']['hists'][histogram], name="zjets", sampleType="background", systematics=None)

    # singletop
    h_hist_singletop = h_hists_dict['singletop-s']['hists'][histogram]
    h_hist_singletop.Add( h_hists_dict['singletop-t']['hists'][histogram] )
    h_hist_singletop.Add( h_hists_dict['singletop-tW']['hists'][histogram] )
    h_hist_singletop.Add( h_hists_dict['singletop-tbar']['hists'][histogram] )
    h_hist_singletop.Add( h_hists_dict['singletop-tbarW']['hists'][histogram] )
    hist.Add(h_hist_singletop, name="singletop", sampleType="background", systematics=None)

    # w+jets
    h_hist_wjets = h_hists_dict['wjets1']['hists'][histogram]
    h_hist_wjets.Add( h_hists_dict['wjets2']['hists'][histogram] )
    h_hist_wjets.Add( h_hists_dict['wjets3']['hists'][histogram] )
    h_hist_wjets.Add( h_hists_dict['wjets4']['hists'][histogram] )
    hist.Add(h_hist_wjets, name="wjets", sampleType="background", systematics=None)

    # ttbar
    if any(ic in histogram for ic in ['ljet_BEST','ljet_pt','ljet_SDmass']):
        h_hist_ttbarFULL  = h_hists_dict['ttbar']['hists'][histogram.replace(selection,"FULL_"+selection)]
        h_hist_ttbarQB    = h_hists_dict['ttbar']['hists'][histogram.replace(selection,"BQ_"+selection)]
        h_hist_ttbarW     = h_hists_dict['ttbar']['hists'][histogram.replace(selection,"W_"+selection)]
        h_hist_ttbarOTHER = h_hists_dict['ttbar']['hists'][histogram.replace(selection,"BONLY_"+selection)]
        h_hist_ttbarOTHER.Add( h_hists_dict['ttbar']['hists'][histogram.replace(selection,"QONLY_"+selection)] )
        h_hist_ttbarOTHER.Add( h_hists_dict['ttbar']['hists'][histogram.replace(selection,"NONE_"+selection)] )
        hist.Add(h_hist_ttbarOTHER, name="ttbar_OTHER", sampleType="background", systematics=None)
        hist.Add(h_hist_ttbarQB,    name="ttbar_BQ",    sampleType="background", systematics=None)
        hist.Add(h_hist_ttbarW,     name="ttbar_W",     sampleType="background", systematics=None)
        hist.Add(h_hist_ttbarFULL,  name="ttbar_FULL",  sampleType="background", systematics=None)
    else:
        hist.Add(h_hists_dict['ttbar']['hists'][histogram], name="ttbar", sampleType="background", systematics=None)

    # data
    hist.Add(h_hists_dict['data']['hists'][histogram], name="data", sampleType="data", systematics=None)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "     ++ Saved plot to "+hist.saveAs+"\n"


## THE END ##
