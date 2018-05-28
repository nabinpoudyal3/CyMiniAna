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
from array import array
from argparse import ArgumentParser

from Analysis.CyMiniAna.hepPlotter.hepPlotterDataMC import HepPlotterDataMC
from Analysis.CyMiniAna.hepPlotter.hepPlotterSystematics import HepPlotterSystematics
import Analysis.CyMiniAna.hepPlotter.hepPlotterTools as hpt
import Analysis.CyMiniAna.hepPlotter.hepPlotterLabels as hpl
import Analysis.CyMiniAna.util as util

metadata = hpt.getMetadata()  # use default metadatafile


def getHistograms(files,histograms):
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
            try:
                h_temp = getattr(f,h)
                h_temp.SetDirectory(0)
                hists[h].Add( h_temp )
            except KeyError:
                hists[h] = getattr(f,h)        # retrieve the histogram
                hists[h].SetDirectory(0)

    return {"hists":hists,"primaryDataset":pd,"name":name}


x_labels      = hpl.variable_labels()
sample_labels = hpl.sample_labels()
extralabel    = {"ejets":sample_labels['ejets'].label,"mujets":sample_labels['mujets'].label}
contain = [
'NONE',
'QONLY',
'BONLY',
'BQ',
'W',
'FULL',
]


parser = ArgumentParser(description="DataMC Plotter")

parser.add_argument('--hists', action='store',dest='listOfHists',
                    default='config/listOfHists.txt',
                    help='Name of file that contains histograms to plot')
parser.add_argument('-o','--outpath', action='store',default='plots/23May2018/',
                    dest='outpath',
                    help='Directory for storing output plots')
parser.add_argument('--selection',action='store',default='mujets',
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
'ttbar-ext',
#'ttbar',
'wjets1',
'wjets2',
'wjets3',
'wjets4',
'data'
]
dataname = 'SingleElectron' if selection == 'ejets' else 'SingleMuon'

data_files  = []
for d in ['C','D','E','F','G','H','Hv3']:
    data_files += util.file2list("config/cyminiana_samples/{0}{1}.txt".format(dataname,d))

filelists = dict( (k,util.file2list('config/cyminiana_samples/{0}.txt'.format(k))) for k in samples if k!='data')
filelists['data'] = data_files

# Load histograms
h_hists_dict = {}
for sample in samples:
    print " Load histograms from ",sample
    h_hists_dict[sample] = getHistograms(filelists[sample],histograms)


## Make the plots
for histogram in histograms:

    histogramName = histogram.replace("_"+selection,"")
    if histogramName.startswith("h_"): histogramName = histogramName[2:]

    if histogram.startswith("h_mu_") and selection=="ejets": continue
    if histogram.startswith("h_el_") and selection=="mujets": continue

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
    hist.saveAs           = outpath+"/datamc_ttbarAC_{0}_{1}".format(selection,histogramName) # save figure with name

    hist.extra_text.Add(extralabel[selection],coords=[0.03,0.90])

    hist.initialize()

    ## -- Add the data from each sample -- group some together
    h_hist_diboson = None
    h_hist_zjets   = None
    h_hist_singletop = None
    h_hist_wjets   = None
    h_hist_ttbar   = None
    h_hist_data    = None

    for sample in samples:

        # diboson
        if h_hists_dict[sample]['name'] == 'diboson':
            try:
                h_hist_diboson.Add( h_hists_dict[sample]['hists'][histogram] )
            except AttributeError:
                h_hist_diboson = h_hists_dict[sample]['hists'][histogram]

        # zjets
        if h_hists_dict[sample]['name'] == 'zjets':
            try:
                h_hist_zjets.Add( h_hists_dict[sample]['hists'][histogram] )
            except AttributeError:
                h_hist_zjets = h_hists_dict[sample]['hists'][histogram]

        # singletop
        if h_hists_dict[sample]['name'] == 'singletop':
            try:
                h_hist_singletop.Add( h_hists_dict[sample]['hists'][histogram] )
            except AttributeError:
                h_hist_singletop = h_hists_dict[sample]['hists'][histogram]

        # w+jets
        if h_hists_dict[sample]['name'] == 'wjets':
            try:
                h_hist_wjets.Add( h_hists_dict[sample]['hists'][histogram] )
            except AttributeError:
                h_hist_wjets = h_hists_dict[sample]['hists'][histogram]

        # ttbar
        if h_hists_dict[sample]['name'] == 'ttbar':
            try:
                h_hist_ttbar.Add( h_hists_dict[sample]['hists'][histogram] )
            except AttributeError:
                h_hist_ttbar = h_hists_dict[sample]['hists'][histogram]

        # data
        if h_hists_dict[sample]['name'] in ['data','ejets','mujets']:
            try:
                h_hist_data.Add( h_hists_dict[sample]['hists'][histogram] )
            except AttributeError:
                h_hist_data = h_hists_dict[sample]['hists'][histogram]

    hist.Add(h_hist_diboson, name="diboson", sampleType="background", systematics=None)
    hist.Add(h_hist_zjets,   name="zjets",   sampleType="background", systematics=None)
    hist.Add(h_hist_singletop, name="singletop", sampleType="background", systematics=None)
    hist.Add(h_hist_wjets, name="wjets",     sampleType="background", systematics=None)
    hist.Add(h_hist_ttbar, name="ttbar",     sampleType="background", systematics=None)
    hist.Add(h_hist_data,  name="data",      sampleType="data",       systematics=None)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "     ++ Saved plot to "+hist.saveAs+"\n"


## THE END


"""
#    if 'ljet' in histo:
#        h_hist_ttbarFULL  = h_hists_ttbar[histo.replace(selection,"FULL_{0}".format(selection))]
#        h_hist_ttbarQB    = h_hists_ttbar[histo.replace(selection,"BQ_{0}".format(selection))]
#        h_hist_ttbarW     = h_hists_ttbar[histo.replace(selection,"W_{0}".format(selection))]
#        h_hist_ttbarOTHER = h_hists_ttbar[histo.replace(selection,"OTHER_{0}".format(selection))]
#        hist.Add(h_hist_ttbarOTHER, name="ttbar_OTHER", sampleType="background", systematics=None)
#        hist.Add(h_hist_ttbarQB,    name="ttbar_BQ",    sampleType="background", systematics=None)
#        hist.Add(h_hist_ttbarW,     name="ttbar_W",     sampleType="background", systematics=None)
#        hist.Add(h_hist_ttbarFULL,  name="ttbar_FULL",  sampleType="background", systematics=None)
#    else:
#        h_hist_ttbar = h_hists_ttbar[histo]
#        hist.Add(h_hist_ttbar, name="ttbar",     sampleType="background", systematics=None)
"""
