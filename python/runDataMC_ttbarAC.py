"""
Created:         6 March 2018
Last Updated:    6 March 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Make data/mc plots

To run:
 python python/runDataMC_ttbarAC.py --hists <histogramNames.txt> -o <output_path>


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
parser.add_argument('-o','--outpath', action='store',default='plots/01May2018/datamc/',
                    dest='outpath',
                    help='Directory for storing output plots')
results = parser.parse_args()

listOfHists = results.listOfHists
outpath     = results.outpath
histograms  = util.file2list(listOfHists)


metadata      = hpt.getMetadata()
x_labels      = hpl.variable_labels()
sample_labels = hpl.sample_labels()
selection     = 'ejets'

samples = [
'diboson-ww',
'diboson-wz',
'diboson-zz',
'zjets',
'singletop-s',
'singletop-t',
'singletop-tW',
'singletop-tbar',
'singletop-tbarW',
'ttbar-ext',
'ttbar',
'wjets1',
'wjets2',
'wjets3',
'wjets4'
]
dataname = 'SingleElectron' if selection == 'ejets' else 'SingleMuon'
samples += [dataname]

filelists = dict( (k,util.file2list('config/cyminiana_samples/{0}.txt'.format(k))) for k in samples)

extralabel = {"ejets":sample_labels['ejets'].label,"mujets":sample_labels['mujets'].label}

# Dictionary for default files
h_hists_dict = dict( (k,{}) for k in samples)
for h in histograms:
    h = h.format(selection)
    for samp in samples:
        h_hists_dict[samp][h] = None

# Load histograms
for samp in samples:
    print "     > Opening data from {0} files".format(samp)

    for file in filelists[samp]:
        f = ROOT.TFile.Open(file)
#        m = f.Get("tree/metadata")
#        m.GetEntry(0)
        sf = 1.0

#        if samp.startswith('ttbar') or samp==dataname : 
#            sf = 1.0
#        else:
#            # Need to normalize the samples properly!
#            samp_metadata = metadata[str(m.primaryDataset)]
#            sf = 1/samp_metadata.sumOfWeights

        for hist in histograms:
            hi = hist.format(selection)

            h_hist = getattr(f,hi)
            h_hist.Scale(sf)
            if h_hists_dict[samp][hi] is None: 
                h_hists_dict[samp][hi] = h_hist.Clone()
                h_hists_dict[samp][hi].SetDirectory(0)
            else: h_hists_dict[samp][hi].Add( h_hist.Clone() )



## Make the plots
for histogram in histograms:

    histogram = histogram.format(selection)

    histogramName = histogram.replace("h_","").replace("_"+selection,"")
    print "  :: Plotting "+histogram

    ## setup histogram
    hist = HepPlotterDataMC()

    xlimits = None
    legendLoc   = 1
    plotLogAxis = False

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
    hist.logplot     = plotLogAxis
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

    ## -- Add the data from each file
    histo = histogram.format(selection)

    # diboson
    h_hist_diboson = h_hists_dict['diboson-ww'][histo]
    h_hist_diboson.Add( h_hists_dict['diboson-wz'][histo] )
    h_hist_diboson.Add( h_hists_dict['diboson-zz'][histo] )

    # zjets
    h_hist_zjets = h_hists_dict['zjets'][histo]

    # singletop
    h_hist_stop  = h_hists_dict['singletop-s'][histo]
    h_hist_stop.Add( h_hists_dict['singletop-t'][histo] )
    h_hist_stop.Add( h_hists_dict['singletop-tbar'][histo] )
    h_hist_stop.Add( h_hists_dict['singletop-tW'][histo] )
    h_hist_stop.Add( h_hists_dict['singletop-tbarW'][histo] )

    # w+jets
    h_hist_wjets = h_hists_dict['wjets1'][histo]
    h_hist_wjets.Add( h_hists_dict['wjets2'][histo] )
    h_hist_wjets.Add( h_hists_dict['wjets3'][histo] )
    h_hist_wjets.Add( h_hists_dict['wjets4'][histo] )

    # ttbar
    h_hist_ttbar = h_hists_dict['ttbar-ext'][histo]
#    h_hist_ttbar.Add( h_hists_dict['ttbar'][histo] )

    h_hist_data  = h_hists_dict[dataname][histo]


    hist.Add(h_hist_diboson, name="diboson", sampleType="background", systematics=None)
    hist.Add(h_hist_zjets, name="zjets",     sampleType="background", systematics=None)
    hist.Add(h_hist_stop,  name="singletop", sampleType="background", systematics=None)
    hist.Add(h_hist_wjets, name="wjets",     sampleType="background", systematics=None)
    hist.Add(h_hist_ttbar, name="ttbar",     sampleType="background", systematics=None)
    hist.Add(h_hist_data,  name="data",      sampleType="data",       systematics=None)


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

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "     ++ Saved plot to "+hist.saveAs+"\n"


## THE END
