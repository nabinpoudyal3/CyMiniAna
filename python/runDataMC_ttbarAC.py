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


parser = ArgumentParser(description="DataMC Plotter")

parser.add_argument('--hists', action='store',dest='listOfHists',
                    default='examples/share/listOfHistsDataMC.txt',
                    help='Name of file that contains histograms to plot')
parser.add_argument('-o','--outpath', action='store',default='./',
                    dest='outpath',
                    help='Directory for storing output plots')
results = parser.parse_args()

listOfHists = results.listOfHists
outpath = results.outpath

histograms = util.file2list(listOfHists)

rebins   = {
"jet_pt":10,\
"ljet_tau32wta":10,\
"deltay":50,\
"beta_ttbar":array('d',[0.,0.3,0.6,1.0]),\
"pt_ttbar":array('d',[0.0,50.,100,200.]),\
"m_ttbar":array('d',[0.0,500,800,1000,1200,1500,2000,2500,4000,5000]),
}

x_labels = hpl.variable_labels()

stop_files  = util.file2list('config/plots/listOfSingleTopFiles.txt')
wjets_files = util.file2list('config/plots/listOfWjetsFiles.txt')
ttbar_files = util.file2list('config/plots/listOfTtbarFiles.txt')
data_files  = util.file2list('config/plots/listOfDataFiles.txt')



# For DataMC plots, one histogram with each sample goes on a single plot
# so I have this structured to loop over the histograms, then the ROOT files
# with lots of histograms, it is probably better to open root files only once,
# load all the histograms, then do the plotting
for histogram in histograms:

    histogramName = histogram.replace("h_","").replace("_eventVars","")
    print "  :: Plotting "+histogram+"\n"

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
    hist.rebin       = rebins[histogramName]  # rebin per histogram
    hist.logplot     = plotLogAxis
    hist.x_label     = x_labels[histogramName].label
    hist.y_label     = "Events"
    hist.y_ratio_label = "Data/MC"
    hist.lumi          = '{0:.1f}'.format( (24799.9+3212.96)*1e-3 )
    hist.CMSlabel    = 'top left'  # 'top left', 'top right'; hack code for something else
    hist.numLegendColumns = 1
    hist.CMSlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
    hist.format           = 'pdf'       # file format for saving image
    hist.saveAs           = outpath+"/datamc_ttbarAC_"+histogramName # save figure with name

    hist.initialize()

    ## -- Add the data from each file
    h_hist_stop  = None
    h_hist_wjets = None
    h_hist_ttbar = None
    h_hist_data  = None


    print "     > Opening data from stop files"
    for file in stop_files:
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)
        h_hist = getattr(f,histogram)
        h_hist.Scale( 1 )
        if h_hist_stop is None: 
            h_hist_stop = h_hist.Clone()
            h_hist_stop.SetDirectory(0)
        else: h_hist_stop.Add( h_hist.Clone() )

    print "     > Opening data from wjets files"
    for file in wjets_files:
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)
        h_hist = getattr(f,histogram)
        h_hist.Scale( 1 )
        if h_hist_wjets is None:
            h_hist_wjets = h_hist.Clone()
            h_hist_wjets.SetDirectory(0)
        else: h_hist_wjets.Add( h_hist.Clone() )

    print "     > Opening data from ttbar files"
    for file in ttbar_files:
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)
        h_hist = getattr(f,histogram)
        h_hist.Scale( 1 )
        if h_hist_ttbar is None: 
            h_hist_ttbar = h_hist.Clone()
            h_hist_ttbar.SetDirectory(0)
        else: h_hist_ttbar.Add( h_hist.Clone() )

    print "     > Opening data from data files"
    for file in data_files:
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)
        h_hist = getattr(f,histogram)
        h_hist.Scale( 1 )
        if h_hist_data is None: 
            h_hist_data = h_hist.Clone()
            h_hist_data.SetDirectory(0)
        else: h_hist_data.Add( h_hist.Clone() )


    hist.Add(h_hist_stop,  name="singletop", sampleType="background", systematics=None)
    hist.Add(h_hist_wjets, name="wjets",     sampleType="background", systematics=None)
    hist.Add(h_hist_ttbar, name="ttbar",     sampleType="background", systematics=None)
    hist.Add(h_hist_data,  name="data",      sampleType="data",       systematics=None)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


## THE END
