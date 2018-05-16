"""
Created:         1 September 2016
Last Updated:   21 September 2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109
-----
Steering script for making Data/MC plots.
Primarily want to do this from histograms (faster to make those in C++ & 
we feed histograms into the limit-setting framework).

This can be modified or extended by whomever.

To run:
python python/runDataMC.py --files <files.txt> --hists <histogramNames.txt> -o <output_path>
"""
import sys
import ROOT
from array import array
from argparse import ArgumentParser

from hepPlotterDataMC import HepPlotterDataMC
from hepPlotterSystematics import HepPlotterSystematics
import hepPlotterTools as hpt
import hepPlotterLabels as hpl


parser = ArgumentParser(description="DataMC Plotter")

parser.add_argument('-f','--files', action='store',dest='listOfFiles',
                    default='examples/share/listOfFilesDataMC.txt',
                    help='Name of file that contains root files from which to plot')
parser.add_argument('--hists', action='store',dest='listOfHists',
                    default='examples/share/listOfHistsDataMC.txt',
                    help='Name of file that contains histograms to plot')
parser.add_argument('--systs', action='store',dest='listOfSysts',
                    default='share/listOfSytsDataMC.txt',
                    help='Name of file that contains detector systematics to plot')
parser.add_argument('-o','--outpath', action='store',default='./',
                    dest='outpath',
                    help='Directory for storing output plots')
parser.add_argument('-s','--selection', action='store',default='',
                    dest='selection',
                    help='Selection used in CyMiniAna to make histograms')
results = parser.parse_args()

listOfFiles = results.listOfFiles
listOfHists = results.listOfHists
listOfDetectorSysts = results.listOfSysts  # just detector systematics
outpath = results.outpath

files      = open(listOfFiles,"r").readlines()
histograms = open(listOfHists,"r").readlines()
detectorSystematics = open(listOfDetectorSysts,"r").readlines()

# setup for examples choices
if not results.selection:
    selection = util.read_config("config/cmaConfig.txt")["selection"]    # histograms are named based on selection
else:
    selection = results.selection.split(",")
variables = hpl.variable_labels()   # label and binning


# For DataMC plots, one histogram with each sample goes on a single plot
# so I have this structured to loop over the histograms, then the ROOT files
# with lots of histograms, it is probably better to open root files only once,
# load all the histograms, then do the plotting, but oh well
for histogram in histograms:

    histogram = histogram.rstrip('\n')

    # CyMiniAna saves histograms with extra text in the name, remove that here
    # "jet_pt" -> "h_jet_pt_SELECTION"
    histogramName = histogram.replace("h_","")
    for sel in selection:
        histogramName = histogramName.replace("_"+sel,"")
    histogramName = histogramName.replace( "_"+"-".join(selection),"")
    print "  :: Plotting "+histogram+"\n"

    ## setup histogram
    hist = HepPlotterDataMC()

    hist.drawStatUncertainty = True      # draw stat uncertainty separately
    hist.drawSystUncertainty = False     # draw syst uncertainty separately
    hist.drawStatSystUncertainty = False # draw stat+syst uncertainty
    hist.legendLoc   = 1
    hist.stackSignal = False
    hist.blind       = False
    hist.xlim        = None
    hist.ratio_plot  = True        # plot a ratio of things [Data/MC]
    hist.ratio_type  = "ratio"     # "ratio"
    hist.stacked     = True        # stack backgrounds
    hist.binning     = None
    hist.rebin       = variables[histogramName].binning # rebin per histogram
    hist.logplot     = {"y":False,"x":False,"data":False}
    hist.x_label     = variables[histogramName].label
    hist.y_label     = "Events"
    hist.lumi        = 'XY.Z'
    hist.CMSlabel    = 'top left'  # 'top left', 'top right'; hack code for something else
    hist.CMSlabelStatus   = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
    hist.numLegendColumns = 1
    hist.y_ratio_label    = "Data/Pred."
    hist.format           = 'pdf'       # file format for saving image
    hist.saveAs           = outpath+"/datamc_newHists_zprime_"+histogramName # save figure with name

    hist.initialize()

    ## -- Add the data from each file
    for file in files:
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)

        print "     > Opening data from ",file

        h_name       = file.split("/")[-1].split(".root")[0]  # retrieve the name ('ttbar','wjets',etc.)
        h_sampleType = hpt.getSampleType(h_name)  # retrieve the sample type ('background','data','signal')

        h_hist = None
        if h_sampleType=='background' and h_name=='multijet':
            h_hist = getattr(f,histogram.replace("nominal","signal_nominal"))  #calcQCD(f,histogram)
        else:
            h_hist = getattr(f,"h_"+histogram.replace("nominal","signal_nominal"))  # retrieve the histogram


        ## -- Get all the systematics for this histogram
        totsyst = None

        if (hist.drawSystUncertainty or hist.drawStatSystUncertainty) and sampleType=='background':
            hpSysts = HepPlotterSystematics()
            hpSysts.sampleName = h_name
            hpSysts.variable = histogramName
            hpSysts.outpath  = "./"                  # directory for individual syst plots
            hpSysts.rebin    = rebins[histogramName] # rebin histogram
            hpSysts.plotSystematics = True           # make individual plots of systematics

            hpSysts.initialize(h_hist)               # setup some systematic variables

            # loop over detector systematics (e.g., JES,JER,btagging,etc.)
            for detSyst in detectorSystematics:
                detSyst = detSyst.rstrip('\n')
                detSyst = detSyst.split(",")    # for systematics that are "up,down"; 
                                                # returns 1 item list if only one systematics
                h_systs = [getattr( f,d ) for d in detSyst] # stored in same file as nominal
                hpSysts.execute(histogram=h_systs,systematic=detSyst)

            # loop over other uncertainties here
            # may need to specify different root file (modeling, XSection, theory, etc.)

            totsyst = hpSysts.getTotalSystematicUncertainty()

        hist.Add(h_hist,name=h_name,sampleType=h_sampleType,file=file,systematics=totsyst)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


## THE END
