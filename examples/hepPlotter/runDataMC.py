"""
Created:         1 September 2016
Last Updated:   21 September 2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
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
from argparse import ArgumentParser

from Analysis.CyMiniAna.hepPlotter.hepPlotterDataMC import HepPlotterDataMC
from Analysis.CyMiniAna.hepPlotter.hepPlotterSystematics import HepPlotterSystematics
import Analysis.CyMiniAna.hepPlotter.hepPlotterTools as hpt
import Analysis.CyMiniAna.hepPlotter.hepPlotterLabels as hpl


parser = ArgumentParser(description="DataMC Plotter")

parser.add_argument('-f','--files', action='store',dest='listOfFiles',
                    default='examples/share/listOfFilesDataMC.txt',
                    help='Name of file that contains root files from which to plot')
parser.add_argument('--hists', action='store',dest='listOfHists',
                    default='examples/share/listOfHistsDataMC.txt',
                    help='Name of file that contains histograms to plot')
parser.add_argument('--systs', action='store',dest='listOfSysts',
                    default='examples/share/listOfSystsDataMC.txt',
                    help='Name of file that contains detector systematics to plot')
parser.add_argument('-o','--outpath', action='store',default='./',
                    dest='outpath',
                    help='Directory for storing output plots')
results = parser.parse_args()

listOfFiles = results.listOfFiles
listOfHists = results.listOfHists
listOfDetectorSysts = results.listOfSysts  # just detector systematics
outpath = results.outpath

files      = open(listOfFiles,"r").readlines()
histograms = open(listOfHists,"r").readlines()
detectorSystematics = open(listOfDetectorSysts,"r").readlines()


# setup for examples choices
x_labels = hpl.variable_labels()
rebins   = {"leptonicT_m":40,"ST":200}

# For DataMC plots, one histogram with each sample goes on a single plot
# so I have this structured to loop over the histograms, then the ROOT files
# with lots of histograms, it is probably better to open root files only once,
# load all the histograms, then do the plotting
for histogram in histograms:

    histogram = histogram.rstrip('\n')
    histogramName = histogram.replace("_nominal","").replace("h_","")
    print "  :: Plotting "+histogram+"\n"

    ## setup histogram
    hist = HepPlotterDataMC()

    hist.drawStatUncertainty = True      # draw stat uncertainty separately
    hist.drawSystUncertainty = False     # draw syst uncertainty separately
    hist.drawStatSystUncertainty = True  # draw stat+syst uncertainty
    hist.stackSignal = False
    hist.blind       = False
    hist.ratio_plot  = True        # plot a ratio of things [Data/MC]
    hist.ratio_type  = "ratio"     # "ratio"
    hist.stacked     = True        # stack plots
    hist.binning     = 20
    hist.rebin       = rebins[histogramName] # rebin per histogram
    hist.logplot     = False       # plot on log scale
    hist.x_label     = x_labels[histogramName].label
    hist.y_label     = "Events"
    hist.y_ratio_label = "Data/Pred."
    hist.lumi         = '1'           # in /fb
    hist.CMSlabel       = 'top left'  # 'top left', 'top right'; hack code for something else
    hist.CMSlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
    hist.numLegendColumns = 1
    hist.format           = 'png'       # file format for saving image
    hist.saveAs           = outpath+"datamc_"+histogramName # save figure with name

    hist.extra_text.Add("Extra text",coords=[0.5,0.5]) # see hepPlotterLabels.py for exact use of extra_text (PlotText() objects)

    hist.initialize()

    ## -- Add the data from each file
    for file in files:
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)

        print "     > Opening data from ",file

        h_hist       = getattr(f,histogram)       # retrieve the histogram
        h_name       = hpt.getName(file)          # retrieve the name ('ttbar','wjets',etc.)
        h_sampleType = hpt.getSampleType(h_name)  # retrieve the sample type ('background','data','signal')

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

        ## -- Add the histogram to the figure
        hist.Add(h_hist,name=h_name,sampleType=h_sampleType,file=file,systematics=totsyst)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


## THE END
