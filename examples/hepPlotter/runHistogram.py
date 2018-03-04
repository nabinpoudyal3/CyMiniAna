"""
Created:         3 September 2016
Last Updated:    3 September 2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Steering script for making simple histograms.
Primarily want to do this from root histograms (faster to make those in c++)

This can be modified or extended by whomever.

To run:
python python/runHistogram.py --files <files.txt> --hists <histogramNames.txt> -o <output_path>
"""
import sys
import ROOT
from argparse import ArgumentParser

from hepPlotter import HepPlotter
import hepPlotterTools as hpt
import hepPlotterLabels as hpl

parser = ArgumentParser(description="Histogram Plotter")

parser.add_argument('-f','--files', action='store',
                    dest='listOfFiles',default='examples/share/listOfFiles.txt',
                    help='Name of file that contains root files to plot')
parser.add_argument('--hists', action='store',
                    dest='listOfHists',default='examples/share/listOfHists.txt',
                    help='Name of file that contains histograms to plot')
parser.add_argument('-o','--outpath', action='store',
                    dest='outpath',default='examples/hepPlotter/',
                    help='Directory for storing output plots')
results = parser.parse_args()

listOfFiles = results.listOfFiles
listOfHists = results.listOfHists
outpath     = results.outpath

files      = open(listOfFiles,"r").readlines()
histograms = open(listOfHists,"r").readlines()

betterColors = hpt.betterColors()['linecolors']

x_labels = hpl.text_dicts()['variables'] # labels based on histogram plotted (e.g., 'Jet pT')
labels   = hpl.text_dicts()['samples']   # labels based on filetype being plotted (e.g., 'ttbar')

numberOfHists = 0
# Access data -- assumes you are plotting histograms from multiple sources in one figure
for hi,histogram in enumerate(histograms):

    histogram     = histogram.strip('\n')
    histogramName = histogram.split("h_")[1].split("_pre")[0]
    print "  :: Plotting "+histogram+"\n"

    ## setup histogram
    hist = HepPlotter("histogram",1)

    hist.ratio_plot  = False       # plot a ratio of things [Data/MC]
    hist.ratio_type  = "ratio"     # "ratio"
    hist.stacked     = True        # stack plots
    hist.rebin       = 10
    hist.logplot     = False       # plot on log scale
    hist.x_label     = x_labels[histogramName]["label"]
    hist.y_label     = "Events"
    hist.y_ratio_label = ""
    hist.lumi          = 'XY.Z'   # in /fb
    hist.format        = 'png'       # file format for saving image
    hist.saveAs        = outpath+"/hist_"+histogram # save figure with name
    hist.CMSlabel       = 'top left'  # 'top left', 'top right'; hack code for something else
    hist.CMSlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
#    hist.extra_text.Add("text here",coords=[x,y]) # see hepPlotter for exact use of extra_text (PlotText() objects)

    hist.initialize()

    ## Add the data from each file
    for fi,file in enumerate(files):
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)
        filename = file.split("/")[-1].split(".")[0]

        print "     > Opening data from ",file

        h_hist       = getattr(f,histogram)       # retrieve the histogram

        hist.Add(h_hist,name=filename+"_"+histogram,linecolor=betterColors[numberOfHists],
                 draw='step',label=labels[filename]['label'])
        numberOfHists+=1

    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


#### 2D Histogram example

histograms = ["h_jets_n_v_HT_pre_nominal"]

for hi,histogram in enumerate(histograms):

    histogram     = histogram.strip('\n')
    histogramName = histogram.split("h_")[1].split("_pre")[0]
    print "  :: Plotting "+histogram+"\n"
    ## setup histogram
    hist = HepPlotter("histogram",2)

    hist.ratio_plot  = False       # plot a ratio of things [Data/MC]
    hist.ratio_type  = "ratio"     # "ratio"
    hist.stacked     = False       # stack plots
    hist.rebin       = 1
    hist.logplot     = False       # plot on log scale
    hist.x_label     = x_labels[histogramName.split("_v_")[0]]["label"]
    hist.y_label     = x_labels[histogramName.split("_v_")[1].split("_")[0]]["label"]
    hist.y_ratio_label = ""
    hist.lumi          = 'XY.Z'   # in /fb
    hist.format        = 'png'       # file format for saving image
    hist.saveAs        = outpath+"/hist_"+histogram # save figure with name
    hist.colormap      = 'default'
    hist.CMSlabel       = 'outer'     # 'outer' for 2D plots
    hist.CMSlabelStatus = 'Simulation Internal'  # ('Simulation')+'Internal' || 'Preliminary'

    hist.initialize()

    ## Add the data from each file
    for fi,file in enumerate(files):
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)
        filename = file.split("/")[-1].split(".")[0]

        print "     > Opening data from ",file

        h_hist       = getattr(f,histogram)       # retrieve the histogram

        hist.Add(h_hist,name=filename+"_"+histogram,linecolor=betterColors[numberOfHists],
                 draw='step',label=labels[filename]['label'])
        numberOfHists+=1

    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


## THE END

