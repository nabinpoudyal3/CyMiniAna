"""
Created:        3  September 2016
Last Updated:   9  March     2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Steering script for making simple efficiency plots from TEfficiency objects.
This can be modified or extended by whomever.

To run:
 python python/runEfficiency.py --files <files.txt> --hists <histogramNames.txt> -o <output_path>
"""
import sys
import ROOT
import util
from argparse import ArgumentParser

from Analysis.CyMiniAna.hepPlotter.hepPlotter import HepPlotter
import Analysis.CyMiniAna.hepPlotter.hepPlotterTools as hpt
import Analysis.CyMiniAna.hepPlotter.hepPlotterLabels as hpl


parser = ArgumentParser(description="Efficiency Plotter")

parser.add_argument('-f','--files', action='store',default=None,
                    dest='listOfFiles',
                    help='Name of file that contains root files to plot')
parser.add_argument('--hists', action='store',default=None,
                    dest='listOfHists',
                    help='Name of file that contains histograms to plot')
parser.add_argument('-o','--outpath', action='store',default=None,
                    dest='outpath',
                    help='Directory for storing output plots')
results = parser.parse_args()

outpath    = results.outpath
files      = util.file2list(results.listOfFiles)  # ROOT files to read
histograms = util.file2list(results.listOfEffs)   # TEfficiencies/Histograms to plot

labels    = hpl.variable_labels()
variable  = 'jet_pt'                              # Setup to plot multiple efficiencies/hists as a function of one variable
betterColors = hpt.betterColors()['linecolors']

## Add the data from each file
## Assume the data is structured such that you want to plot
## multiple efficiencies from the same file in one plot
## -> change to your desired structure / plot
##    e.g., to plot efficiencies from two sources (files) on 1 plot: 
##          switch order of file & hist loops
##          To plot multiple kinds of variables on different plots, 
##          you'll need another loop
for file in files:
    f = ROOT.TFile.Open(file)
    filename = file.split("/")[-1].split(".")[0]

    print "  > Opening data from ",filename

    ## setup histogram
    hist = HepPlotter("efficiency",1)

    hist.drawEffDist = True    # draw the physics distribution for efficiency (jet_pt for jet trigger)
    hist.rebin       = 1
    hist.x_label     = labels[variable].label
    hist.y_label     = "Efficiency"
    hist.format      = 'png'       # file format for saving image
    hist.saveAs      = outpath+"eff_"+filename # save figure with name
    hist.CMSlabel       = 'top left'  # 'top left', 'top right'; hack code for something else
    hist.CMSlabelStatus = 'Simulation Internal'  # ('Simulation')+'Internal' || 'Preliminary' 

    #hist.extra_text.Add(text,coords[0,0])

    hist.initialize()

    # loop over variables to put on one plot
    for hi,histogram in enumerate(histograms):

        print "    :: Plotting "+histogram

        h_hist = getattr(f,histogram)       # retrieve the histogram

        # if it is a histogram, just plot black errorbars
        if isinstance(h_hist,ROOT.TH1):
            lineColor = 'k'
        else:
            lineColor = betterColors[hi]

        hist.Add(h_hist,name=histogram,draw='errorbar',label=labels[histogram],
                 linecolor=lineColor)

    plot = hist.execute()
    hist.savefig()
    print "  > Saved plot to "+hist.saveAs+"\n"


## THE END
