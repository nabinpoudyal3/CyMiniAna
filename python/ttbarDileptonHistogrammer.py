"""
Created:         7 September 2017
Last Updated:    7 September 2017

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
from array import array
from argparse import ArgumentParser

from hepPlotter import HepPlotter
import hepPlotterTools as hpt
import hepPlotterLabels as hpl


def getFileName(file):
    return file.split("/")[-1].split(".root")[0]



parser = ArgumentParser(description="Histogram Plotter")

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
outpath = results.outpath

histograms = open(results.listOfHists,"r").readlines()

files      = open(results.listOfFiles,"r").readlines()
files      = [i.rstrip('\n') for i in files]
signal_files = [getFileName(i) for i in files if '/GluGlu' in i]
bckg_files   = [getFileName(i) for i in files if '/TT' in i]

betterColors = hpt.betterColors()['linecolors']

rebins   = {
'h_amwt_t':array('d',[-1,-0.5,0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,1]),
'h_amwt_ES_t':array('d',[-1,-0.5,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
}
xlims    = {"h_amwt_t":(-1,0.2),"h_amwt_ES_t":(-1,1)}
x_labels = {"h_amwt_t":"AMWT Weight","h_amwt_ES_t":"Event Shape Weight"} # labels based on histogram plotted (e.g., 'Jet pT')

labels   = {
"GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_Summer16MiniAODv2_v5.0.1+80X_HHAnalysis_2017-03-01.v0_histos_none":r"m$_\text{G}$ = 300 GeV",
"GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_Summer16MiniAODv2_v5.0.1+80X_HHAnalysis_2017-03-01.v0_histos_none":r"m$_\text{G}$ = 500 GeV",
"GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_Summer16MiniAODv2_v5.0.1+80X_HHAnalysis_2017-03-01.v0_histos_none":r"m$_\text{G}$ = 700 GeV",
"GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_Summer16MiniAODv2_v5.0.1+80X_HHAnalysis_2017-03-01.v0_histos_none":r"m$_\text{G}$ = 900 GeV",
"TTTo2L2Nu_13TeV-powheg_Summer16MiniAODv2_v5.0.1+80X_HHAnalysis_2017-03-01.v0_histos_none":r't$\bar{\text{t}}$'}
           # labels based on filetype being plotted (e.g., 'ttbar')

# Access data -- assumes you are plotting histograms from multiple sources in one figure
for hi,histogram in enumerate(histograms):

    histogram = histogram.strip('\n')
    print "  :: Plotting "+histogram+"\n"

    ## setup histogram
    hist = HepPlotter("histogram",1)

    hist.ratio_plot  = True        # plot a ratio of things [Data/MC]
    hist.ratio_type  = "significance"     # "ratio"
    hist.normed      = True
    hist.stacked     = False       # stack plots
    hist.rebin       = rebins[histogram]
    hist.logplot     = False       # plot on log scale
    hist.xlim        = xlims[histogram]
    hist.x_label     = x_labels[histogram]
    hist.y_label     = "Arbitrary Units"
    hist.y_ratio_label  = r"S/$\sqrt{\text{B}}$"
    hist.lumi           = ''   # in /fb
    hist.format         = 'png'       # file format for saving image
    hist.saveAs         = outpath+"hist_"+histogram # save figure with name
    hist.CMSlabel       = 'top right'  # 'top left', 'top right'; hack code for something else
    hist.CMSlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
    hist.numLegendColumns = 1
    hist.legendLoc        = 2

    hist.extra_text.Add(r"G$\rightarrow$H(b$\bar{\text{b}}$)H(WW$^\text{*}$)",coords=[.97,.83],ha='right')

    hist.initialize()

    numberOfHists = 0
    ## Add the data from each file
    for fi,file in enumerate(files):
        f = ROOT.TFile.Open(file)
        filename = getFileName(file)

        print "     > Opening data from ",file

        h_hist = getattr(f,histogram)       # retrieve the histogram

        if filename.startswith("TT"): 
            linecolor = 'k'
            ratio_num = False
            ratio_den = True
            ratio_partners = [i+"_"+histogram for i in signal_files]
        else:
            linecolor = betterColors[numberOfHists*2]
            numberOfHists+=1
            ratio_num = True
            ratio_den = False
            ratio_partners = bckg_files[0]+"_"+histogram

        hist.Add(h_hist,name=filename+"_"+histogram,linecolor=linecolor,
                 draw='step',label=labels[filename],ratio_den=ratio_den,ratio_num=ratio_num,
                 ratio_partner=ratio_partners)

    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


## THE END
