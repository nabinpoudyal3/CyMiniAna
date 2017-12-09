"""
Created:         1 August 2017
Last Updated:    5 August 2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----
Example for making Data/MC plots with systematic uncertainties

To see options, use:
$ python python/datamcSysts.py --help
"""
import sys
import ROOT
import collections
from array import array
from argparse import ArgumentParser

from hepPlotter import HepPlotter
from hepPlotterDataMC import HepPlotterDataMC
from hepPlotterSystematics import HepPlotterSystematics
import hepPlotterTools as hpt
import hepPlotterLabels as hpl


parser = ArgumentParser(description="DataMC Plotter")

parser.add_argument('-f','--files', action='store',dest='listOfFiles',
                    default='examples/share/listOfFilesDataMC.txt',
                    help='Name of file that contains root files from which to plot')
parser.add_argument('--systs', action='store',dest='listOfSysts',
                    default='share/listOfSystsDataMC.txt',
                    help='Name of file that contains detector systematics to plot')
parser.add_argument('-o','--outpath', action='store',default='plots/',
                    dest='outpath',
                    help='Directory for storing output plots')
parser.add_argument('-dnn',action='store',default='AustralianShepherd',
                    dest='dnn',
                    help='Neural network used for top ID')
results = parser.parse_args()

listOfFiles = results.listOfFiles
listOfDetectorSysts = results.listOfSysts  # just detector systematics
outpath = results.outpath
dnn     = results.dnn

detectorSystematics = open(listOfDetectorSysts,"r").readlines()  # list of detector-related systematics
files_all    = open(listOfFiles,"r").readlines()                 # All files to plot
files_others = open("share/listOfFiles2Plot_NNs_mergedOthers.txt","r").readlines() # merged 'others' backgrounds for >=2 b-tags

# remove parts of filename to easily identify the file 
fileEnding   = "_TtbarRes_{0}_ljets_none_output" 
fileEnding2  = fileEnding.replace("_ljets_","_output_")

DNNCONFIG    = {'AustralianShepherd':'Australian Shepherd',
                'GoldenRetriever':'Golden Retriever'}

# setup for examples choices
x_labels = {
"ljet_DNN_htop":r"Large-R Jet Candidate DNN"}

rebins   = {
"ljet_DNN_htop":10}


region_names = collections.OrderedDict()
region_names["1bep"] = r"1 b-tag, e$^{+}$"
region_names["1bem"] = r"1 b-tag, e$^{-}$"
region_names["1bmp"] = r"1 b-tag, $\mu^{+}$"
region_names["1bmm"] = r"1 b-tag, $\mu^{-}$"
region_names["2bep"] = r"$\ge$2 b-tags, e$^{+}$"
region_names["2bem"] = r"$\ge$2 b-tags, e$^{-}$"
region_names["2bmp"] = r"$\ge$2 b-tags, $\mu^{+}$"
region_names["2bmm"] = r"$\ge$2 b-tags, $\mu^{-}$"
# region_names['1p']  = r"$\ge$1 b-tag, $\ell^{+}$"
# region_names['1m']  = r"$\ge$1 b-tag, $\ell^{-}$"
regions = ["2bep"] #region_names.keys()

print regions

histogram     = "ljet_DNN_htop_{0}_nominal"
histogramName = histogram.split("_{0}")[0]


print "  :: Plotting "+histogram+"\n"

xlimits     = None
plotLogAxis = False
xcoord      = 0.03
halign      = "left"
legendLoc   = 1 
atlasLabel  = 'top left'

# Make one plot per region for each distribution
for region in regions:

    print "  > Region ",region
    saveAs = "{0}/datamc_detsysts_{1}_{2}".format(outpath,histogramName,region)

    if region.endswith("m"):
        xcoord = 0.97
        halign = "right"
        legendLoc  = 2
        atlasLabel = 'top right'
    else:
        xcoord = 0.03
        halign = "left"
        legendLoc  = 1
        atlasLabel = 'top left'

    ## setup histogram
    hist = HepPlotterDataMC()

    hist.drawStatUncertainty = False    # draw stat uncertainty separately
    hist.drawSystUncertainty = False    # draw syst uncertainty separately
    hist.drawStatSystUncertainty = True # draw stat+syst uncertainty
    hist.legendLoc   = legendLoc
    hist.stackSignal = False
    hist.blind       = False
    hist.xlim        = xlimits
    hist.ratio_plot  = True        # plot a ratio of things [Data/MC]
    hist.ratio_type  = "ratio"     # "ratio"
    hist.stacked     = True        # stack plots
    hist.rebin       = rebins[histogramName] # rebin per histogram
    hist.logplot     = plotLogAxis
    hist.x_label     = x_labels[histogramName]
    hist.y_label     = "Events"
    hist.y_ratio_label = "Data/Pred."
    hist.lumi          = '36.1'
    hist.ATLASlabel    = atlasLabel  # 'top left', 'top right'; hack code for something else
    hist.numLegendColumns = 1
    hist.ATLASlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
    hist.format           = 'png'       # file format for saving image
    hist.saveAs           = saveAs
    hist.extra_text.Add( region_names[region], coords=[xcoord,0.80], ha=halign )
    hist.extra_text.Add( DNNCONFIG[dnn], coords=[xcoord,0.74], ha=halign )

    hist.initialize()

    ## -- Add the data from each file
    if region.startswith("2"): files = files_others
    else: files = files_all

    for file in files:
        file = file.format(dnn,region).rstrip("\n")
        f = ROOT.TFile.Open(file)

        print "     > Opening data from ",file

        # retrieve the name ('ttbar','wjets',etc.)
        h_name       = file.split("/")[-1].split(".root")[0].replace(fileEnding.format(region),"").replace(fileEnding2.format(region),"")
        h_sampleType = hpt.getSampleType(h_name)  # retrieve the sample type ('background','data','signal')

        h_histName = "h_"+histogram.format("none")
        h_hist = None
        h_hist = getattr(f,h_histName)  # retrieve the histogram

        ## -- Get all the systematics for this histogram
        totsyst = None

        if (hist.drawSystUncertainty or hist.drawStatSystUncertainty) and h_sampleType=='background':
            print "       -- systematics "
            hpSysts = HepPlotterSystematics()
            hpSysts.binning  = 1
            hpSysts.sampleName = h_name
            hpSysts.variable = histogramName
            hpSysts.outpath  = "{0}/systs/{1}/".format(outpath,region) # directory for individual syst plots
            hpSysts.rebin    = rebins[histogramName]    # rebin histogram
            hpSysts.plotSystematics = True              # make individual plots of systematics

            hpSysts.initialize(h_hist.Clone())          # setup some systematic variables

            # loop over detector systematics (e.g., JES,JER,btagging,etc.)
            for detSyst in detectorSystematics:
                # h_ljet_DNN_htop_none_nominal
                # h_ljet_DNN_htop_none_LARGERJET_Medium_JET_Comb_Tracking_Kin__1up

                detSyst = detSyst.rstrip('\n')
                detSyst = detSyst.split(",")    # for systematics that are "up,down"; 
                                                # returns 1 item list if only one systematics
                h_systs = [getattr( f,h_histName.replace("nominal",d) ) for d in detSyst] # stored in same file as nominal
                hpSysts.execute(histogram=h_systs,systematic=detSyst)

            # loop over other uncertainties here
            # may need to specify different root file (modeling, XSection, theory, etc.)

            totsyst = hpSysts.getTotalSystematicUncertainty()

        ## -- Add the histogram to the figure
        hist.Add(h_hist,name=h_name,sampleType=h_sampleType,file=file,systematics=totsyst)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "   * Saved plot to "+saveAs+" *\n"

## THE END

