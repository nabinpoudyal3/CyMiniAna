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


def calcQCD(rfile,histogram):
    """Calculate the QCD histogram"""
    # 0.5*(yields["G"]/yields["A"] + yields["H"]/yields["B"])*yields["C"]
    h_A = getattr(rfile,histogram.replace("nominal","0b0t_nominal"))
    h_B = getattr(rfile,histogram.replace("nominal","1b0t_nominal"))
    h_C = getattr(rfile,histogram.replace("nominal","2b0t_nominal"))
    h_G = getattr(rfile,histogram.replace("nominal","0b2t_nominal"))
    h_H = getattr(rfile,histogram.replace("nominal","1b2t_nominal"))

    qcd_hist = h_C.Clone()
    qcd_hist.Scale(0.5)

    h_G.Divide(h_A)
    h_H.Divide(h_B)

    h_G.Add(h_H)
    qcd_hist.Multiply(h_G)

    return qcd_hist



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
results = parser.parse_args()

listOfFiles = results.listOfFiles
listOfHists = results.listOfHists
listOfDetectorSysts = results.listOfSysts  # just detector systematics
outpath = results.outpath

files      = open(listOfFiles,"r").readlines()
histograms = open(listOfHists,"r").readlines()
detectorSystematics = open(listOfDetectorSysts,"r").readlines()
fileEnding = "_pre_signal"

# setup for examples choices
x_labels = {
"n_btags":"Number of b-tags (77\%)",\
"n_trackjets":"Number of Track Jets",\
"n_topCandidates":"Number of Large-R Jet Top Candidates",\
"n_ljets":"Number of Large-R Jets",\
"deltaR_t1_t2":r"$\Delta$R(top$_\text{1}$,top$_\text{2}$)",\
"ljet_m_t1":"Leading Large-R Jet Candidate Mass [GeV]",\
"ljet_m_t2":"Sub-leading Large-R Jet Candidate Mass [GeV]",\
"ljet_eta_t1":r"Leading Large-R Jet Candidate $\eta$",\
"ljet_eta_t2":r"Sub-leading Large-R Jet Candidate $\eta$",\
"ljet_pt_t1":r"Leading Large-R Jet Candidate p$_\text{T}$",\
"ljet_pt_t2":r"Sub-leading Large-R Jet Candidate p$_\text{T}$",\
"ljet_tau32wta_t1":r'Leading Large-R Jet Candidate $\tau_{\text{32}}^{\text{wta}}$',\
"ljet_tau32wta_t2":r'Sub-leading Large-R Jet Candidate $\tau_{\text{32}}^{\text{wta}}$',\
"ljet_charge_t1":r"Leading Large-R Jet Candidate Charge",\
"ljet_charge_t2":r"Sub-leading Large-R Jet Candidate Charge",\
"ljet_btagged_tjet_charge_t1":r"Leading Large-R Jet Candidate b-tagged Track Subjets Charge",\
"ljet_btagged_tjet_charge_t2":r"Sub-leading Large-R Jet Candidate b-tagged Track Subjets Charge",\
"ljet_nbtagged_tjet_charge_t1":r"Leading Large-R Jet Candidate Non-b-tagged Track Subjets Charge",\
"ljet_nbtagged_tjet_charge_t2":r"Sub-leading Large-R Jet Candidate Non-b-tagged Track Subjets Charge",\
"ljet_deltaQ_btag_nbtag_t1":r"Leading Large-R Jet Candidate Track Subjets $\Delta$(Q$_\text{b-tag}$,Q$_\text{non-b-tag}$)",\
"ljet_deltaQ_btag_nbtag_t2":r"Sub-leading Large-R Jet Candidate Track Subjets $\Delta$(Q$_\text{b-tag}$,Q$_\text{non-b-tag}$)",\
"ljet_subjet_pt0_charge_t1":r"Leading Large-R Jet Candidate Leading p$_\text{T}$ Track Subjet Charge",\
"ljet_subjet_pt0_charge_t2":r"Sub-leading Large-R Jet Candidate Leading p$_\text{T}$ Track Subjet Charge",\
"ljet_subjet_pt0_btag_charge_t1":r"Leading Large-R Jet Candidate Leading p$_\text{T}$ b-tagged Track Subjet Charge",\
"ljet_subjet_pt0_btag_charge_t2":r"Sub-leading Large-R Jet Candidate Leading p$_\text{T}$ b-tagged Track Subjet Charge",\
"ljet_subjet_pt0_nbtag_charge_t1":r"Leading Large-R Jet Candidate Leading p$_\text{T}$ Non-b-tagged Track Subjet Charge",\
"ljet_subjet_pt0_nbtag_charge_t2":r"Sub-leading Large-R Jet Candidate Leading p$_\text{T}$ Non-b-tagged Track Subjet Charge",\
"deltaR_t_tbar":r"$\Delta$R(top,anti-top)",\
"ljet_m_t":"Top Quark Candidate Mass [GeV]",\
"ljet_m_tbar":"Top Antiquark Candidate Mass [GeV]",\
"ljet_eta_t":r"Top Quark Candidate $\eta$",\
"ljet_eta_tbar":r"Top Antiquark Candidate $\eta$",\
"ljet_pt_t":r"Top Quark Candidate p$_\text{T}$",\
"ljet_pt_tbar":r"Top Antiquark Candidate p$_\text{T}$",\
"ljet_DNN_t":r"Top Quark DNN",\
"ljet_DNN_tbar":r"Top Antiquark DNN",\
"ljet_charge_t":r"Top Quark Charge",\
"ljet_charge_tbar":r"Top Antiquark Charge",\
"ljet_btagged_tjet_charge_t":r"Top Quark b-tagged Track Subjets Charge",\
"ljet_btagged_tjet_charge_tbar":r"Top Antiquark b-tagged Track Subjets Charge",\
"ljet_nbtagged_tjet_charge_t":r"Top Quark Non-b-tagged Track Subjets Charge",\
"ljet_nbtagged_tjet_charge_tbar":r"Top Antiquark Non-b-tagged Track Subjets Charge",\
"ljet_deltaQ_btag_nbtag_t":r"Top Quark Track Subjets $\Delta$(Q$_\text{b-tag}$,Q$_\text{non-b-tag}$)",\
"ljet_deltaQ_btag_nbtag_tbar":r"Top Antiquark Track Subjets $\Delta$(Q$_\text{b-tag}$,Q$_\text{non-b-tag}$)",\
"ljet_subjet_pt0_charge_t":r"Top Quark Leading p$_\text{T}$ Track Subjet Charge",\
"ljet_subjet_pt0_charge_tbar":r"Top Antiquark Leading p$_\text{T}$ Track Subjet Charge",\
"ljet_subjet_pt0_btag_charge_t":r"Top Quark Leading p$_\text{T}$ b-tagged Track Subjet Charge",\
"ljet_subjet_pt0_btag_charge_tbar":r"Top Antiquark Leading p$_\text{T}$ b-tagged Track Subjet Charge",\
"ljet_subjet_pt0_nbtag_charge_t":r"Top Quark Leading p$_\text{T}$ Non-b-tagged Track Subjet Charge",\
"ljet_subjet_pt0_nbtag_charge_tbar":r"Top Antiquark Leading p$_\text{T}$ Non-b-tagged Track Subjet Charge",\
"ljet_tau32wta_t":r'Top Quark Candidate $\tau_{\text{32}}^{\text{wta}}$',\
"ljet_tau32wta_tbar":r'Top Antiquark Candidate $\tau_{\text{32}}^{\text{wta}}$',\
"ljet_t50":r"Large-R Jet Top Tagger (50\%)",\
"ljet_t80":r"Large-R Jet Top Tagger (80\%)",\
"ljet_m":"Large-R Jet Mass [GeV]",\
"ljet_eta":r"Large-R Jet $\eta$",\
"ljet_pt":r"Large-R Jet p$_\text{T}$",\
"ljet_tau32wta":r'Large-R Jet $\tau_{\text{32}}^{\text{wta}}$',\
"deltay":r"$\Delta |\text{y}|$",\
"beta_ttbar":r"$\beta_{z,\text{t}\bar{\text{t}}}$",\
"pt_ttbar":r"p$_{\text{T},\text{t}\bar{\text{t}}}$ [GeV]",\
"m_ttbar":r"m$_{\text{t}\bar{\text{t}}}$ [GeV]"}


rebins   = {
"n_btags":1,\
"n_trackjets":1,\
"n_topCandidates":1,\
"n_ljets":1,\
"deltaR_t1_t2":5,\
"ljet_m_t1":14,\
"ljet_m_t2":14,\
"ljet_eta_t1":4,\
"ljet_eta_t2":4,\
"ljet_pt_t1":140,\
"ljet_pt_t2":140,\
"ljet_tau32wta_t1":10,\
"ljet_tau32wta_t2":10,\
"ljet_charge_t1":10,\
"ljet_charge_t2":10,\
"ljet_btagged_tjet_charge_t1":10,\
"ljet_btagged_tjet_charge_t2":10,\
"ljet_nbtagged_tjet_charge_t1":10,\
"ljet_nbtagged_tjet_charge_t2":10,\
"ljet_deltaQ_btag_nbtag_t1":20,\
"ljet_deltaQ_btag_nbtag_t2":20,\
"ljet_subjet_pt0_charge_t1":10,\
"ljet_subjet_pt0_charge_t2":10,\
"ljet_subjet_pt0_btag_charge_t1":10,\
"ljet_subjet_pt0_btag_charge_t2":10,\
"ljet_subjet_pt0_nbtag_charge_t1":10,\
"ljet_subjet_pt0_nbtag_charge_t2":10,\
"deltaR_t_tbar":5,\
"ljet_m_t":14,\
"ljet_m_tbar":14,\
"ljet_eta_t":4,\
"ljet_eta_tbar":4,\
"ljet_pt_t":140,\
"ljet_pt_tbar":140,\
"ljet_tau32wta_t":10,\
"ljet_tau32wta_tbar":10,\
"ljet_DNN_t":10,\
"ljet_DNN_tbar":10,\
"ljet_charge_t":10,\
"ljet_charge_tbar":10,\
"ljet_btagged_tjet_charge_t":10,\
"ljet_btagged_tjet_charge_tbar":10,\
"ljet_nbtagged_tjet_charge_t":10,\
"ljet_nbtagged_tjet_charge_tbar":10,\
"ljet_deltaQ_btag_nbtag_t":20,\
"ljet_deltaQ_btag_nbtag_tbar":20,\
"ljet_subjet_pt0_charge_t":10,\
"ljet_subjet_pt0_charge_tbar":10,\
"ljet_subjet_pt0_btag_charge_t":10,\
"ljet_subjet_pt0_btag_charge_tbar":10,\
"ljet_subjet_pt0_nbtag_charge_t":10,\
"ljet_subjet_pt0_nbtag_charge_tbar":10,\
"ljet_t50":1,\
"ljet_t80":1,\
"ljet_m":14,\
"ljet_eta":4,\
"ljet_pt":140,\
"ljet_tau32wta":10,\
"deltay":50,\
"beta_ttbar":array('d',[0.,0.3,0.6,1.0]),\
"pt_ttbar":array('d',[0.0,50.,100,200.]),\
"m_ttbar":array('d',[0.0,500,800,1000,1200,1500,2000,2500,4000,5000])}


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

    xlimits = None
    if histogramName.startswith('deltaR'):
        xlimits = (1.0,5.0)
    elif histogramName.startswith('ljet_m'):
        xlimits = (50,300)
    elif histogramName.startswith('ljet_pt'):
        xlimits = (200,3000)
    elif 'tau32' in histogramName:
        xlimits = (0,1.0)
    elif 'deltaQ' in histogramName or 'charge' in histogramName:
        xlimits = (-6.,6.)

    if '_eta_' in histogramName or 'beta_ttbar' in histogramName or 'pt_ttbar' in histogramName:
        hist.ymaxScale = 1.9
    else:
        hist.ymaxScale = None

    legendLoc   = 6 if (histogramName.endswith('50') or histogramName.endswith('80') or histogramName.endswith('DNN_t')) else 1
    plotLogAxis = True if histogramName.startswith('m_ttbar') else False

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
    hist.binning     = 20
    hist.rebin       = rebins[histogramName] # rebin per histogram
    hist.logplot     = plotLogAxis
    hist.x_label     = x_labels[histogramName]
    hist.y_label     = "Events"
    hist.y_ratio_label = "Data/Pred."
    hist.lumi          = '{0:.1f}'.format( (24799.9+3212.96)*1e-3 )
    hist.ATLASlabel    = 'top left'  # 'top left', 'top right'; hack code for something else
    hist.numLegendColumns = 1
    hist.ATLASlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary' 
    hist.format           = 'png'       # file format for saving image
    hist.saveAs           = outpath+"/datamc_newHists_zprime_"+histogramName # save figure with name

    hist.initialize()

    ## -- Add the data from each file
    for file in files:
        file = file.rstrip("\n")
        f = ROOT.TFile.Open(file)

        print "     > Opening data from ",file

        h_name       = file.split("/")[-1].split(".root")[0].replace(fileEnding,"")  # retrieve the name ('ttbar','wjets',etc.)
        h_sampleType = hpt.getSampleType(h_name)  # retrieve the sample type ('background','data','signal')

        h_hist = None
        if h_sampleType=='background' and h_name=='multijet':
            h_hist = getattr(f,histogram.replace("nominal","signal_nominal"))  #calcQCD(f,histogram)
        else:
            h_hist = getattr(f,"h_"+histogram.replace("nominal","signal_nominal"))  # retrieve the histogram


        if h_sampleType=='background' and h_name!='multijet':
            h_hist.Scale( (24799.9+3212.96) / 14688.86 )


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
        if "DNN_t_" in histogram and h_name=='ttbar_had':
            h_hist_t_top   = getattr(f,"h_ljet_DNN_t_truth_top_signal_nominal")
            h_hist_t_tbar  = getattr(f,"h_ljet_DNN_t_truth_tbar_signal_nominal")
            h_hist_t_other = getattr(f,"h_ljet_DNN_t_other_signal_nominal")
            hist.Add(h_hist_t_other,name=h_name+"_other",sampleType=h_sampleType,file=file,systematics=totsyst)
            hist.Add(h_hist_t_tbar, name=h_name+"_tbar" ,sampleType=h_sampleType,file=file,systematics=totsyst)
            hist.Add(h_hist_t_top,  name=h_name+"_top"  ,sampleType=h_sampleType,file=file,systematics=totsyst)
        elif "DNN_tbar_" in histogram and h_name=='ttbar_had':
            h_hist_tbar_top   = getattr(f,"h_ljet_DNN_tbar_truth_top_signal_nominal")
            h_hist_tbar_tbar  = getattr(f,"h_ljet_DNN_tbar_truth_tbar_signal_nominal")
            h_hist_tbar_other = getattr(f,"h_ljet_DNN_tbar_other_signal_nominal")
            hist.Add(h_hist_tbar_other,name=h_name+"_other",sampleType=h_sampleType,file=file,systematics=totsyst)
            hist.Add(h_hist_tbar_top,  name=h_name+"_top"  ,sampleType=h_sampleType,file=file,systematics=totsyst)
            hist.Add(h_hist_tbar_tbar, name=h_name+"_tbar" ,sampleType=h_sampleType,file=file,systematics=totsyst)
        else:
            hist.Add(h_hist,name=h_name,sampleType=h_sampleType,file=file,systematics=totsyst)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


## THE END
