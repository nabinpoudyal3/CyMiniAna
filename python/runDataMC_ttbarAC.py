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

selection = 'mujets'
xsections = {
'ttbar':831.76,
'ttbar-ext':831.76,
'wjets1':676.3,
'wjets2':23.94,
'wjets3':3.031,
'wjets4':0.4524,
'singletop-s':3.36,
'singletop-t':136.02,
'singletop-tW':80.95,
'singletop-tbarW':35.6,
'singletop-tbar':35.6,
}

nevents = {
'ttbar':77229341, #92925926.0,
'ttbar-ext':77709140,
'wjets1':10089661.0,
'wjets2':1001250.0,
'wjets3':248178.0,
'wjets4':304306.0,
'singletop-s':2989199.0,
'singletop-t':5993676.0,
'singletop-tW':998276.0,
'singletop-tbarW':992024.0,
'singletop-tbar':3928063.0,
}
LUMI = 35870/36074.56



def getHist(histname,file,h_hists,weight):
    """Get the histogram from the ROOT file and store in a dictionary"""
#    if 'OTHER' in histname:
#        h_hist = getattr(file,histname.replace("OTHER","BONLY"))
#        h_hist.Add( getattr(file,histname.replace("OTHER","QONLY")) )
#        h_hist.Add( getattr(file,histname.replace("OTHER","NONE")) )
#    else:
#        h_hist = getattr(file,histname)
    h_hist = getattr(file,histname)
    h_hist.Scale( weight )
    if h_hists[histname] is None: 
        h_hists[histname] = h_hist.Clone()
        h_hists[histname].SetDirectory(0)
    else: h_hists[histname].Add( h_hist.Clone() )

    return h_hists


def histfromfile(files,name,hists,isMC=True):
    """load histogram from files"""
    for file in files:
        f = ROOT.TFile.Open(file)
        sf = xsections[name]*LUMI/nevents[name] if isMC else 1.
        for hist in histograms:
            hi = hist.format(selection)
            hists = getHist( hi,f,hists,sf )
    return hists


parser = ArgumentParser(description="DataMC Plotter")

parser.add_argument('--hists', action='store',dest='listOfHists',
                    default='config/listOfHists.txt',
                    help='Name of file that contains histograms to plot')
parser.add_argument('-o','--outpath', action='store',default='plots/29April2018/datamc/',
                    dest='outpath',
                    help='Directory for storing output plots')
results = parser.parse_args()

listOfHists = results.listOfHists
outpath     = results.outpath
histograms  = util.file2list(listOfHists)

x_labels = hpl.variable_labels()
sample_labels = hpl.sample_labels()

stop_s_files     = util.file2list('config/files2plot/singletop-s.txt')
stop_t_files     = util.file2list('config/files2plot/singletop-t.txt')
stop_tbar_files  = util.file2list('config/files2plot/singletop-tbar.txt')
stop_tW_files    = util.file2list('config/files2plot/singletop-tW.txt')
stop_tbarW_files = util.file2list('config/files2plot/singletop-tbarW.txt')
wjets1_files = util.file2list('config/files2plot/wjets1.txt') 
wjets2_files = util.file2list('config/files2plot/wjets2.txt')
wjets3_files = util.file2list('config/files2plot/wjets3.txt')
wjets4_files = util.file2list('config/files2plot/wjets4.txt')
ttbar_files     = util.file2list('config/files2plot/ttbar.txt')
ttbar_ext_files = util.file2list('config/files2plot/ttbar-ext.txt')
data_files  = util.file2list('config/files2plot/SingleElectron.txt') + util.file2list('config/files2plot/SingleMuon.txt')

extralabel = {"ejets":sample_labels['ejets'].label,"mujets":sample_labels['mujets'].label}

contain = [
'NONE',
'QONLY',
'BONLY',
'BQ',
'W',
'FULL',
]

h_hists_stop = {}
h_hists_wjets = {}
h_hists_ttbar = {}
h_hists_data = {}

for h in histograms:
    h = h.format(selection)
    h_hists_stop[h] = None
    h_hists_wjets[h] = None
    h_hists_data[h] = None
#    if 'ljet' in h:
#        h_hists_ttbar[h.replace("_"+selection,"_FULL_{0}".format(selection))]  = None
#        h_hists_ttbar[h.replace("_"+selection,"_W_{0}".format(selection))]     = None
#        h_hists_ttbar[h.replace("_"+selection,"_BQ_{0}".format(selection))]    = None
#        h_hists_ttbar[h.replace("_"+selection,"_OTHER_{0}".format(selection))] = None
#    else:
#        h_hists_ttbar[h] = None
    h_hists_ttbar[h] = None



print "     > Opening data from stop files"
h_hists_stop = histfromfile(stop_s_files,'singletop-s',h_hists_stop)
h_hists_stop = histfromfile(stop_s_files,'singletop-t',h_hists_stop)
h_hists_stop = histfromfile(stop_s_files,'singletop-tbar',h_hists_stop)
h_hists_stop = histfromfile(stop_s_files,'singletop-tW',h_hists_stop)
h_hists_stop = histfromfile(stop_s_files,'singletop-tbarW',h_hists_stop)

print "     > Opening data from wjets files"
h_hists_wjets = histfromfile(wjets1_files,'wjets1',h_hists_wjets)
h_hists_wjets = histfromfile(wjets2_files,'wjets2',h_hists_wjets)
h_hists_wjets = histfromfile(wjets3_files,'wjets3',h_hists_wjets)
h_hists_wjets = histfromfile(wjets4_files,'wjets4',h_hists_wjets)

print "     > Opening data from ttbar files"
h_hists_ttbar = histfromfile(ttbar_files,'ttbar',h_hists_ttbar)
h_hists_ttbar = histfromfile(ttbar_ext_files,'ttbar-ext',h_hists_ttbar)

print "     > Opening data from data files"
h_hists_data = histfromfile(data_files,'data',h_hists_data,isMC=False)





for histogram in histograms:

    histogram = histogram.format(selection)

    histogramName = histogram.replace("h_","").replace("_"+selection,"")
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
    hist.saveAs           = outpath+"/datamc_ttbarAC_contain-log_{0}_{1}".format(selection,histogramName) # save figure with name

    hist.extra_text.Add(extralabel[selection],coords=[0.03,0.90])

    hist.initialize()

    ## -- Add the data from each file
    histo = histogram.format(selection)
    h_hist_stop  = h_hists_stop[histo]
    h_hist_wjets = h_hists_wjets[histo]
    h_hist_data  = h_hists_data[histo]

    hist.Add(h_hist_stop,  name="singletop", sampleType="background", systematics=None)
    hist.Add(h_hist_wjets, name="wjets",     sampleType="background", systematics=None)


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

    h_hist_ttbar = h_hists_ttbar[histo]
    hist.Add(h_hist_ttbar, name="ttbar",     sampleType="background", systematics=None)
    hist.Add(h_hist_data,  name="data",      sampleType="data",       systematics=None)

    # make the plot
    plot = hist.execute()
    hist.savefig()
    print "  :: Saved plot to "+hist.saveAs+"\n"


## THE END
