import ROOT
from multijetEstimation import MultijetEstimation
from array import *




def fillHists(output_histograms,h_tmp,histname,region):
    """Fill histograms for output"""
    if "_ttbar_deltay_" in histname:
        print " 2D histograms: {0}".format(h_tmp.Integral())
        for bin in range(1,h_tmp.GetNbinsX()+1):
            for ybin in range(h_tmp.GetNbinsY()+1):
                binNumber = h_tmp.GetBin(bin,ybin)
                output_histograms[histname.format(region)].AddBinContent(binNumber,h_tmp.GetBinContent(binNumber))
    else:
        for bin in range(1,h_tmp.GetNbinsX()+1):
            output_histograms[histname.format(region)].AddBinContent(bin,h_tmp.GetBinContent(bin))

    return


## List of histograms
hists = [
"h_n_btags_{0}_nominal",\
"h_n_trackjets_{0}_nominal",\
"h_n_topCandidates_{0}_nominal",\
"h_n_ljets_{0}_nominal",\
"h_deltaR_t1_t2_{0}_nominal",\
"h_ljet_m_t1_{0}_nominal",\
"h_ljet_m_t2_{0}_nominal",\
"h_ljet_eta_t1_{0}_nominal",\
"h_ljet_eta_t2_{0}_nominal",\
"h_ljet_pt_t1_{0}_nominal",\
"h_ljet_pt_t2_{0}_nominal",\
"h_ljet_tau32wta_t1_{0}_nominal",\
"h_ljet_tau32wta_t2_{0}_nominal",\
"h_ljet_split23_t1_{0}_nominal",\
"h_ljet_split23_t2_{0}_nominal",\
"h_ljet_charge_t1_{0}_nominal",\
"h_ljet_charge_t2_{0}_nominal",\
"h_ljet_btagged_tjet_charge_t1_{0}_nominal",\
"h_ljet_btagged_tjet_charge_t2_{0}_nominal",\
"h_ljet_nbtagged_tjet_charge_t1_{0}_nominal",\
"h_ljet_nbtagged_tjet_charge_t2_{0}_nominal",\
"h_ljet_deltaQ_btag_nbtag_t1_{0}_nominal",\
"h_ljet_deltaQ_btag_nbtag_t2_{0}_nominal",\
"h_ljet_subjet_pt0_charge_t1_{0}_nominal",\
"h_ljet_subjet_pt0_charge_t2_{0}_nominal",\
"h_ljet_subjet_pt0_btag_charge_t1_{0}_nominal",\
"h_ljet_subjet_pt0_btag_charge_t2_{0}_nominal",\
"h_ljet_subjet_pt0_nbtag_charge_t1_{0}_nominal",\
"h_ljet_subjet_pt0_nbtag_charge_t2_{0}_nominal",\
"h_deltaR_t_tbar_{0}_nominal",\
"h_ljet_m_t_{0}_nominal",\
"h_ljet_m_tbar_{0}_nominal",\
"h_ljet_eta_t_{0}_nominal",\
"h_ljet_eta_tbar_{0}_nominal",\
"h_ljet_pt_t_{0}_nominal",\
"h_ljet_pt_tbar_{0}_nominal",\
"h_ljet_tau32wta_t_{0}_nominal",\
"h_ljet_tau32wta_tbar_{0}_nominal",\
"h_ljet_split23_t_{0}_nominal",\
"h_ljet_split23_tbar_{0}_nominal",\
"h_ljet_DNN_t_{0}_nominal",\
"h_ljet_DNN_tbar_{0}_nominal",\
"h_ljet_charge_t_{0}_nominal",\
"h_ljet_charge_tbar_{0}_nominal",\
"h_ljet_btagged_tjet_charge_t_{0}_nominal",\
"h_ljet_btagged_tjet_charge_tbar_{0}_nominal",\
"h_ljet_nbtagged_tjet_charge_t_{0}_nominal",\
"h_ljet_nbtagged_tjet_charge_tbar_{0}_nominal",\
"h_ljet_deltaQ_btag_nbtag_t_{0}_nominal",\
"h_ljet_deltaQ_btag_nbtag_tbar_{0}_nominal",\
"h_ljet_subjet_pt0_charge_t_{0}_nominal",\
"h_ljet_subjet_pt0_charge_tbar_{0}_nominal",\
"h_ljet_subjet_pt0_btag_charge_t_{0}_nominal",\
"h_ljet_subjet_pt0_btag_charge_tbar_{0}_nominal",\
"h_ljet_subjet_pt0_nbtag_charge_t_{0}_nominal",\
"h_ljet_subjet_pt0_nbtag_charge_tbar_{0}_nominal",\
"h_ljet_t50_{0}_nominal",\
"h_ljet_t80_{0}_nominal",\
"h_ljet_m_{0}_nominal",\
"h_ljet_eta_{0}_nominal",\
"h_ljet_pt_{0}_nominal",\
"h_ljet_tau32wta_{0}_nominal",\
"h_ljet_split23_{0}_nominal",\
"h_deltay_{0}_nominal",\
"h_beta_ttbar_{0}_nominal",\
"h_pt_ttbar_{0}_nominal",\
"h_m_ttbar_{0}_nominal",\
"h_beta_ttbar_deltay_{0}_nominal",\
"h_pt_ttbar_deltay_{0}_nominal",\
"h_m_ttbar_deltay_{0}_nominal",\
]



## ABCD method regions
naming   = {"0b0t":"A",\
            "0b1t":"D",\
            "0b2t":"G",\
            "1b0t":"B",\
            "1b1t":"E",\
            "1b2t":"H",\
            "2b0t":"C",\
            "2b1t":"F",\
            "2b2t":"S"}
qcd_regions = naming.keys()
# qcd_regions.sort() # need the regions sorted in order to process regions
                   # other than validation and signal
regions  = qcd_regions+["validation","signal"]
yields   = {}
outHists = {}
qcdSignal = MultijetEstimation()  # signal region calculation


## Load files to estimate QCD
filetype   = "qcd"
filedir    = "/home/majersky/AllHadAC/CyMiniAnaAC/output_qcd/"

f_data     = ROOT.TFile.Open(filedir+"data_"+filetype+".root")
f_stop     = ROOT.TFile.Open(filedir+"singletop_"+filetype+".root")
f_ttbarhad = ROOT.TFile.Open(filedir+"ttbar_had_"+filetype+".root")
f_ttbarlep = ROOT.TFile.Open(filedir+"ttbar_lep_"+filetype+".root")


## output
outFile = ROOT.TFile("multijet_"+filetype+".root","recreate")  # QCD file

# uneven binning for some histos
beta_ttbar_bins = array('d', [0., 0.3, 0.6, 1.0])
pt_ttbar_bins = array('d', [0., 50., 100., 200.])
m_ttbar_bins = array('d', [0., 500., 800., 1000., 1200., 1500., 2000., 2500., 4000., 5000.])

for region in regions:
    outHists["h_deltay_"+region+"_nominal"]     = ROOT.TH1D("deltay_"+region+"_nominal","deltay_"+region+"_nominal",          12, -3.0,   3.0)
    # outHists["h_beta_ttbar_"+region+"_nominal"] = ROOT.TH1D("beta_ttbar_"+region+"_nominal","beta_ttbar_"+region+"_nominal",      100,  0.0,   1.0)
    outHists["h_beta_ttbar_"+region+"_nominal"] = ROOT.TH1D("beta_ttbar_"+region+"_nominal","beta_ttbar_"+region+"_nominal", len(beta_ttbar_bins)-1, beta_ttbar_bins)
    # outHists["h_pt_ttbar_"+region+"_nominal"]   = ROOT.TH1D("pt_ttbar_"+region+"_nominal","pt_ttbar_"+region+"_nominal",       2000,  0.0, 200.0)
    outHists["h_pt_ttbar_"+region+"_nominal"]   = ROOT.TH1D("pt_ttbar_"+region+"_nominal","pt_ttbar_"+region+"_nominal", len(pt_ttbar_bins)-1, pt_ttbar_bins)
    outHists["h_m_ttbar_"+region+"_nominal"]    = ROOT.TH1D("m_ttbar_"+region+"_nominal","m_ttbar_"+region+"_nominal", len(m_ttbar_bins)-1, m_ttbar_bins)
    outHists["h_beta_ttbar_deltay_"+region+"_nominal"] = ROOT.TH2D("beta_ttbar_deltay_"+region+"_nominal","beta_ttbar_deltay_"+region+"_nominal",  100,  0.0,   1.0, 600, -3.0, 3.0)
    outHists["h_pt_ttbar_deltay_"+region+"_nominal"]   = ROOT.TH2D("pt_ttbar_deltay_"+region+"_nominal","pt_ttbar_deltay_"+region+"_nominal",   2000,  0.0, 200.0, 600, -3.0, 3.0)
    outHists["h_m_ttbar_deltay_"+region+"_nominal"]    = ROOT.TH2D("m_ttbar_deltay_"+region+"_nominal","m_ttbar_deltay_"+region+"_nominal",    5000,  0.0,5000.0, 600, -3.0, 3.0)
    outHists["h_n_btags_"+region+"_nominal"]           = ROOT.TH1D("n_btags_"+region+"_nominal","n_btags_"+region+"_nominal",          11, -0.5,  10.5)
    outHists["h_n_trackjets_"+region+"_nominal"]       = ROOT.TH1D("n_trackjets_"+region+"_nominal","n_trackjets_"+region+"_nominal",      21, -0.5,  20.5)
    outHists["h_n_topCandidates_"+region+"_nominal"]   = ROOT.TH1D("n_topCandidates_"+region+"_nominal","n_topCandidates_"+region+"_nominal",   6, -0.5,   5.5)
    outHists["h_n_ljets_"+region+"_nominal"]       = ROOT.TH1D("n_ljets_"+region+"_nominal","n_ljets_"+region+"_nominal",          11, -0.5,  10.5)

    outHists["h_deltaR_t1_t2_"+region+"_nominal"]  = ROOT.TH1D("deltaR_t1_t2_"+region+"_nominal","deltaR_t1_t2_"+region+"_nominal",    20,  0.0,   5.0)
    outHists["h_ljet_m_t1_"+region+"_nominal"]     = ROOT.TH1D("ljet_m_t1_"+region+"_nominal","ljet_m_t1_"+region+"_nominal",       25, 50.0, 400.0)
    outHists["h_ljet_m_t2_"+region+"_nominal"]     = ROOT.TH1D("ljet_m_t2_"+region+"_nominal","ljet_m_t2_"+region+"_nominal",       25, 50.0, 400.0)
    outHists["h_ljet_eta_t1_"+region+"_nominal"]   = ROOT.TH1D("ljet_eta_t1_"+region+"_nominal","ljet_eta_t1_"+region+"_nominal",      10, -2.0,   2.0)
    outHists["h_ljet_eta_t2_"+region+"_nominal"]   = ROOT.TH1D("ljet_eta_t2_"+region+"_nominal","ljet_eta_t2_"+region+"_nominal",      10, -2.0,   2.0)
    outHists["h_ljet_pt_t1_"+region+"_nominal"]    = ROOT.TH1D("ljet_pt_t1_"+region+"_nominal","ljet_pt_t1_"+region+"_nominal",     20,200.0,3000.0)
    outHists["h_ljet_pt_t2_"+region+"_nominal"]    = ROOT.TH1D("ljet_pt_t2_"+region+"_nominal","ljet_pt_t2_"+region+"_nominal",     20,200.0,3000.0)
    outHists["h_ljet_tau32wta_t1_"+region+"_nominal"]  = ROOT.TH1D("ljet_tau32wta_t1_"+region+"_nominal","ljet_tau32wta_t1_"+region+"_nominal",100,0.0,10.0)
    outHists["h_ljet_tau32wta_t2_"+region+"_nominal"]  = ROOT.TH1D("ljet_tau32wta_t2_"+region+"_nominal","ljet_tau32wta_t2_"+region+"_nominal",100,0.0,10.0)
    outHists["h_ljet_split23_t1_"+region+"_nominal"]  = ROOT.TH1D("ljet_split23_t1_"+region+"_nominal","ljet_split23_t1_"+region+"_nominal",10,0.0,100.0)
    outHists["h_ljet_split23_t2_"+region+"_nominal"]  = ROOT.TH1D("ljet_split23_t2_"+region+"_nominal","ljet_split23_t2_"+region+"_nominal",10,0.0,100.0)
    outHists["h_ljet_charge_t1_"+region+"_nominal"]    = ROOT.TH1D("ljet_charge_t1_"+region+"_nominal", "ljet_charge_t1_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_charge_t2_"+region+"_nominal"]    = ROOT.TH1D("ljet_charge_t2_"+region+"_nominal", "ljet_charge_t2_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_btagged_tjet_charge_t1_"+region+"_nominal"]  = ROOT.TH1D("ljet_btagged_tjet_charge_t1_"+region+"_nominal", "ljet_btagged_tjet_charge_t1_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_btagged_tjet_charge_t2_"+region+"_nominal"]  = ROOT.TH1D("ljet_btagged_tjet_charge_t2_"+region+"_nominal", "ljet_btagged_tjet_charge_t2_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_nbtagged_tjet_charge_t1_"+region+"_nominal"] = ROOT.TH1D("ljet_nbtagged_tjet_charge_t1_"+region+"_nominal","ljet_nbtagged_tjet_charge_t1_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_nbtagged_tjet_charge_t2_"+region+"_nominal"] = ROOT.TH1D("ljet_nbtagged_tjet_charge_t2_"+region+"_nominal","ljet_nbtagged_tjet_charge_t2_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_deltaQ_btag_nbtag_t1_"+region+"_nominal"] = ROOT.TH1D("ljet_deltaQ_btag_nbtag_t1_"+region+"_nominal","ljet_deltaQ_btag_nbtag_t1_"+region+"_nominal", 100, -20.,20.);
    outHists["h_ljet_deltaQ_btag_nbtag_t2_"+region+"_nominal"] = ROOT.TH1D("ljet_deltaQ_btag_nbtag_t2_"+region+"_nominal","ljet_deltaQ_btag_nbtag_t2_"+region+"_nominal", 100, -20.,20.);
    outHists["h_ljet_subjet_pt0_charge_t1_"+region+"_nominal"] = ROOT.TH1D("ljet_subjet_pt0_charge_t1_"+region+"_nominal","ljet_subjet_pt0_charge_t1_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_subjet_pt0_charge_t2_"+region+"_nominal"] = ROOT.TH1D("ljet_subjet_pt0_charge_t2_"+region+"_nominal","ljet_subjet_pt0_charge_t2_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_subjet_pt0_btag_charge_t1_"+region+"_nominal"]  = ROOT.TH1D("ljet_subjet_pt0_btag_charge_t1_"+region+"_nominal","ljet_subjet_pt0_btag_charge_t1_"+region+"_nominal",  100, -10.,10.);
    outHists["h_ljet_subjet_pt0_btag_charge_t2_"+region+"_nominal"]  = ROOT.TH1D("ljet_subjet_pt0_btag_charge_t2_"+region+"_nominal","ljet_subjet_pt0_btag_charge_t2_"+region+"_nominal",  100, -10.,10.);
    outHists["h_ljet_subjet_pt0_nbtag_charge_t1_"+region+"_nominal"] = ROOT.TH1D("ljet_subjet_pt0_nbtag_charge_t1_"+region+"_nominal","ljet_subjet_pt0_nbtag_charge_t1_"+region+"_nominal",  100, -10.,10.);
    outHists["h_ljet_subjet_pt0_nbtag_charge_t2_"+region+"_nominal"] = ROOT.TH1D("ljet_subjet_pt0_nbtag_charge_t2_"+region+"_nominal","ljet_subjet_pt0_nbtag_charge_t2_"+region+"_nominal", 100, -10.,10.);

    outHists["h_deltaR_t_tbar_"+region+"_nominal"] = ROOT.TH1D("deltaR_t_tbar_"+region+"_nominal","deltaR_t_tbar_"+region+"_nominal",   20,  0.0,   5.0)
    outHists["h_ljet_m_t_"+region+"_nominal"]      = ROOT.TH1D("ljet_m_t_"+region+"_nominal","ljet_m_t_"+region+"_nominal",        25, 50.0, 400.0)
    outHists["h_ljet_m_tbar_"+region+"_nominal"]   = ROOT.TH1D("ljet_m_tbar_"+region+"_nominal","ljet_m_tbar_"+region+"_nominal",  25, 50.0, 400.0)
    outHists["h_ljet_eta_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_eta_t_"+region+"_nominal","ljet_eta_t_"+region+"_nominal",       10, -2.0,   2.0)
    outHists["h_ljet_eta_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_eta_tbar_"+region+"_nominal","ljet_eta_tbar_"+region+"_nominal", 10, -2.0,   2.0)
    outHists["h_ljet_pt_t_"+region+"_nominal"]     = ROOT.TH1D("ljet_pt_t_"+region+"_nominal","ljet_pt_t_"+region+"_nominal",       20,200.0,3000.0)
    outHists["h_ljet_pt_tbar_"+region+"_nominal"]  = ROOT.TH1D("ljet_pt_tbar_"+region+"_nominal","ljet_pt_tbar_"+region+"_nominal", 20,200.0,3000.0)
    outHists["h_ljet_tau32wta_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_tau32wta_t_"+region+"_nominal","ljet_tau32wta_t_"+region+"_nominal",      100,0.0,10.0)
    outHists["h_ljet_tau32wta_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_tau32wta_tbar_"+region+"_nominal","ljet_tau32wta_tbar_"+region+"_nominal",100,0.0,10.0)
    outHists["h_ljet_split23_t_"+region+"_nominal"]     = ROOT.TH1D("ljet_split23_t_"+region+"_nominal","ljet_split23_t_"+region+"_nominal",      10,0.0,100.0)
    outHists["h_ljet_split23_tbar_"+region+"_nominal"]  = ROOT.TH1D("ljet_split23_tbar_"+region+"_nominal","ljet_split23_tbar_"+region+"_nominal",10,0.0,100.0)
    outHists["h_ljet_DNN_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_DNN_t_"+region+"_nominal","ljet_DNN_t_"+region+"_nominal",        10,  0.0, 1.0)
    outHists["h_ljet_DNN_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_DNN_tbar_"+region+"_nominal","ljet_DNN_tbar_"+region+"_nominal",  10,  0.0, 1.0)

    outHists["h_ljet_charge_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_charge_t_"+region+"_nominal", "ljet_charge_t_"+region+"_nominal",       100, -10.,10.);
    outHists["h_ljet_charge_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_charge_tbar_"+region+"_nominal", "ljet_charge_tbar_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_btagged_tjet_charge_t_"+region+"_nominal"]     = ROOT.TH1D("ljet_btagged_tjet_charge_t_"+region+"_nominal", "ljet_btagged_tjet_charge_t_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_btagged_tjet_charge_tbar_"+region+"_nominal"]  = ROOT.TH1D("ljet_btagged_tjet_charge_tbar_"+region+"_nominal", "ljet_btagged_tjet_charge_tbar_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_nbtagged_tjet_charge_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_nbtagged_tjet_charge_t_"+region+"_nominal","ljet_nbtagged_tjet_charge_t_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_nbtagged_tjet_charge_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_nbtagged_tjet_charge_tbar_"+region+"_nominal","ljet_nbtagged_tjet_charge_tbar_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_deltaQ_btag_nbtag_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_deltaQ_btag_nbtag_t_"+region+"_nominal","ljet_deltaQ_btag_nbtag_t_"+region+"_nominal", 100, -20.,20.);
    outHists["h_ljet_deltaQ_btag_nbtag_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_deltaQ_btag_nbtag_tbar_"+region+"_nominal","ljet_deltaQ_btag_nbtag_tbar_"+region+"_nominal", 100, -20.,20.);
    outHists["h_ljet_subjet_pt0_charge_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_subjet_pt0_charge_t_"+region+"_nominal","ljet_subjet_pt0_charge_t_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_subjet_pt0_charge_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_subjet_pt0_charge_tbar_"+region+"_nominal","ljet_subjet_pt0_charge_tbar_"+region+"_nominal", 100, -10.,10.);
    outHists["h_ljet_subjet_pt0_btag_charge_t_"+region+"_nominal"]     = ROOT.TH1D("ljet_subjet_pt0_btag_charge_t_"+region+"_nominal","ljet_subjet_pt0_btag_charge_t_"+region+"_nominal",  100, -10.,10.);
    outHists["h_ljet_subjet_pt0_btag_charge_tbar_"+region+"_nominal"]  = ROOT.TH1D("ljet_subjet_pt0_btag_charge_tbar_"+region+"_nominal","ljet_subjet_pt0_btag_charge_tbar_"+region+"_nominal",  100, -10.,10.);
    outHists["h_ljet_subjet_pt0_nbtag_charge_t_"+region+"_nominal"]    = ROOT.TH1D("ljet_subjet_pt0_nbtag_charge_t_"+region+"_nominal","ljet_subjet_pt0_nbtag_charge_t_"+region+"_nominal",  100, -10.,10.);
    outHists["h_ljet_subjet_pt0_nbtag_charge_tbar_"+region+"_nominal"] = ROOT.TH1D("ljet_subjet_pt0_nbtag_charge_tbar_"+region+"_nominal","ljet_subjet_pt0_nbtag_charge_tbar_"+region+"_nominal", 100, -10.,10.);

    outHists["h_ljet_t50_"+region+"_nominal"] = ROOT.TH1D("ljet_t50_"+region+"_nominal","ljet_t50_"+region+"_nominal",          2, -0.5,   1.5)
    outHists["h_ljet_t80_"+region+"_nominal"] = ROOT.TH1D("ljet_t80_"+region+"_nominal","ljet_t80_"+region+"_nominal",          2, -0.5,   1.5)
    outHists["h_ljet_m_"+region+"_nominal"] = ROOT.TH1D("ljet_m_"+region+"_nominal","ljet_m_"+region+"_nominal",           25, 50.0, 400.0)
    outHists["h_ljet_eta_"+region+"_nominal"] = ROOT.TH1D("ljet_eta_"+region+"_nominal","ljet_eta_"+region+"_nominal",         10, -2.0,   2.0)
    outHists["h_ljet_pt_"+region+"_nominal"] = ROOT.TH1D("ljet_pt_"+region+"_nominal","ljet_pt_"+region+"_nominal",        20,200.0,3000.0)
    outHists["h_ljet_tau32wta_"+region+"_nominal"]  = ROOT.TH1D("ljet_tau32wta_"+region+"_nominal","ljet_tau32wta_"+region+"_nominal",100,0.0,10.0)
    outHists["h_ljet_split23_"+region+"_nominal"]  = ROOT.TH1D("ljet_split23_"+region+"_nominal","ljet_split23_"+region+"_nominal",10,0.0,100.0)


    ## loop over histograms
    for hist in hists:

        if region=="signal":
            print " REGION-{0} : {2} ".format(region,h_data.Integral(),hist)
            tmp = qcdSignal.calcQCDsignal( outHists,hist )
            fillHists(outHists,tmp,hist,region)
        elif region=="validation":
            print " REGION-{0} : {2} ".format(region,h_data.Integral(),hist)
            tmp = qcdSignal.calcQCDvalidation( outHists,hist )
            fillHists(outHists,tmp,hist,region)
        else:
            h_data = f_data.Get(hist.format(region))
            h_stop = f_stop.Get(hist.format(region))
            h_ttbarhad = f_ttbarhad.Get(hist.format(region))
            h_ttbarlep = f_ttbarlep.Get(hist.format(region))


            ## Total MC background
            h_stop.Add(h_ttbarhad)
            h_stop.Add(h_ttbarlep)

            mc_background = h_stop.Clone()

            ## Subtract from Data
            h_data.Add(mc_background,-1)

            fillHists(outHists,h_data,hist,region)

            print " REGION-{0} : {2} = {1}".format(region,h_data.Integral(),hist)
            yields[ naming[region] ] = h_data.Integral()

##BCKG = 0.5*(yields["G"]/yields["A"] + yields["H"]/yields["B"])*yields["C"]
##print BCKG


outFile.Write()
outFile.Close()
