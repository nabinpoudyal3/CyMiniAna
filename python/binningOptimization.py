import ROOT
from array import array
import numpy as np
import hepPlotterTools as hpt
import csv

np.set_printoptions(precision=3)

outfile  = open("share/signficance_mttbar.csv","wt")
csvwriter = csv.writer(outfile)
csvwriter.writerow( ('mttbar_bin', 'N Bins', 'Significance (avg)', \
                     'Significance (max)', 'Significance', \
                     'N Events (bckg)', 'N Events (sig)', 'Bins') )


filetype = "pre_signal"
filedir  = "/lustre/umt3/user/demarley/TopAC/AT2.4.19/skims/signal/"

f_stop     = ROOT.TFile.Open(filedir+"singletop_"+filetype+".root")
f_ttbarhad = ROOT.TFile.Open(filedir+"ttbar_had_"+filetype+".root")
f_ttbarlep = ROOT.TFile.Open(filedir+"ttbar_lep_"+filetype+".root")
f_multijet = ROOT.TFile.Open(filedir+"multijet_"+filetype+".root")
f_zprime   = ROOT.TFile.Open(filedir+"zprime_1000_"+filetype+".root")


## delta|y|
# [-3,3]
print " --------------------- "
print " Delta|y| Optimization "
print

middle_bins = [0.1*i for i in range(1,20)]
nbins       = len(middle_bins)
#REBIN
#S/sqrt(B)

# deltay_ttbar_had = f_ttbarhad.h_deltay_signal_nominal
# deltay_ttbar_lep = f_ttbarlep.h_deltay_signal_nominal
# deltay_singletop = f_stop.h_deltay_signal_nominal
# deltay_multijet  = f_multijet.deltay_signal_nominal
# deltay_zprime    = f_zprime.h_deltay_signal_nominal

m_ttbar_deltay_ttbar_had = f_ttbarhad.h_m_ttbar_deltay_signal_nominal
m_ttbar_deltay_ttbar_lep = f_ttbarlep.h_m_ttbar_deltay_signal_nominal
m_ttbar_deltay_singletop = f_stop.h_m_ttbar_deltay_signal_nominal
m_ttbar_deltay_multijet  = f_multijet.m_ttbar_deltay_signal_nominal
m_ttbar_deltay_zprime    = f_zprime.h_m_ttbar_deltay_signal_nominal


mttbar_bins = [
[1,1000],
[1001,1100],
[1101,1300],
[1301,1800],
[1801,5000]
]



for b,bin in enumerate(mttbar_bins):

    print bin[0],bin[1]
    deltay_ttbar_had = m_ttbar_deltay_ttbar_had.ProjectionY("ttbar_had_{0}_{1}".format(bin[0],bin[1]), bin[0],bin[1])
    deltay_ttbar_lep = m_ttbar_deltay_ttbar_lep.ProjectionY("ttbar_lep_{0}_{1}".format(bin[0],bin[1]), bin[0],bin[1])
    deltay_singletop = m_ttbar_deltay_singletop.ProjectionY("singletop_{0}_{1}".format(bin[0],bin[1]), bin[0],bin[1])
    deltay_multijet  = m_ttbar_deltay_multijet.ProjectionY("multijet_{0}_{1}".format(bin[0],bin[1]), bin[0],bin[1])
    deltay_zprime    = m_ttbar_deltay_zprime.ProjectionY("zprime_{0}_{1}".format(bin[0],bin[1]), bin[0],bin[1])


    # 4 bins
    for middle_bin in middle_bins:

        middle_bin_name = str(middle_bin)
        newbins = array('d',[-3.0,-middle_bin,0.0,middle_bin,3.0])

        print " Binning : {0}".format(middle_bin)

        new_ttbar_had = deltay_ttbar_had.Clone()
        new_ttbar_had2 = hpt.hist2list(new_ttbar_had,reBin=newbins,name=middle_bin_name)
        new_ttbar_lep = deltay_ttbar_lep.Clone()
        new_ttbar_lep2 = hpt.hist2list(new_ttbar_lep,reBin=newbins,name=middle_bin_name)
        new_stop = deltay_singletop.Clone()
        new_stop2 = hpt.hist2list(new_stop,reBin=newbins,name=middle_bin_name)
        new_qcd = deltay_multijet.Clone()
        new_qcd2 = hpt.hist2list(new_qcd,reBin=newbins,name=middle_bin_name)
        new_zprime = deltay_zprime.Clone()
        new_zprime2 = hpt.hist2list(new_zprime,reBin=newbins,name=middle_bin_name)

        # s/sqrt(b)
        bckg = new_ttbar_had2['data'] + new_ttbar_lep2['data'] + new_stop2['data'] + new_qcd2['data']
        sig  = np.divide( new_zprime2['data'],np.sqrt(bckg) )
        avg_sig = np.average(sig)

        # ('N Bins', 'Significance (avg)', 'Significance (max)', 'Significance', 'Bins')
        if not np.isnan(avg_sig):
            csvwriter.writerow( (b,4, avg_sig, max(sig), sig, bckg, new_zprime2['data'],newbins) )


    # 6 bins
    for m,left_bin in enumerate(middle_bins):
        for n in range( m+1,nbins ):
            right_bin = middle_bins[n]
            newbins = array('d',[-3.0,-right_bin,-left_bin,0.0,left_bin,right_bin,3.0])
            binname = str(left_bin)+":"+str(right_bin)

            new_ttbar_had = deltay_ttbar_had.Clone()
            new_ttbar_had = hpt.hist2list(new_ttbar_had,reBin=newbins,name=binname)
            new_ttbar_lep = deltay_ttbar_lep.Clone()
            new_ttbar_lep = hpt.hist2list(new_ttbar_lep,reBin=newbins,name=binname)
            new_stop = deltay_singletop.Clone()
            new_stop = hpt.hist2list(new_stop,reBin=newbins,name=binname)
            new_qcd = deltay_multijet.Clone()
            new_qcd = hpt.hist2list(new_qcd,reBin=newbins,name=binname)
            new_zprime = deltay_zprime.Clone()
            new_zprime = hpt.hist2list(new_zprime,reBin=newbins,name=binname)

            # s/sqrt(b)
            bckg = new_ttbar_had['data'] + new_ttbar_lep['data'] + new_stop['data'] + new_qcd['data']
            sig  = np.divide( new_zprime['data'],np.sqrt(bckg) )

            avg_sig = np.average(sig)
            print " Binning : {0}".format(newbins.tolist())

            # ('N Bins', 'Significance (avg)', 'Significance (max)', 'Significance', 'Bins')
            if not np.isnan(avg_sig):
                csvwriter.writerow( (b,6, avg_sig, max(sig), sig, bckg, new_zprime['data'], newbins) )

    # 8 bins
    for m,left_bin in enumerate(middle_bins):
        for n in range(m+1,nbins):
            center_bin = middle_bins[n]
            for k in range(n+1,nbins):
                right_bin = middle_bins[k]
                newbins   = array('d',[-3.0,-right_bin,-center_bin,-left_bin,0.0,\
                                      left_bin,center_bin,right_bin,3.0])

                binname    = str(left_bin)+":"+str(center_bin)+":"+str(right_bin)

                new_ttbar_had = deltay_ttbar_had.Clone()
                new_ttbar_had = hpt.hist2list(new_ttbar_had,reBin=newbins,name=binname)
                new_ttbar_lep = deltay_ttbar_lep.Clone()
                new_ttbar_lep = hpt.hist2list(new_ttbar_lep,reBin=newbins,name=binname)
                new_stop = deltay_singletop.Clone()
                new_stop = hpt.hist2list(new_stop,reBin=newbins,name=binname)
                new_qcd = deltay_multijet.Clone()
                new_qcd = hpt.hist2list(new_qcd,reBin=newbins,name=binname)
                new_zprime = deltay_zprime.Clone()
                new_zprime = hpt.hist2list(new_zprime,reBin=newbins,name=binname)

                # s/sqrt(b)
                bckg = new_ttbar_had['data'] + new_ttbar_lep['data'] + new_stop['data'] + new_qcd['data']
                sig  = np.divide( new_zprime['data'],np.sqrt(bckg) )

                print " Binning : {0}".format(newbins.tolist())

                # ('N Bins', 'Significance (avg)', 'Significance (max)', 'Significance', 'Bins')
                avg_sig = np.average(sig)
                if not np.isnan(avg_sig):
                    csvwriter.writerow( (b,8, avg_sig, max(sig), sig, bckg, new_zprime['data'], newbins) )

outfile.close()

"""




## m_ttbar v delta|y|
# ( [0,5000],[-3,3] )
mttbar_bins = [ 
#array('d',[0,750,900,1300,6000]),
#array('d',[0,750,900,1200,6000]),
#array('d',[0,750,900,1300,2000,6000]),
#array('d',[0,750,900,1200,1600,6000]),
#array('d',[0,750,850,950,1050,1200,1800,6000]),
#array('d',[0,750,900,1050,1300,1800,6000]),
#array('d',[0,750,900,1050,1300,1600,6000]),
array('d',[0,900,1050,1300,1600,2000,6000]),
array('d',[0,900,1050,1300,1800,6000]),
]

#m_ttbar_deltay_ttbar_had = f_ttbarhad.h_m_ttbar_deltay_signal_nominal
#m_ttbar_deltay_ttbar_lep = f_ttbarlep.h_m_ttbar_deltay_signal_nominal
#m_ttbar_deltay_singletop = f_stop.h_m_ttbar_deltay_signal_nominal
#m_ttbar_deltay_zprime    = f_zprime.h_m_ttbar_deltay_signal_nominal

m_ttbar_deltay_ttbar_had = f_ttbarhad.h_m_ttbar_signal_nominal
m_ttbar_deltay_ttbar_lep = f_ttbarlep.h_m_ttbar_signal_nominal
m_ttbar_deltay_singletop = f_stop.h_m_ttbar_signal_nominal
m_ttbar_deltay_zprime    = f_zprime.h_m_ttbar_signal_nominal
m_ttbar_deltay_multijet  = f_multijet.m_ttbar_signal_nominal

for m,mttbar_bin in enumerate(mttbar_bins):

    nbins = len(mttbar_bin)-1
    name  = str(nbins)+str(m)
    print " Binning : {0}".format(mttbar_bin.tolist())

    new_ttbar_had = m_ttbar_deltay_ttbar_had.Clone()
    new_ttbar_had2 = hpt.hist2list(new_ttbar_had,reBin=mttbar_bin,name=name)
    new_ttbar_lep = m_ttbar_deltay_ttbar_lep.Clone()
    new_ttbar_lep2 = hpt.hist2list(new_ttbar_lep,reBin=mttbar_bin,name=name)
    new_stop = m_ttbar_deltay_singletop.Clone()
    new_stop2 = hpt.hist2list(new_stop,reBin=mttbar_bin,name=name)
    new_qcd = m_ttbar_deltay_multijet.Clone()
    new_qcd2 = hpt.hist2list(new_qcd,reBin=mttbar_bin,name=name)
    new_zprime = m_ttbar_deltay_zprime.Clone()
    new_zprime2 = hpt.hist2list(new_zprime,reBin=mttbar_bin,name=name)

    # s/sqrt(b)
    bckg = new_ttbar_had2['data'] + new_ttbar_lep2['data'] + new_stop2['data'] + new_qcd2['data']
    sig  = np.divide( new_zprime2['data'],np.sqrt(bckg) )
    avg_sig = np.average(sig)

    # ('N Bins', 'Significance (avg)', 'Significance (max)', 'Significance', 'N Events (bckg)', 'N Events (sig)', 'Bins')
    csvwriter.writerow( (nbins, avg_sig, max(sig), sig, bckg, new_zprime2['data'], mttbar_bin) )


# ## beta_ttbar v delta|y|
# # ( [0,1],[-5,5] )
# beta_ttbar_deltay_ttbar_had = f_ttbarhad.h_beta_ttbar_deltay_signal_nominal
# beta_ttbar_deltay_ttbar_lep = f_ttbarlep.h_beta_ttbar_deltay_signal_nominal
# beta_ttbar_deltay_singletop = f_stop.h_beta_ttbar_deltay_signal_nominal
# beta_ttbar_deltay_zprime    = f_zprime.h_beta_ttbar_deltay_signal_nominal
# 
# ## pT_ttbar v delta|y|
# # ( [0,200],[-5,5] )
# pt_ttbar_deltay_ttbar_had = f_ttbarhad.h_pt_ttbar_deltay_signal_nominal
# pt_ttbar_deltay_ttbar_lep = f_ttbarlep.h_pt_ttbar_deltay_signal_nominal
# pt_ttbar_deltay_singletop = f_stop.h_pt_ttbar_deltay_signal_nominal
# pt_ttbar_deltay_zprime    = f_zprime.h_pt_ttbar_deltay_signal_nominal
"""
