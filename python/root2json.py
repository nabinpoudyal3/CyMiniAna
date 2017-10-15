"""
Created:        11 November  2016
Last Updated:   11 December  2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109
-----

Script to convert information in root file to format for Keras.

Use 'isbtagged_XX' to get b-tagged jets, but just in case:
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/
        BTaggingBenchmarks#MV2c10_tagger_added_on_22th_May

name            weight cut    b-jet efficiency [%]   purity [%]   c RR    tau RR   light RR
FixedCutBEff_70    0.6455     70.00                  93.55        7.09    17.25    119.69
FixedCutBEff_77    0.3706     77.00                  89.81        4.21    8.09     57.90


save large-R jet information with selection:
  pT>300 GeV
  matched to parton
"""
import sys
import json
import ROOT
from math import fabs
from copy import deepcopy
from argparse import ArgumentParser
from time import strftime,localtime

class LorentzVector(object):
    def __init__(self):
        self.p4 = ROOT.TLorentzVector()
        self.charge = 0

        return

class Ljet(LorentzVector):
    def __init__(self):
        LorentzVector.__init__(self)
        self.target = -1
        self.tjets  = []
        self.tjets_b  = []
        self.tjets_nb = []

        return



def deltaR(tlvA,tlvB,dR=0.75):
    """Determine if two objects are matched"""
    result = False
    try:
        result = (tlvA.p4.DeltaR(tlvB.p4)<dR)
    except:
        return False

    return result




parser = ArgumentParser(description="Root2JSON4Keras")

parser.add_argument('-f','--files', action='store',
                    dest='listOfFiles',default='examples/share/listOfFiles.txt',
                    help='Name of file that contains root files to plot')
parser.add_argument('-o','--outpath', action='store',
                    dest='outpath',default='data/',
                    help='Directory for storing output plots')
results = parser.parse_args()

listOfFiles = open( results.listOfFiles,'r').readlines()
outpath     = results.outpath
date        = strftime("%d%b", localtime())


## -- arguments for the ljet selection -- ##
ljet_charge_max = 5.
btag_wkpt   = "77"
ljet_pt_cut = 300000.
tjet_pt_cut = 10000.
deltaR_tru  = 0.75
deltaR_tjet = 0.8     # ljet R = 1.0; tjet R = 0.2
t_index     = 1       # +2/3 charge
tbar_index  = 0       # -2/3 charge
parton_def  = 'afterFSR'
success     = '\x01'
eventLevel  = False   # eventLevel or objectLevel DNN
## ---- ##

# record this information into output file
metadata = {'ljet_pt':ljet_pt_cut,'btag_wkpt':btag_wkpt,'tjet_pt':tjet_pt_cut,\
            'deltaR_ljet_parton':deltaR_tru,'deltaR_ljet_tjet':deltaR_tjet,\
            't_target':t_index,'tbar_target':tbar_index,\
            'parton':parton_def,'ljet_charge_max':ljet_charge_max}


for file in listOfFiles:
    print " Running over ",file.rstrip('\n')
    f   = ROOT.TFile.Open(file.rstrip('\n'))
    nom = f.Get("nominal")
    tru = f.Get("truth")

    # ljet pt,eta,charge
    # ljet -- btagged subjets charge (sum)
    # ljet -- not btagged subjets charge (sum)
    # ljet -- highest pT subjet charge
    data = {'pt':[],'eta':[],
            'charge':[],\
            'subjets_btag_charge': [],\
            'subjets_nbtag_charge':[],\
            'subjets_delta_btag_nbtag':[],\
            'subjet_pt0_charge':[],\
            'subjet_pt0_btag_charge':[],\
            'subjet_pt0_nbtag_charge':[],\
            'reco_m_ttbar':[],
            'truth_m_ttbar':[],\
            'target':[]}

    metadata['file'] = file
    data['metadata'] = metadata

    output = outpath+file.split("/")[-1].split(".root")[0]+"_"+btag_wkpt+"_"+date
    print " > Saving output to "+output


    ## -- Event Loop -- ##
    tru_entry    = 0
    nEntries_tru = tru.GetEntries()

    print " Processing "+str(nom.GetEntries())+" entries"
    for entry in xrange(nom.GetEntries()):
        nom.GetEntry(entry)
        if not entry%10000: print "   > Entry ",entry

        ## -- Build Truth-level objects
        skip_event = False
        finding_parton_entry = True
        while finding_parton_entry:
            if tru_entry > nEntries_tru:
                break
            tru.GetEntry(tru_entry)
            if tru.eventNumber == nom.eventNumber: 
                finding_parton_entry = False
            else:
                tru_entry+=1
        tru.GetEntry(tru_entry)

        # If both tops decay leptonically, we can skip the event
        # -- we want to look at the tops that decay hadronically (allhad & l+jets)!
        if fabs( tru.MC_Wdecay1_from_tbar_pdgId ) > 5 and fabs( tru.MC_Wdecay1_from_t_pdgId ) > 5:
            continue

        # after FSR: W & b
        MC_b_from_t = ROOT.TLorentzVector()
        MC_W_from_t = ROOT.TLorentzVector()
        truth_t     = None
        MC_b_from_tbar = ROOT.TLorentzVector()
        MC_W_from_tbar = ROOT.TLorentzVector()
        truth_tbar     = None
        if parton_def == 'afterFSR':
            # only save information if the W decays hadronically
            tmp_t = LorentzVector()
            if fabs(tru.MC_Wdecay1_from_t_pdgId)<=5 and tru.MC_b_from_t_pt>1e-6 and tru.MC_W_from_t_pt>1e-6:
                MC_b_from_t.SetPtEtaPhiM( tru.MC_b_from_t_pt,tru.MC_b_from_t_eta,tru.MC_b_from_t_phi,tru.MC_b_from_t_m )
                MC_W_from_t.SetPtEtaPhiM( tru.MC_W_from_t_pt,tru.MC_W_from_t_eta,tru.MC_W_from_t_phi,tru.MC_W_from_t_m )
                tmp_t.p4 = MC_W_from_t + MC_b_from_t
                truth_t  = deepcopy(tmp_t)

            # only save information if the W decays hadronically
            tmp_tbar = LorentzVector()
            if fabs(tru.MC_Wdecay1_from_tbar_pdgId)<=5 and tru.MC_b_from_tbar_pt>1e-6 and tru.MC_W_from_tbar_pt>1e-6:
                MC_b_from_tbar.SetPtEtaPhiM( tru.MC_b_from_tbar_pt,tru.MC_b_from_tbar_eta,tru.MC_b_from_tbar_phi,tru.MC_b_from_tbar_m )
                MC_W_from_tbar.SetPtEtaPhiM( tru.MC_W_from_tbar_pt,tru.MC_W_from_tbar_eta,tru.MC_W_from_tbar_phi,tru.MC_W_from_tbar_m )
                tmp_tbar.p4 = MC_W_from_tbar + MC_b_from_tbar
                truth_tbar = deepcopy(tmp_tbar)

        if truth_t is None and truth_tbar is None:
            print " Truth tops are None "
            continue
        truth_m_ttbar = tru.MC_ttbar_afterFSR_m



        ## -- Build the reco-level large-R jets
        ljet_t    = None
        ljet_tbar = None
        for l,lj in enumerate(nom.ljet_pt):
            if lj < ljet_pt_cut: continue        # pT cut

            ljet = Ljet()
            ljet.p4.SetPtEtaPhiE(lj,nom.ljet_eta[l],nom.ljet_phi[l],nom.ljet_e[l])

            if fabs( nom.ljet_charge[l] ) > ljet_charge_max:
                continue
            ljet.charge = nom.ljet_charge[l]

            if deltaR( truth_t, ljet, deltaR_tru ):
                ljet.target = t_index
                ljet_t = deepcopy(ljet)
            elif deltaR( truth_tbar, ljet, deltaR_tru ):
                ljet.target = tbar_index
                ljet_tbar = deepcopy(ljet)

        # Need at least one jet truth matched
        if ljet_t is None and ljet_tbar is None:
            print " Both tops are None "
            continue

        try:
            reco_m_ttbar = (ljet_t.p4+ljet_tbar.p4).M()
        except AttributeError:
            reco_m_ttbar = -1



        ## -- Build the reco-level track jets
        tjets = getattr( nom,"tjet_isbtagged_"+btag_wkpt )

        for tj,charge in enumerate(nom.tjet_charge):
            if nom.tjet_pt[tj]<tjet_pt_cut: continue  # add pT cut because ntuples are wrong

            tjet = LorentzVector()
            tjet.p4.SetPtEtaPhiE( nom.tjet_pt[tj], nom.tjet_eta[tj], nom.tjet_phi[tj], nom.tjet_e[tj] )
            tjet.charge = charge

            # Track Jet matched to large-R jet (top)
            if deltaR( ljet_t, tjet, deltaR_tjet):
                ljet_t.tjets.append(tjet)

                # b-tagging
                if tjets[tj] == success:
                    ljet_t.tjets_b.append(tjet)
                else:
                    ljet_t.tjets_nb.append(tjet)


            # Track Jet matched to large-R jet (antitop)
            if deltaR( ljet_tbar, tjet, deltaR_tjet):
                ljet_tbar.tjets.append(tjet)

                # b-tagging
                if tjets[tj] == success:
                    ljet_tbar.tjets_b.append(tjet)
                else:
                    ljet_tbar.tjets_nb.append(tjet)


        ### SIMPLE BINARY CLASSIFICATION ###
        ## Top
        if ljet_t is not None and len(ljet_t.tjets_b)>0:
            subjets_btag_Q  = sum([i.charge for i in ljet_t.tjets_b])
            subjets_nbtag_Q = sum([i.charge for i in ljet_t.tjets_nb])

            data['pt'].append( ljet_t.p4.Pt())
            data['eta'].append(ljet_t.p4.Eta())
            data['charge'].append(ljet_t.charge)
            data['subjets_btag_charge'].append(  subjets_btag_Q  )
            data['subjets_nbtag_charge'].append( subjets_nbtag_Q  )
            data['subjets_delta_btag_nbtag'].append(subjets_btag_Q-subjets_nbtag_Q)
            data['subjet_pt0_charge'].append(       ljet_t.tjets[0].charge)
            data['subjet_pt0_btag_charge'].append(  ljet_t.tjets_b[0].charge)
            if len(ljet_t.tjets_nb)>0:
                data['subjet_pt0_nbtag_charge'].append( ljet_t.tjets_nb[0].charge)
            else:
                data['subjet_pt0_nbtag_charge'].append( 0.0 )
            data['reco_m_ttbar'].append(  reco_m_ttbar )
            data['truth_m_ttbar'].append( truth_m_ttbar )
            data['target'].append(ljet_t.target)
#        else:
#            print " No top candidate: {0}".format(ljet_t)

        ## Anti-top
        if ljet_tbar is not None and len(ljet_tbar.tjets_b)>0:
            subjets_btag_Q  = sum([i.charge for i in ljet_tbar.tjets_b])
            subjets_nbtag_Q = sum([i.charge for i in ljet_tbar.tjets_nb])

            data['pt'].append( ljet_tbar.p4.Pt())
            data['eta'].append(ljet_tbar.p4.Eta())
            data['charge'].append(ljet_tbar.charge)
            data['subjets_btag_charge'].append(  subjets_btag_Q  )
            data['subjets_nbtag_charge'].append( subjets_nbtag_Q  )
            data['subjets_delta_btag_nbtag'].append(subjets_btag_Q-subjets_nbtag_Q)
            data['subjet_pt0_charge'].append(       ljet_tbar.tjets[0].charge  )
            data['subjet_pt0_btag_charge'].append(  ljet_tbar.tjets_b[0].charge  )
            if len(ljet_tbar.tjets_nb)>0:
                data['subjet_pt0_nbtag_charge'].append( ljet_tbar.tjets_nb[0].charge  )
            else:
                data['subjet_pt0_nbtag_charge'].append( 0.0 )
            data['reco_m_ttbar'].append(  reco_m_ttbar )
            data['truth_m_ttbar'].append( truth_m_ttbar )
            data['target'].append(ljet_tbar.target)
#        else:
#            print " No antitop candidate: {0}".format(ljet_tbar)

    print " Writing data "+output+'.json'
    with open(output+'.json', 'w') as outfile:
        json.dump(data, outfile)


## THE END ##
