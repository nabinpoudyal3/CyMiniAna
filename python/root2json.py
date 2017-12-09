"""
Created:        11 November  2016
Last Updated:   11 December  2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Script to convert information in root file to JSON format for Keras.

save large-R jet information with selection:

features:
"""
import sys
import json
import ROOT
import util
from math import fabs
from copy import deepcopy
from time import strftime,localtime


class Root2json(object):
    """Class for converting ROOT information into JSON format for NN"""
    def __init__(self):

        # Setup verbose output
        self.verbose_level = "INFO"

        # Basic values needed
        self.outpath     = "./"
        self.listOfFiles = "examples/share/listOfFiles.txt"
        self.metadata    = {} # set in initialize()
        self.date        = strftime("%d%b", localtime())
        self.outdata     = {} # set in clear_output()
        self.nEntries    = -1

        ## -- arguments for the ljet selection -- ##
        self.ljet_charge_max = 5.
        self.btag_wkpt    = "77"
        self.ljet_pt_cut  = 300000.
        self.ljet_eta_cut = 2.
        self.tjet_pt_cut  = 10000.
        self.deltaR_tru   = 0.75
        self.deltaR_tjet  = 0.8     # ljet R = 1.0; tjet R = 0.2
        self.t_index      = 1       # +2/3 charge
        self.tbar_index   = 0       # -2/3 charge
        self.parton_def   = 'afterFSR'
        self.success      = '\x01'
        self.nsubjets     = 3
        ## ---- ##


    def initialize(self):
        """Setup a few things"""
        # record this information into output file
        self.metadata = {'ljet_pt':self.ljet_pt_cut,'ljet_eta':self.ljet_eta_cut,\
                         'btag_wkpt':self.btag_wkpt,'tjet_pt':self.tjet_pt_cut,\
                         'deltaR_ljet_parton':self.deltaR_tru,\
                         'deltaR_ljet_tjet':self.deltaR_tjet,\
                         't_target':self.t_index,'tbar_target':self.tbar_index,\
                         'parton':self.parton_def,'ljet_charge_max':self.ljet_charge_max}

        self.verbose = util.VERBOSE()
        self.verbose.level = self.verbose_level

        self.clear_output()

        return



    def clear_output(self):
        """Clear the dictionary output for each file. These are the features to save."""
        self.outdata = {'pt':[],'eta':[],\
                        'charge':[],\
                        'reco_m_ttbar':[],\
                        'truth_m_ttbar':[],\
                        'target':[],\
                        'n_tjets':[],\
                        'deltaQ':[],\
                        'mv2c10_100_charge':[],\
                        'mv2c10_85_charge':[],\
                        'mv2c10_77_charge':[],\
                        'mv2c10_70_charge':[],\
                        'mv2c10_60_charge':[]}
        for i in range(self.nsubjets):
            self.outdata['tjet_{0}_charge'.format(i)] = []
            self.outdata['tjet_{0}_mv2c10'.format(i)] = []
            self.outdata['tjet_{0}_pt'.format(i)]     = []
            self.outdata['tjet_{0}_eta'.format(i)]    = []
            self.outdata['tjet_{0}_phi'.format(i)]    = []

        return



    def execute(self):
        """Execute the code"""
        listOfFiles = util.file2list(self.listOfFiles)

        for file in listOfFiles:
            self.verbose.INFO("ROOT2JSON : Running over {0}".format(file))
            f   = ROOT.TFile.Open(file)
            nom = f.Get("nominal")
            tru = f.Get("truth")

            ## Physics data for NN
            self.clear_output()

            self.metadata['file']    = file
            self.outdata['metadata'] = self.metadata

            output_filename = file.split("/")[-1].split(".root")[0]
            self.output = self.outpath+"/"+output_filename+"_"+self.btag_wkpt+"_"+self.date
            self.verbose.INFO("ROOT2JSON : > Saving output to "+self.output)


            ## -- Event Loop -- ##
            tru_entry    = 0
            nEntries_tru = tru.GetEntries()
            nEntries_nom = nom.GetEntries()

            nEntries = 0
            if self.nEntries<0 or self.nEntries>nEntries_nom:
                nEntries = nEntries_nom
            else:
                nEntries = self.nEntries

            self.verbose.INFO("ROOT2JSON : Processing {0} entries".format(nEntries))
            for entry in xrange(nEntries):
                nom.GetEntry(entry)
                if not entry%10000: self.verbose.INFO("ROOT2JSON :   > Entry {0}".format(entry))
                if nEntries<10 and entry>0:
                    self.verbose.INFO("ROOT2JSON :   > Entry {0}".format(entry))

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

                # after FSR: truth tops
                truth_t    = None
                truth_tbar = None

                # only save information if the W decays hadronically (and exists)
                if fabs(tru.MC_Wdecay1_from_t_pdgId)<=5 and tru.MC_b_from_t_pt>1e-6 and tru.MC_W_from_t_pt>1e-6:
                    truth_t    = util.getTruthTopQuark('t',ttree=tru)
                if fabs(tru.MC_Wdecay1_from_tbar_pdgId)<=5 and tru.MC_b_from_tbar_pt>1e-6 and tru.MC_W_from_tbar_pt>1e-6:
                    truth_tbar = util.getTruthTopQuark('tbar',ttree=tru)

                if truth_t is None and truth_tbar is None:
                    self.verbose.DEBUG("ROOT2JSON : Truth tops are None ")
                    self.verbose.DEBUG("ROOT2JSON : Continuing to next event ")
                    continue

                truth_m_ttbar = tru.MC_ttbar_afterFSR_m

                ## -- Build the reco-level large-R jets -- ##
                ljet_t    = None
                ljet_tbar = None

                for l,lj in enumerate(nom.ljet_pt):
                    # Quality criteria (pT and charge)
                    if lj < self.ljet_pt_cut:
                        continue
                    if fabs( nom.ljet_eta[l]) > self.ljet_eta_cut:
                        continue
                    if fabs( nom.ljet_charge[l] ) > self.ljet_charge_max:
                        continue

                    # Build the object
                    ljet = self.build_largeRjet(nom,l)

                    # truth match
                    if util.deltaR( truth_t, ljet, self.deltaR_tru ):
                        ljet.target = self.t_index
                        ljet_t = deepcopy(ljet)
                    elif util.deltaR( truth_tbar, ljet, self.deltaR_tru ):
                        ljet.target = self.tbar_index
                        ljet_tbar = deepcopy(ljet)


                # Need at least one jet truth matched
                if ljet_t is None and ljet_tbar is None:
                    self.verbose.WARNING("ROOT2JSON : Reco tops not truth matched. Continuing.")
                    continue

                try:
                    reco_m_ttbar = (ljet_t.p4+ljet_tbar.p4).M()
                except AttributeError:
                    reco_m_ttbar = -1

                # Save top/anti-top information
                # - Top
                if ljet_t is not None:
                    ljet_t.reco_mttbar  = reco_m_ttbar
                    ljet_t.truth_mttbar = truth_m_ttbar
                    self.saveData(ljet_t)
                # - Anti-top
                if ljet_tbar is not None:
                    ljet_tbar.reco_mttbar  = reco_m_ttbar
                    ljet_tbar.truth_mttbar = truth_m_ttbar
                    self.saveData(ljet_tbar)


            self.verbose.INFO("ROOT2JSON : Writing data {0}.json".format(self.output))
            with open(self.output+'.json', 'w') as outfile:
                json.dump(self.outdata, outfile)


    def saveData(self,top):
        """Save information to dictionary"""
        self.outdata['reco_m_ttbar'].append(top.reco_mttbar)
        self.outdata['truth_m_ttbar'].append(top.truth_mttbar)
        self.outdata['charge'].append(top.charge)
        self.outdata['pt'].append(top.p4.Pt())
        self.outdata['eta'].append(top.p4.Eta())
        self.outdata['target'].append(top.target)
        self.outdata['n_tjets'].append(top.n_tjets)

        tmp_sub_tjets = []
        for i in range(self.nsubjets):
            self.outdata['tjet_{0}_charge'.format(i)].append(top.tjets[i].charge)
            self.outdata['tjet_{0}_mv2c10'.format(i)].append(top.tjets[i].mv2c10)
            self.outdata['tjet_{0}_pt'.format(i)].append(top.tjets[i].p4.Pt())
            self.outdata['tjet_{0}_eta'.format(i)].append(top.tjets[i].p4.Eta())
            self.outdata['tjet_{0}_phi'.format(i)].append(top.tjets[i].p4.Phi())
            tmp_sub_tjets.append( deepcopy(top.tjets[i]) )

        deltaQ = self.calculate_deltaQ(tmp_sub_tjets)
        self.outdata['deltaQ'].append( deltaQ )

        mv2c10_charges = self.calculate_mv2c10_charges(top.tjets)
        self.outdata['mv2c10_100_charge'].append(mv2c10_charges[100])
        self.outdata['mv2c10_85_charge'].append(mv2c10_charges[85])
        self.outdata['mv2c10_77_charge'].append(mv2c10_charges[77])
        self.outdata['mv2c10_70_charge'].append(mv2c10_charges[70])
        self.outdata['mv2c10_60_charge'].append(mv2c10_charges[60])

        ## print some stuff for debugging
        self.verbose.DEBUG("ROOT2JSON : {0}".format(top))
        self.verbose.DEBUG("ROOT2JSON :    - mv2c10 charges    = {0}".format(mv2c10_charges))
        self.verbose.DEBUG("ROOT2JSON :    - track jet index   = {0}".format([i.index for i in top.tjets]))
        self.verbose.DEBUG("ROOT2JSON :    - track jet mv2c10  = {0}".format([i.mv2c10 for i in top.tjets]))
        self.verbose.DEBUG("ROOT2JSON :    - track jet charges = {0}".format([i.charge for i in top.tjets]))

        return


    def build_largeRjet(self,ttree,index):
        """
        Construct the large-R jet object.
        
        @param ttree  TTree to read from (likely the 'nominal' one)
        @param index  Index of large-R jet in the vector
        """
        ljet = util.Ljet()
        ljet.index  = index
        ljet.charge = ttree.ljet_charge[index]
        ljet.p4.SetPtEtaPhiE(ttree.ljet_pt[index],ttree.ljet_eta[index],
                             ttree.ljet_phi[index],ttree.ljet_e[index])

        # ghost-associated track jets
        ghostTjets = [i for i in ttree.ljet_trackjet_index[index]]
        ghostTjets.sort() # put leading pT first (tjet_* are pT-sorted, 0=leading,1=subleading,etc.)

        numTjets = 0
        for i,ind in enumerate(ghostTjets):
            if ind<0: continue  # protection against non-existent ghost-associated track jets
            mv2c10s = util.getBtagging(ttree.tjet_mv2c10[ind])
            tjet = util.LorentzVector()
            tjet.index = ind
            tjet.p4.SetPtEtaPhiE( ttree.tjet_pt[ind],ttree.tjet_eta[ind],ttree.tjet_phi[ind],ttree.tjet_e[ind])
            tjet.mv2c10     = mv2c10s[0]
            tjet.mv2c10_bin = mv2c10s[1]
            tjet.charge = ttree.tjet_charge[ind]
            ljet.tjets.append( tjet )
            numTjets += 1

        if numTjets==0:
            return None
        ljet.n_tjets = numTjets

        # Add dummy track jets to fill out vector (need fixed size)
        # dummy tjets: pT=0,eta=ljet_eta,phi=ljet_phi,e=0;mv2c10=-1;charge=0
        if numTjets < self.nsubjets:
            dummy_tjet = util.LorentzVector()
            dummy_tjet.p4.SetPtEtaPhiE( 1e-6,ljet.p4.Eta(),ljet.p4.Phi(),1e-6)
            dummy_tjet.index  = -100
            dummy_tjet.mv2c10 = -1
            dummy_tjet.mv2c10_bin = 0
            dummy_tjet.charge = 0

            missing_tjets = self.nsubjets - numTjets
            for _ in range(missing_tjets):
                ljet.tjets.append( deepcopy(dummy_tjet) )

        if len(ljet.tjets)<self.nsubjets: 
            self.verbose.ERROR("ROOT2JSON : Not enough track jets: {0},{1},{2}".format(numTjets,len(ljet.tjets),ghostTjets))

        return ljet


    def calculate_deltaQ(self,sub_tjets):
        """Difference between charge of highest and lowest mv2c10"""
        sub_tjets.sort(key=lambda x: x.mv2c10, reverse=True)
        max_mv2c10_charge = sub_tjets[0].charge  # first index should be the largest mv2c10 value
        min_mv2c10_charge = sub_tjets[-1].charge # last index should be the lowest mv2c10 value
        return max_mv2c10_charge - min_mv2c10_charge

    def calculate_mv2c10_charges(self,track_jets):
        """Calculate the mv2c10 charges"""
        mv2c10s = {100:0,85:0,77:0,70:0,60:0,0:0}
        for i in track_jets:
            mv2c10s[i.mv2c10_bin] += i.charge

        return mv2c10s


## THE END ##

