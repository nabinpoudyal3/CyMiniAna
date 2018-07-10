"""
Created:        26 April 2018
Last Updated:   26 April 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Base class for performing unfolding
"""
import json
import datetime

import uproot
import numpy as np
import pandas as pd

import util
from unfoldingPlotter import UnfoldingPlotter


# fix random seed for reproducibility
seed = 2018
np.random.seed(seed)


class Unfolding(object):
    """Unfolding class"""
    def __init__(self):
        self.date = datetime.date.today().strftime('%d%b%Y')

        ## Handling NN objects and data -- set in the class
        self.output_dir  = ""           # directory for storing outputs

        self.variables   = []           # variables to load from the dataframe ('deltay','mtt',etc.)
        self.backgrounds = []           # names of background samples
        self.stat_only   = False        # only consider the statistical uncertainty
        self.xsections   = None         # normalization uncertainties, e.g., util.read_config('xsection.txt')
        self.objsysts    = None         # systematic uncertainties for each detector-related object

        ## -- Adjust unfolding parameters
        self.nMCMC       = 1000000
        self.monitoring  = False
        self.nThin       = 1


    def initialize(self):   #,config):
        """Initialize a few parameters after they've been set by user"""
        self.msg_svc       = util.VERBOSE()
        self.msg_svc.level = self.verbose_level
        self.msg_svc.initialize()

        self.load_hep_data()            # Load the physics data

        ## PyFBU
        self.pyfbu = PyFBU()            # Python framework for Fully Bayesian Unfolding

        self.pyfbu.nMCMC      = self.nMCMC                 # Sampling
        self.pyfbu.monitoring = self.monitoring            # Monitoring
        self.pyfbu.nThin      = self.nThin

        self.pyfbu.lower = [t/5 for t in self.truth]       # lower bound of hyperbox for sampling
        self.pyfbu.upper = [t*2 for t in self.truth]       # upper bound of hyperbox for sampling

        self.pyfbu.data       = self.data.numpy[0]
        self.pyfbu.response   = self.resmat
        self.pyfbu.background = self.bckgs
# backgrounds = {'bckg1':bckg1[0].tolist(),
#                'bckg2':bckg2[0].tolist()} # bin content for each background
# backgrounds_systs = {'bckg1':0.5,     # 50% normalization uncertainty
#                      'bckg2':0.04}    #  4% normalization uncertainty
# myfbu.objsyst = {
#     'signal':{
#             'syst1':[0.,0.03,0.02,0.01],
#             'syst2':[0.,0.01,0.00,0.01]},
#     'background':{
#             'syst1':{'bckg1':[0.,0.,0.,0.],
#                      'bckg2':[0.1,0.1,0.1,0.1]},
#             'syst2':{'bckg1':[0.,0.01,0.01,0.],
#                      'bckg2':[0.,0.,0.,0.]}
#     }
# }

        ## -- Plotting framework
        print " >> Store output in ",self.output_dir
        self.plotter = UnfoldingPlotter()  # class for plotting relevant NN information
        self.plotter.output_dir   = self.output_dir
        self.plotter.image_format = 'pdf'           # must use .pdf at the LPC

        return


    def execute(self):
        """
        Perform unfolding
           Normalization uncertainties -- keep stored in a text file for easy, global access
        """
        if self.stat_only: 
            self.pyfbu.backgroundsyst = {'bckg':0.} # stat only
            self.pyfbu.objsyst = {'signal':{},'background':{}}
        else: 
            self.pyfbu.backgroundsyst = dict( (k,self.xsection[k]) for k in self.backgrounds )   # {'qcd':self.xsection['qcd']}

            ## load systematics from ROOT file
            ## > Update LUMI uncertainty depending on data-driven / floating backgrounds
            systematics_dict = {'signal':{},'background':{}}

            signalsystdict = dict([(str(k),v) for k,v in systematics_dict['signal'].iteritems() if k in self.systematics])
            signalsystdict['LUMI'] = [self.xsection['LUMI'] for _ in pyfbu.data]

            bckgsystdict   = systematics_dict['background']
            bckgsystdict['LUMI'] = dict(zip(self.backgrounds,[[0. if k in self.data_driven else self.xsection['LUMI'] for _ in pyfbu.data] for k in self.backgrounds]))

            self.pyfbu.objsyst = {'signal':signalsystdict,'background':bckgsystdict}

        ## Run the algorithm
        self.pyfbu.run()

        return


    def load_hep_data(self):
        """Load the physics data (histograms) for unfolding"""
        file = uproot.open(self.hep_data)

        self.data   = file['data']
        self.truth  = file['truth_deltay']
        self.resmat = file['resmat']

        self.backgrounds = []
        self.backgrounds_syst = []

        self.systs = {}

        return


    def build_resmat(self):
        """
        Build the response matrix
           save truth vs reco histogram & divide by truth (acceptance)
        """
        pass


    def save_resmat(self):
        """Save the resmat, if it was made here, for use later"""
        pass


    def save_posteriors(self):
        """Save posteriors to access later"""
        bckgkeys = pyfbu.backgroundsyst.keys()
        unfbins  = pyfbu.trace

        np.save(outdir+'unfolded',unfbins)        ## save unfolded result
        for nuis in bckgkeys+listOfSystematics:
            trace = pyfbu.nuisancestrace[nuis]
            np.save(outdir+nuis,trace)            ## save posteriors

        return


    def diagnostics(self,preTraining=False,postTraining=False):
        """Diagnostic tests of the NN"""
        self.msg_svc.INFO("DL : Diagnostics")

        # Diagnostics before the unfolding
        if preUnfolding:
            self.msg_svc.INFO("FBU : -- pre-unfolding")
            self.plotter.prefit_yields()     # prefit data/mc
            self.plotter.response_matrix()   # draw response matrix

        # post unfolding
        if postUnfolding:
            self.msg_svc.INFO("FBU : -- post-unfolding")
            self.plotter.trace()             # trace for each unfolded bin
            self.plotter.nuisances()         # prior and posteriors for nuisances
            self.plotter.pulls()             # pulls
            self.plotter.postfit_yields()    # postfits

        return


## THE END ##
