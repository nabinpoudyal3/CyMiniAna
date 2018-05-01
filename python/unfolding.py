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

import matplotlib
matplotlib.use('PDF')   # png not supported at LPC; do this before anything else tries to set the backend
import uproot
import numpy as np
import pandas as pd

import util
from unfoldingPlotter import unfoldingPlotter


# fix random seed for reproducibility
seed = 2018
np.random.seed(seed)


class Unfolding(object):
    """Unfolding class"""
    def __init__(self):
        self.date = datetime.date.today().strftime('%d%b%Y')

        ## Handling NN objects and data -- set in the class
        self.df          = None         # dataframe containing physics information
        self.output_dir  = ""           # directory for storing outputs

        self.variables   = []           # variables to load from the dataframe ('deltay','mtt',etc.)
        self.backgrounds = []           # names of background samples
        self.stat_only   = False        # only consider the statistical uncertainty
        self.xsections   = None         # normalization uncertainties, e.g., util.read_config('xsection.txt')
        self.exp_systs   = None         # systematic uncertainties for each detector-related object

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

        self.pyfbu.lower = [t/5 for t in self.truth_data]  # lower bound of hyperbox for sampling
        self.pyfbu.upper = [t*2 for t in self.truth_data]  # upper bound of hyperbox for sampling

        self.pyfbu.data       = self.df['data']
        self.pyfbu.response   = self.df['h_resmat']
        self.pyfbu.background = self.df['bckg']

        ## -- Plotting framework
        print " >> Store output in ",self.output_dir
        self.plotter = UnfoldingPlotter()  # class for plotting relevant NN information
        self.plotter.output_dir   = self.output_dir
        self.plotter.image_format = 'pdf'           # must use .pdf at the LPC

        return


    def execute(self):
        """Perform unfolding"""

        ## Normalization uncertainties -- keep stored in a text file for easy, global access

        if self.stat_only: 
            self.pyfbu.backgroundsyst = {'bckg':0.} # stat only
            self.pyfbu.objsyst = {'signal':{},'background':{}}
        else: 
            self.pyfbu.backgroundsyst = dict( (k,self.xsection[k]) for k in self.backgrounds )   # {'qcd':xsection['qcd']}

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



    def load_hep_data(self,variables2plot=[]):
        """
        Load the physics data (flat ntuple) for NN using uproot
        Convert to DataFrame for easier slicing 

        @param variables2plot    If there are extra variables to plot, include them here
        """
        file = uproot.open(self.hep_data)
        data = file[self.treename]
        self.df = data.pandas.df( self.variables+variables2plot )

        self.truth_data = self.df['truth_deltay']
        self.load_resmat()

        self.metadata = file['metadata']   # names of samples, target values, etc.

        return


    def load_resmat(self):
        """Load existing model to make plots or predictions"""
        self.resmat = self.df['resmat']
        return


    def build_resmat(self):
        """Build the response matrix"""
        # save truth vs reco histogram & divide by truth (acceptance)
        return


    def save_resmat(self):
        """Save the model for use later"""
        return


    def save_posteriors(self):
        """
        Save the features to a json file to load via lwtnn later
        Hard-coded scale & offset; must change later if necessary
        """
        bckgkeys = pyfbu.backgroundsyst.keys()
        unfbins = pyfbu.trace

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
            #prefit data/mc
            #draw response matrix

        # post unfolding
        if postUnfolding:
            self.msg_svc.INFO("FBU : -- post-unfolding")
            #pulls
            #postfits

        return


## THE END ##
