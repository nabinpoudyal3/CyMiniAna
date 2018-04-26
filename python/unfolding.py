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
import util
import datetime

import matplotlib
matplotlib.use('PDF')   # png not supported at LPC; do this before anything else tries to set the backend
import uproot
import numpy as np
import pandas as pd

from unfoldingPlotter import unfoldingPlotter


# fix random seed for reproducibility
seed = 2018
np.random.seed(seed)


class Unfolding(object):
    """Unfolding class"""
    def __init__(self):
        self.date = datetime.date.today().strftime('%d%b%Y')

        ## Handling NN objects and data -- set in the class
        self.df  = None          # dataframe containing physics information

        self.pyfbu = PyFBU()


    def initialize(self):   #,config):
        """Initialize a few parameters after they've been set by user"""
        self.msg_svc       = util.VERBOSE()
        self.msg_svc.level = self.verbose_level
        self.msg_svc.initialize()


        ## -- Plotting framework
        print " >> Store output in ",self.output_dir
        self.plotter = UnfoldingPlotter()  # class for plotting relevant NN information
        self.plotter.output_dir   = self.output_dir
        self.plotter.image_format = 'pdf'


        ## -- Adjust unfolding parameters
        self.pyfbu.nMCMC      = config['nMCMC']                     # default = 1000000    # Sampling
        self.pyfbu.monitoring = util.str2bool( config['monitor'] )  # default = False      # Monitoring
        self.pyfbu.nThin      = int( config['nThin'] )              # default = 1

        self.load_hep_data()
        truth = ttbar['truth_deltay']
        self.pyfbu.lower = [t/5 for t in truth]      # lower bound of hyperbox for sampling
        self.pyfbu.upper = [t*2 for t in truth]      # upper bound of hyperbox for sampling

        self.pyfbu.data = uproot.open(data)[hist]
        self.pyfbu.response   = ttbar['h_resmat']
        self.pyfbu.background = uproot.open(bckg)[hist]

        return


    def runUnfolding(self):
        """Perform unfolding"""

## Normalization uncertainties -- keep stored in a text file for easy, global access
## > pyfbu.backgroundsyst = {'bckg':0.} # stat only
xsection = util.read_config('xsection.txt')
self.pyfbu.backgroundsyst = {
    'qcd':xsection['qcd'],
}

## load systematics from ROOT file
## > Update LUMI uncertainty depending on data-driven / floating backgrounds
exp_systs = uproot.open(objsyst)
systematics_dict = {'signal':{},'background':{}}#

signalsystdict = dict([(str(k),v) for k,v in systematics_dict['signal'].iteritems() if k in listOfSystematics])
signalsystdict['LUMI'] = [xsection['LUMI'] for _ in pyfbu.data]

bckgsystdict   = systematics_dict['background']
bckgsystdict['LUMI'] = dict(zip(bckgkeys,[[0. if k in data_driven else xsection['LUMI'] for _ in pyfbu.data] for k in bckgkeys]))

self.pyfbu.objsyst = {'signal':signalsystdict,'background':bckgsystdict}

## Run the algorithm
self.pyfbu.run()

        return



    def load_hep_data(self,variables2plot=[]):
        """
        Load the physics data (flat ntuple) for NN using uproot
        Convert to DataFrame for easier slicing 

        @param variables2plot    If there are extra variables to plot, 
                                 that aren't features of the NN, include them here
        """
        file = uproot.open(self.hep_data)
        data = file[self.treename]
        dataframe = data.pandas.df( self.features+['target']+variables2plot )

        self.metadata = file['metadata']   # names of samples, target values, etc.

        return


    def load_resmat(self):
        """Load existing model to make plots or predictions"""
        return


    def build_resmat(self):
        """Build the response matrix"""
        return


    def save_resmat(self):
        """Save the model for use later"""
        output = self.output_dir+'/'+self.model_name

        return


    def save_posteriors(self):
        """
        Save the features to a json file to load via lwtnn later
        Hard-coded scale & offset; must change later if necessary
        """
## Save output posteriors
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
            #prefit
            #response matrix

        # post unfolding
        if postUnfolding:
            self.msg_svc.INFO("DL : -- post-unfolding")
            #pulls
            #postfits

        return


## THE END ##
