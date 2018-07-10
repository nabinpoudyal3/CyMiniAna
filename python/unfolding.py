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
np.random.seed(2018)


class Unfolding(object):
    """Unfolding class"""
    def __init__(self):
        self.date = datetime.date.today().strftime('%d%b%Y')

        ## Handling unfolding objects and data -- set in the class
        self.output_dir  = ""           # directory for storing outputs
        self.plotter     = UnfoldingPlotter()  # class for plotting relevant unfolding information

        self.variables   = []           # variables to load from the dataframe ('deltay','mtt',etc.)
        self.backgrounds = []           # names of background samples (e.g., w+jets, z+jets, etc.)
        self.stat_only   = False        # only consider the statistical uncertainty
        self.xsections   = None         # normalization uncertainties, e.g., util.read_config('xsection.txt')
        self.objsysts    = None         # systematic uncertainties for each detector-related object
        self.regularization = None      # regularization to use


    def initialize(self):   #,config):
        """Initialize a few parameters after they've been set by user"""
        self.msg_svc       = util.VERBOSE()
        self.msg_svc.level = self.verbose_level
        self.msg_svc.initialize()

        ## -- Plotting framework
        self.plotter.output_dir   = self.output_dir
        self.plotter.image_format = 'pdf'           # must use .pdf at the LPC

        return


    def execute(self):
        """
        Perform unfolding
            Normalization uncertainties -- keep stored in a text file for easy, global access
        """
        pass


    def load_hep_data(self):
        """
        Load the physics data (histograms) for unfolding
            ROOT I/O
        """
        pass


    def diagnostics(self,pre=False,post=False):
        """Diagnostic tests of the Unfolding"""
        pass


## THE END ##
