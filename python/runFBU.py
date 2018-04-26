"""
Created:         17 April 2018
Last Updated:    26 April 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Perform Unfolding
Process histograms filled with nominal values and uncertainties

To run:
   python python/runFBU.py config.txt

where <config.txt> is the configuration file.
See 'config/fbuconfig.txt' as an example.
"""
import sys
import os
import json
import numpy as np

import util
import Analysis.uproot as uproot
from Analysis.fbu import PyFBU


if len(sys.argv)<2: 
    util.HELP()

config = util.read_config(sys.argv[1])

output_path = config['output_path']
listOfSystematics = util.file2list( config['listOfSystematics'] )

## Input information (histograms of delta|y| distribution)
hist  = config['histogram']           # name of histogram (e.g., 'deltay' or 'deltay_mttbar')
data  = config['data']                # real/pseudo data
bckg  = config['bckg']                # background predictions
ttbar = uproot.open(config['ttbar'])  # truth delta|y| distribution & response matrix

objsyst = config['systematics']       # detector-level systematics for signal (ttbar) and backgrounds


unf.runUnfolding()

unf.save()

## THE END ##
