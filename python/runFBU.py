"""
Created:         17 April 2018
Last Updated:    17 April 2018

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


## FBU
pyfbu = PyFBU()

pyfbu.nMCMC      = 1000000    # Sampling
pyfbu.monitoring = False      # Monitoring
pyfbu.nThin      = 1

truth = ttbar['truth_deltay']
pyfbu.lower      = [t/5 for t in truth]      # lower bound of hyperbox for sampling
pyfbu.upper      = [t*2 for t in truth]      # upper bound of hyperbox for sampling

pyfbu.data       = uproot.open(data)[hist]
pyfbu.response   = ttbar['h_resmat']
pyfbu.background = uproot.open(bckg)[hist]

## Normalization uncertainties
## > pyfbu.backgroundsyst = {'bckg':0.} # stat only
pyfbu.backgroundsyst = {
    'qcd':0.5,
}

## load systematics
exp_systs = uproot.open(objsyst)
systematics_dict = #
bckgkeys = pyfbu.backgroundsyst.keys()

signalsystdict = dict([(str(k),v) for k,v in systematics_dict['signal'].iteritems() if k in listOfSystematics])
signalsystdict['LUMI'] = [0. for _ in pyfbu.data]

bckgsystdict   = dict([(str(k),v) for k,v in systematics_dict['background'].iteritems() if k in listOfSystematics])
bckgsystdict['LUMI'] = dict(zip(bckgkeys,[[0. if (mcStat or 'wjets' in k or 'qcd' in k) else 0.028 for _ in pyfbu.data] for k in bckgkeys]))

pyfbu.objsyst = {'signal':signalsystdict,'background':bckgsystdict}

## Run the algorithm
pyfbu.run()



## Save output posteriors
bckgkeys = pyfbu.backgroundsyst.keys()
unfbins = pyfbu.trace

np.save(outdir+'unfolded',unfbins)        ## save unfolded result
for nuis in bckgkeys+listOfSystematics:
    trace = pyfbu.nuisancestrace[nuis]
    np.save(outdir+nuis,trace)            ## save posteriors

## THE END ##
