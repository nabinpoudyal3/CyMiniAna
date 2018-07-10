"""
Created:         17 April 2018
Last Updated:    30 April 2018

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
import uproot
import fbu
import util


if len(sys.argv)<2: 
    util.HELP()


print
print " ------------------------------ "
print " *    Unfolding with PyFBU    * "
print " ------------------------------ "
print


config = util.read_config(sys.argv[1])

date   = strftime("%d%b", localtime())
vb     = util.VERBOSE()

## Set configuration options ##
if len(sys.argv)<2:
    vb.HELP()
    sys.exit(1)

vb.level = config['verbose_level']
vb.initialize()


## Set output directory
output_dir  = "nMCMC{0}_".format(config.nMCMC)
output_dir += "nThin{0}_".format(config.nThin)

hep_data_name = config.hep_data.split('/')[-1].split('.')[0]
listOfSystematics = util.file2list( config['listOfSystematics'] )

## Setup Deep Learning class
fbu = Unfolding()

fbu.nMCMC = config.nMCMC
fbu.nThin = config.nThin
fbu.monitoring = config.monitoring
fbu.output_dir = output

#fbu.stat_only = config['stat_only']
#data files

if not os.path.isdir(output_dir):
    vb.INFO("RUN : '{0}' does not exist ".format(output))
    vb.INFO("RUN :       Creating the directory. ")
    os.system( 'mkdir -p {0}'.format(output_dir) )
else:
    vb.INFO("RUN :  Saving output to {0}".format(output_dir))


## Setup
fbu.initialize()
fbu.execute()

## THE END ##
