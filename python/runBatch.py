"""
Created:        21 September 2016
Last Updated:   16 July      2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109
-----
Script for using the BatchSubmission() class to submit batch jobs
"""
import sys
from batchSubmission import BatchSubmission

batch = BatchSubmission()


## -- set arguments for the batch jobs -- ##
batch.username         = 'dmarley'
batch.executable       = 'run'    # executable to run for the program
batch.mode             = 'slurm'  # ['condor','lxbatch','qsub']   == kind of batch job
batch.test             = False    # just submit one job if testing
batch.submit           = True     # do/don't submit jobs, just testing code
batch.verbose          = True     # lots of print statements
batch.queue            = ''

## Configurations
batch.config   = 'config/cmaConfig.txt'    # configuration file to use
batch.file     = 'config/listOfDataMCSamples.txt'
batch.mc_files = 'config/listOfDataMCSamples.txt'  # not necessary for now

## Submit
batch.execute()

# END #
