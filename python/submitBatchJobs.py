"""
Created:        21 September 2016
Last Updated:   28 January   2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----
Script for using the BatchSubmission() class to submit batch jobs
This is setup for running at the LPC.
If you are running in a different environment, you will need
to modify this script accordingly (mostly concerning the output storage).
This is also designed to submit one batch job at a time.

To run:
$ python python/runBatch.py batchConfig.txt
where <batchConfig.txt> contains configuration options
"""
import os
import sys
import util
import time
import commands
from batch.batchSubmission import BatchSubmission

# setup config
config = sys.argv[1]
cfg    = util.read_config(config)

# reset boolean values
cfg['test']    = util.str2bool( cfg['test'] )
cfg['submit']  = util.str2bool( cfg['submit'] )
cfg['verbose'] = util.str2bool( cfg['verbose'] )

# set 'global' options
date  = time.strftime("%d%b%Y") if cfg['date']=='today' else cfg['date']
files = util.file2list(cfg["files"])                                      # set of files to process (primary dataset names)
eos_path_base = cfg['eos_path_base'].replace("${DATE}",date)              # '/store/user/demarley/susy/semi-resolved-tagging/'+date+'/{0}/'
cmaConfig = util.read_config( cfg['config'] )                             # get the selection from the cmaConfig file

# loop over datasets to process, e.g., ttbar or qcd
for file in files:

    input_files = 'config/samples/{0}.txt'.format(file)  # individual root files to process

    batch = BatchSubmission()

    ## -- set arguments for the batch jobs -- ##
    batch.username   = cfg['username']
    batch.executable = cfg['executable']       # executable to run for the program
    batch.mode       = cfg['mode']             # ['condor','lxbatch','qsub']   == kind of batch job
    batch.test       = cfg['test']             # just submit one job if testing
    batch.submit     = cfg['submit']           # do/don't submit jobs, just testing code
    batch.verbose    = cfg['verbose']          # lots of print statements
    batch.queue      = cfg['queue']            # batch system queue ('condor_submit')
    batch.config     = cfg['config']           # configuration file to use
    batch.fileDir    = file                    # unique directory for each file
    batch.file       = input_files             # individual root files to process

    ## Setup output
    eos_path       = eos_path_base.format(file)
    eos_path_full  = eos_path+"/"+cmaConfig['selection']  # files saved in a directory named after the selection
                                                          # as defined in run.cxx/runML.cxx -- update if changed
    batch.eos_path = eos_path
    batch.eos_tarball_path = eos_path_base.replace('/{0}','')

    # create eos directory, if it doesn't exist
    os_err = commands.getoutput("eos root://cmseos.fnal.gov/ mkdir -p {0}".format(eos_path_full))
    if os_err:
        print "RUNBATCH :: INFO : Attemp to make directory {0}".format(eos_path_full)
        print "RUNBATCH :: INFO : Message = {0}".format(os_err)
    else:
        print "RUNBATCH :: INFO : Created directory {0}".format(eos_path_full)

    # define the directory to write the output
    if cfg['output_dir']=='eos':
        batch.output_dir = 'root://cmseos.fnal.gov/'+eos_path
    else:
        batch.output_dir = cfg['output_dir']

    ## Submit
    batch.execute()

## THE END ##

