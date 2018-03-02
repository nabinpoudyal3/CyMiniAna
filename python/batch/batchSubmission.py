"""
Created:        21 September 2016
Last Updated:   16 July      2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109
-----
Base class for submitting batch jobs.

Setup for lxplus and condor batch systems (only ones the author has used).
"""
import os
import sys
import commands


class BatchSubmission(object):
    """Class for submitting batch jobs in different environments"""
    def __init__(self):
        """Initialize some variables for the class"""
        self.username   = 'dmarley'
        self.executable = 'run'
        self.baseDir    = os.environ['PWD']+'/'
        self.mode       = 'lxbatch' # 'condor','lxbatch'
        self.test       = False     # just submit one job if testing
        self.submit     = True      # don't submit jobs, just testing code
        self.verbose    = False     # print output to screen
        self.queue      = ''        # lxbatch:'8nh', condor:'condor_submit'

        self.config     = 'config/cmaConfig.txt' # CyMiniAna configuration file
        self.file       = 'share/miniSL_ntuples.txt'
        self.mc_files   = 'share/miniSL_ntuples.txt' # where all MC files you want to run over are listed

        self.cfg_filename = "batch/batchConfig_{0}.txt" # CyMiniAna config file per job
        self.filename     = ''

        self.split_by_syst = False
        self.fileEnds      = {'lxbatch':'sh','condor':'condor','slurm':'slrm'}

        self.writeBatchScriptTemplates()


    def initialize(self):
        """Setup a few variables based on class attributes"""
        self.unique_id_base = [self.tmp_ntuple]
        self.unique_id_name = "_".join(self.unique_id_base)
        self.unique_id_path = "/".join(self.unique_id_base)
        self.unique_id_batch_path = "{0}/batch/{1}".format(self.baseDir,self.unique_id_path)

        self.config_options = {}
        config = open(self.config,'r').readlines()
        for i in config:
            j = i.rstrip('\n').split(' ')
            key,item = j
            if item.startswith('config/'):
                item = self.baseDir+item
            self.config_options[key] = item

        return



    def writeConfiguration(self):
        """Write the CyMiniAnaAC configuration file"""
        config_text = ""
        for key in self.config_options:
            value = self.config_options[key]
            if key == "inputfile":
                config_text += "inputfile {0}\n".format(self.filename)
            elif key == "sumWeightsFiles":
                config_text += "sumWeightsFiles {0}\n".format(self.filename) # not necessary in current incarnation
            else:
                config_text += "{0} {1}\n".format(key,value)

        cfg_file = open(self.cfg_filename,"w")
        cfg_file.write(config_text)
        cfg_file.close()

        return


    def writeFile(self,minisl_ntuple):
        """
           Write the name of the root file to a text file.
           - CyMiniAnaAC uses a text file to read the root files that it processes.
        """
        tmp_file = open(self.filename,"w")
        tmp_file.write(minisl_ntuple)
        tmp_file.close()

        return


    def writeBatchScript(self):
        """Write the batch submission script."""
        print " BATCH SUBMISSION : Write batch scripts "

        lxbatch_template = self.lxbatch_template % {'baseDir':self.baseDir,
            'cfg_filename':self.cfg_filename,'executable':self.executable,
            'ntuple':self.tmp_ntuple,'unique_id':self.unique_id_path}

        slurm_template = self.slurm_template % {'baseDir':self.baseDir,'unique_id':self.unique_id_batch_path,
            'cfg_filename':self.cfg_filename,'executable':self.executable,'jobID':'%j'}

        # condor file executes a .sh file in the batch/ directory
        condorExecFileName = "{0}/run_condor.sh".format(self.unique_id_batch_path)
        condor_executable  = self.condor_exec_template % {'baseDir':self.baseDir,
            'executable':self.executable,
            'ntuple':self.tmp_ntuple,'cfg_filename':self.cfg_filename,
            'unique_id':self.unique_id_path,'username':self.username}
        condor_template = self.condor_template % {'baseDir':self.baseDir,
            'condorExec':condorExecFileName,'unique_id':self.unique_id_path}

        fileEnd  = ''
        template = None

        if self.mode=='condor':
            template = condor_template
            # need executable that sets up the environment properly
            condorExecFile     = open(condorExecFileName,'w')
            condorExecFile.write(condor_executable)
            commands.getoutput("chmod +x {0}".format(condorExecFileName))
        elif self.mode=='slurm':
            template = slurm_template
        else:
            template = lxbatch_template

        batchFile = open(self.batchFileName,'w')
        batchFile.write(template)

        return



    def execute(self):
        """Execute the script"""
        # -- Print the configuration
        if self.verbose:
            self.printConfiguration()


        # -- Make some directories, if needed
        if not os.path.exists('batch'):
            os.makedirs('batch')


        # -- Submit jobs
        if not self.file.startswith("/"):
            self.file = self.baseDir+"/"+self.file

        miniSL_ntuples = open(self.file,'r').readlines() # files to run over

        for l,line in enumerate(miniSL_ntuples):
            minisl_ntuple = line.split()[0]

            # -- simplify the filename for reference later
            self.tmp_ntuple = minisl_ntuple
            if minisl_ntuple.split('/')[-1].startswith('user.'):
                # check if submitting jobs straight from analysisTop
                cnd_name   = minisl_ntuple.split('/')[-1].split('.')
                self.tmp_ntuple = cnd_name[2]+cnd_name[3]
            else:
                # normal file from "_pre" selection or something similar
                self.tmp_ntuple = minisl_ntuple.split('/')[-1].split('.')[0]


            # -- Setup the unique ID for making unique files/configurations/etc.
            self.initialize()
            if not os.path.exists(self.unique_id_batch_path):
                os.makedirs(self.unique_id_batch_path)

            output_directory = self.config_options['output_path']
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)


            ## Put the single root file into a text file
            self.filename = "batch/{0}/listOfFiles.txt".format(self.unique_id_path)
            self.writeFile(minisl_ntuple)

            if self.verbose:
                print " "
                print " BATCH SUBMISSION : {0}".format(minisl_ntuple)
                print " BATCH SUBMISSION : >> {0}".format(self.tmp_ntuple)


            ## Add functionality to run 1 systematic at a time (save resources; merge later)
            if self.split_by_syst:
                syst_file = []
                try:
                    syst_file = open(self.systematics,"r").readlines()
                except TypeError:
                    syst_file = self.systematics # maybe this is a list already?

                for systematic in syst_file:
                    self.systematic = systematic.rstrip('\n')

                    self.unique_id_name = "_".join(self.unique_id_base)+"_{0}".format(self.systematic)
                    self.unique_id_path = "/".join(self.unique_id_base)+"/{0}".format(self.systematic)

                    self.unique_id_batch_path = "{0}/batch/{1}".format(self.baseDir,self.unique_id_path)

                    if not os.path.exists(self.unique_id_batch_path):
                        os.makedirs(self.unique_id_batch_path)

                    # write the systematic to a text file (to read)
                    self.syst_filename = "batch/{0}/treenames.txt".format(self.unique_id_path)
                    self.syst_file = open( self.syst_filename, 'w')
                    self.syst_file.write(systematic)
                    self.syst_file.close()

                    self.submit_job()
            else:
                self.submit_job()

        return


    def submit_job(self):
        """Submit the job"""
        # -- Write the CyMiniAnaAC configuration file
        self.cfg_filename = "batch/{0}/cmaConfig.txt".format(self.unique_id_path)
        self.writeConfiguration()
        print " BATCH SUBMISSION : Config filename = {0}".format(self.cfg_filename)

        # -- Write the script executed by the batch system
        fileEnd = self.fileEnds[self.mode]
        self.batchFileName = "{0}/run_batch.{1}".format(self.unique_id_batch_path,fileEnd)
        self.writeBatchScript()
        print " BATCH SUBMISSION : Batch filename  = {0}".format(self.batchFileName)


        # only submit the jobs if you want to (True by default)
        if self.submit:
            if self.mode=='lxbatch':
                commands.getoutput("chmod +x {0}".format(self.batchFileName))
                batch_command_str = "bsub -q {0} {1}".format(self.lxbatch_queue,self.batchFileName)
            elif self.mode=='slurm':
                commands.getoutput("chmod +x {0}".format(self.batchFileName))
                batch_command_str = "sbatch {0}".format(self.batchFileName)
            else:
                batch_command_str = "{0} {1}".format(self.condor_queue,self.batchFileName)

            batch_command = commands.getoutput( batch_command_str )

            if self.verbose:
                print " BATCH SUBMISSION : {0}".format(batch_command_str)
                print " BATCH SUBMISSION : {0}".format(batch_command)

        # if testing job submission, submit only one and exit
        if self.test:
            print " BATCH SUBMISSION : Test complete. "
            print " BATCH SUBMISSION : Exiting. "
            sys.exit(1)

        return



    def writeBatchScriptTemplates(self):
        """Write templates for the batch scripts"""
        self.lxbatch_template = """#!/bin/sh
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'

CURRENTDIR=${PWD}
NEWDIR=%(unique_id)s
mkdir -p ${NEWDIR}
cd ${NEWDIR}

cp -r %(baseDir)s/batch/ .
cp -r %(baseDir)s/share/ .
cp -r %(baseDir)s/util/ .
cp -r %(baseDir)s/Root/ .
cp -r %(baseDir)s/CyMiniAna/ .
cp -r %(baseDir)s/lwtnn/ .
cp -r %(baseDir)s/boost/ .
cp -r %(baseDir)s/python/ .
cp %(baseDir)s/setup.sh .
cp %(baseDir)s/Makefile* .

source setup.sh
make -f Makefile

./%(executable)s %(cfg_filename)s
# python python/batch2eos.py

echo " Job finished "
cd ${CURRENTDIR}
echo " Delete contents "
rm -rf %(ntuple)s
"""

        self.condor_exec_template = """#!/bin/bash
export CERN_USER=%(username)s
export RUCIO_ACCOUNT=$CERN_USER
export CYMINIANADIR=%(baseDir)s

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'

NEWDIR=/tmp/demarley/%(unique_id)s
mkdir -p ${NEWDIR}
cd ${NEWDIR}
echo " Currently in "
echo ${PWD}

cp -r %(baseDir)s/batch/ .
cp -r %(baseDir)s/share/ .
cp -r %(baseDir)s/util/ .
cp -r %(baseDir)s/Root/ .
cp -r %(baseDir)s/CyMiniAna/ .
cp -r %(baseDir)s/lwtnn/ .
cp -r %(baseDir)s/boost/ .
cp -r %(baseDir)s/python/ .
cp %(baseDir)s/setup.sh .
cp %(baseDir)s/Makefile* .

source setup.sh
make -f Makefile

## Environment setup, now execute the script
echo " > %(executable)s "
./%(executable)s %(cfg_filename)s

echo " - Done - "
cd /tmp/demarley/
echo " Delete contents "
rm -rf ${NEWDIR}
"""

        self.condor_template = """
Universe   = vanilla
Executable = %(condorExec)s

Log        = batch/%(unique_id)s/info.log
Error	   = batch/%(unique_id)s/info.error
Output     = batch/%(unique_id)s/info.out

Queue
"""


        self.slurm_template = """#!/bin/bash
#SBATCH -J cma-%(jobID)s
#SBATCH -p background
#SBATCH -n1
#SBATCH -o %(unique_id)s/cma-%(jobID)s.out
#SBATCH -e %(unique_id)s/cma-%(jobID)s.err

echo " > Starting SLURM job at `date` on `hostname`"
echo "   SLURM_JOBID=$SLURM_JOBID"

# Run the vmtest application
echo " > Setup the environment "
source %(baseDir)s/setup.sh

echo " > Run the program "
%(executable)s %(cfg_filename)s

echo " > Ended at `date` on `hostname`"
exit 0
"""

        return


    def printConfiguration(self):
        """Print the configuration of arguments for batch submission"""
        print ""
        print " BATCH SUBMISSION : Configuration "
        print " -------------------------------- "
        print ""

        attributes = ['username','executable','config','baseDir',
        'mode','test','submit','verbose','queue']

        for attr in attributes:
            print " %-*s : %s" % (16,attr,getattr(self,attr))

        print ""
        print " -------------------------------- "
        print ""

        return



## THE END
