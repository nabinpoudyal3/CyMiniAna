"""
1 July 2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Execute the root2json class to convert
ROOT data into JSON format for neural network
training/testing/etc. in python.
"""
from info import VERBOSE
from root2json import Root2json

vb = VERBOSE()
vb.level = "INFO"

vb.INFO("RUN > Set up the root2json object")
r2j = Root2json()

## Define properties (can put this into config file later, if wanted)
vb.INFO("RUN > Define properties for convertin ROOT to JSON")
r2j.verbose_level = "INFO"  # Setup verbose output
r2j.outpath     = "./"      # where to store output
r2j.listOfFiles = "share/listOfFiles_testNN.txt" # ROOT files to process
r2j.nEntries    = 5

# Properties for large-R jets and such
r2j.ljet_charge_max = 5.
r2j.ljet_pt_cut  = 300000.
r2j.ljet_eta_cut = 2.
r2j.tjet_pt_cut  = 10000.
r2j.deltaR_tru   = 0.75
r2j.deltaR_tjet  = 0.8        # ljet R = 1.0; tjet R = 0.2
r2j.t_index      = 1          # +2/3 charge
r2j.tbar_index   = 0          # -2/3 charge
r2j.nsubjets     = 3          # number of ghost-associated track jets to save
r2j.parton_def   = 'afterFSR' # truth parton definition
r2j.success      = '\x01'     # if something is b-tagged (type char)
r2j.btag_wkpt    = "77"       # this isn't necessarily needed anymore for actual selection

# Setup and run the code
vb.INFO("RUN > Initialize")
r2j.initialize()

vb.INFO("RUN > Execute")
r2j.execute()

vb.INFO("RUN > Finished")

## THE END ##

