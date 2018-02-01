14 October 2017  
Dan Marley  


# CyMiniAna

#### C++/Python Analysis Framework

This class is derived from the [CyMiniAna](https://gitlab.cern.ch/dmarley/CyMiniAna) base class.  
Dependencies (available in CMSSW):

    - ROOT6
    - python2.7 (matplotlib & numpy)


## Overview

CyMiniAna is built for (and designed around) two main tasks in the ttbar charge asymmetry analyses:

  1. Event Loop (c++-based framework for speed)
     - Make smaller root files based on a selection criteria ("skimming" or "slimming")
     - Make histograms from existing root files (with or without applying an event selection)
     - Make efficiency curves from existing root files (with or without applying an event selection)
  2. Data/MC Plots (Python-based framework for improved aesthetics & ease-of-use)
     - Make stacked histogram plots that compare the (MC/data-driven) prediction with Data
     - Generic histograms (with support for ratio plots)
     - General efficiency curves (with support for drawing the underlying distribution)

The CyMiniAna Analysis Framework is structured as follows:

Directory  | About
---------  | ---------
python/    | plotting, neural network, and running scripts
src/       | `*.cxx` files
interface/ | `*.h` files
examples/  | Example scripts demonstrating how to use framework or general coding
config/    | text files and generic information
bin/       | steering macros
BuildFile.xml | File used by `scram` to compile software in CMSSW
Makefile   | compiles c++ code -- not needed in CMSSW framework!
README.md  | This file
setup.sh   | Setup script for cmssw environment



## Getting Started

Before running an analysis, check out the relevant tag for your studies or make a new branch to
perform your work.

**DO NOT USE THE MASTER BRANCH**

_The master branch is only used to make new tags!_

To get started, you first need a proper CMSSW release environment and other packages to work with CyMiniAna.  

```shell
## setup CMSSW (slc6_amd64_gcc530)
cmsrel CMSSW_8_0_24_patch1
cd CMSSW_8_0_24_patch1/src/
cmsenv
git cms-init
```

Next, it is necessary to add CyMiniAna and the external packages on which it depends.  
The lwtnn package, used for neural network inference in c++, can be included by 
following the steps [here](https://github.com/demarley/lwtnn/tree/CMSSW_8_0_X-compatible#cmssw-compatibility).

```
# Add lwtnn (using instructions above)
# Add necessary packages for dilepton reconstruction
git cms-addpkg CondFormats/JetMETObjects  # not explicitly used, but needed for Matrix Weighting
git cms-addpkg DataFormats/JetReco        # not explicitly used, available for using official "Jet" objects
git cms-addpkg JetMETCorrections/Modules  # Jet resolution smearing functions
mkdir cms-jet && cd cms-jet
git clone https://github.com/cms-jet/JRDatabase.git  # values for Jet resolution smearing
cd ../

## Add our code - the cms-ttbarAC packages!
git clone https://github.com/cms-ttbarAC/CyMiniAna.git Analysis/CyMiniAna
```

Once everything is checked out, compile it all!  
```
scram b -j8
```

To set the environment anytime you open a new session:  
```shell
source setup.csh   # ALWAYS DO THIS FIRST! (initializes everything)
```

Anytime you modify `*.cxx` or `*.h` code, be sure to recompile:
```shell
scram b -j8
```


## Running an Analysis

The following sections describes the steps needed to process an entire analysis,
from initial root files all the way to final histograms/plots for a note.

A single configuration file, e.g., `config/cmaConfig.txt`, sets the options
provided by the user to direct the program.  These options are listed here:


Option | About
------ | ------
NEvents       | number of events to process (a number < 0 defaults to 'ALL' events)
useJets       | use small-r jets in analysis (default `true`)
useLeptons    | use leptons in analysis (default `true`)
useLargeRJets | use large-R jets in analysis (default `false`)
useRCJets     | use re-clustered jets (fixed radius) in analysis (default `false`)
useNeutrinos  | use neutrinos in analysis (default `true`)
useTruth      | use & save truth information related to truth objects, e.g., jets (default `false`)
jet_btag_wkpt    | working point for jet b-tagging (from ATLAS, might be applicable in CMS)
buildNeutrinos   | reconstruct the neutrino solutions
makeNewFile      | generate a new root file(s) based on some selection from existing root file(s) (`true` or `false`)
makeHistograms   | generate histograms of some variables (with or without a selection) from existing root files (`true` or `false`)
makeEfficiencies | generate efficiency curves (with or without a selection) from existing root files (`true` or `false`)
calcWeightSystematics | Calculate the systematics associated with scale factors (e.g., b-tagging) (`true` or `false`)
weightSystematicsFile       | File containing the names of systematics stored as weights in `nominal` tree (from ATLAS)
weightVectorSystematicsFile | File containing the names of systematics stored as vectors of weights in `nominal` tree (from ATLAS)
input_selection  | level of selection applied to input files (Raw files may have different branches than files already processed by CyMiniAna) (`true` or `false`)
selection        | name for the level of selection to apply (used in `Root/eventSelection.cxx`)
output_path      | path for where to store the output
cutsfile         | text file that lists cuts' names and values (this is *not* used to apply cuts right now, just get the cut names for the cutflow histogram)
treenames        | the names of TTrees to process (e.g., "nominal", others for systematic uncertainties)
inputfile        | text file that lists all of the root files to process
sumWeightsFiles  | text file that lists all of the files needed to calculate sum of weights (monte carlo files only)
verboseLevel     | level to set the amount of `std::cout` statements to the console (options: DEBUG,INFO,WARNING,ERROR)
customFileEnding | Add a string to the end of filename for output files (uniquely identify outputs)
dnn_file         | The `*.json` file that contains DNN information for `lwtnn` tool (default `config/keras_ttbar_DNN.json` -- from ATLAS, not updated yet)
getDNN           | Calculate the DNN value (default `false`; if `false`, grab value from TTree)
getHME           | Calculate the HME value (default `false`; if `false`, grab value from TTree)
doRecoEventLoop  | Loop over reconstructed events (default `true`)
doTruthEventLoop | Loop over truth-level events (default `false`)
NJetSmear        | Number of smearings to perform for jet resolution (default 2)
NMassPoints      | Number of mass points to use in ttbar reconstruction (default `1`; more needed in jet mass measurement)
massMin          | Minimum top mass to use in ttbar reconstruction (default 172.5)


If these options aren't specified in the configuration file, 
default values will be chosen from `interface/configuration.h`.


### Event Loop

The c++-based framework within CyMiniAna builds the event for each entry in the ROOT file.
Physics objects (leptons, jets, etc.) are represented as structs within the framework (`interface/physicsObjects.h`).
The `Event` object is passed between classes (`histogrammer`, `efficiency`, etc.) to share the event information.
Other classes, such as those that build the dilepton ttbar system, are initialized from the `Event` class, provided with structs of physics objects, and return information back to the `Event` class.

Setup `config/cmaConfig.txt` (or your custom configuration file) and ensure any and 
all text files are also setup correctly, e.g., text file that points to the list of
root files you would like to process.  
To run the event selection code:

```shell
$ run config/cmaConfig.txt
```

#### Analysis Flow

1. The steering macro (usable "example" is `util/run.cxx`) first initializes and sets the configurations
    - Declare objects that are 'global' to all files being processed
2. File loop
    - Prepare output that is file-specific
        - Cutflow histograms, initialize output file, etc.
3. TTree loop
    - In ATLAS, the systematic uncertainties were stored as separate TTrees.
        - Looping over the TTrees lets you process all of the systematic uncertainties
    - Declare objects that are 'global' to all events (the event object, output ttree, histograms, and efficiencies)
4. Event Loop
    - Build the evnet object (jets, leptons, etc. for a given event)
    - Apply a selection (if needed)
    - Save information to TTree, histograms, or efficiencies. 


Different classes are used to achieve this information, 
and each one can be changed by the user.

Class     | About
--------- | ---------
configuration  | Class that contains all information for organization.  Multiple functions that return basic information as well
efficiency     | Class for generating efficiency curves (interface between TEfficiency and CyMiniAnaAC)
Event          | Class that contains all of the information from the event -> loads information from TTree and re-organizes information into structs & functions, calculates weights, etc.
eventSelection | Class to apply custom event selection (defined by user)
histogrammer   | Class for generating histograms (interface between TH1/TH2 and CyMiniAnaAC)
miniTree       | Class for writing & filling new ttree to output file
tools          | Collection of functions for doing simple tasks common to different aspects of CyMiniAnaAC
AMWT           | Class for applying the matrix weight technique in dilepton ttbar events
MassSolver     | Class for solving the ttbar dilepton equations

The steering macro is placed in the `bin/` directory.  This can be modified or extended by
the user -- preferably the user writes their own macro with similar functionality to `bin/run.cxx`.

If you add directories to the framework, ensure they will be compiled by checking `BuildFile.xml` and `bin/BuildFile.xml`.



### HepPlotter

To produce basic histograms, efficiency curves, or data/mc plots, the hepPlotter portion
of CyMiniAnaAC offers a simple interface to translate data (ROOT histograms/efficiencies or python lists/arrays) into plots.
In the `examples/hepPlotter` directory, there are multiple scripts at the user's disposal 
for making these plots.  *Python is preferred over using ROOT to make plots because 
of matplotlib's improved aesthetics and ease of use*.


To run the hepPlotter code (substitute `runHistogram.py` for your own script:

`$ python examples/hepPlotter/runHistogram.py --files examples/share/listOfFiles.txt --hists examples/share/listOfHists.txt -o ./`


#### General Histograms/Efficiencies

The `hepPlotter` class makes simple 1D and 2D plots with CMS formatting.  

Interfaces demonstrating how to make histograms or efficiency curves are shown in 
`runHistogram.py` and `runEfficiency.py`.  These two scripts differ in the way
the information is presented.  For efficiencies, the underlying distribution 
can also be drawn (e.g., jet pT spectrum to accompany jet trigger efficiency).
Internally to hepPlotter, TEfficiencies and TH*Ds are treated differently
only to access the data.  For further details, please see the functions 
hist2list(), hist2list2D(), TEfficiency2list(), and TEfficiency2list2D() in `python/hepPlotterTools.py`.

A typical setup to make a histogram plot requires the following steps:

```python
# Declare the plot
hist = HepPlotter(plotType,Ndims) # plotType = "histogram" or "efficiency"; Ndims = 1 or 2

# Modify some plot attributes (change default settings, see '__init__()' in hepPlotter class)
hist.XYZ = some_value

# Initialize
hist.initialize()

# Add data to the plot (this can be inside a for-loop if you want to add many histograms to same plot)
hist.Add(*args,**kwargs) # Can add TH*Ds, TEfficiencies (1D or 2D), python lists/arrays (not recommended!)

# Make the plot
hist.execute()

# Save the figure
hist.savefig()
```

#### Data/MC

To specifically make data/mc plots, use the class `hepPlotterDataMC`, which inherits from 
`hepPlotter` and removes all of the modularity parts to specifically draw data/mc plots.  
The example script for setting making data/mc plots is `examples/hepPlotter/runDataMC.py`.
The two frames (normal distributions of data and monte carlo and sub-frame that shows the ratio) are
drawn in the same plot.

*NOTE: HepPlotter can use numpy arrays (or python lists) or ROOT histograms. 
Histograms are __preferred__ because you can make them more quickly with CyMiniAna than 
looping over events with python & histogramming with numpy/matplotlib.*


# Contact

If there are questions, concerns, or bugs found in this package, please contact the author.  
If you would like contribute to the project, it is recommended to follow the 
[feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow).  
See `CONTRIBUTING.md` for more information.
