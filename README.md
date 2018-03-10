10 March 2018  
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
bin/       | steering macros
config/    | text files and generic information
examples/  | Example scripts demonstrating how to use framework or general coding
interface/ | `*.h` files
python/    | plotting, neural network, and running scripts
src/       | `*.cxx` files
BuildFile.xml | File used by `scram` to compile software in CMSSW
Makefile   | compiles c++ code -- not needed in CMSSW framework!
README.md  | This file
setup.sh   | Setup script for cmssw environment



## Getting Started

Before running an analysis, make a new branch to perform your work.  
It is recommend that you make a new branch from a specific tag.

**DO NOT USE THE MASTER BRANCH FOR DEVELOPMENT OR RUNNING AN ANALYSIS**  
_The master branch is only used to make new tags!_

To get started, you first need a proper CMSSW release environment and other packages to work with CyMiniAna.  

```shell
## setup CMSSW (slc6_amd64_gcc530)
cmsrel CMSSW_9_0_1
cd CMSSW_9_0_1/src/
cmsenv
git cms-init
```

Next, it is necessary to add CyMiniAna and the external packages on which it depends.  
The lwtnn package, used for neural network inference in c++, can be included by 
following the steps [here](https://github.com/demarley/lwtnn/tree/CMSSW_8_0_X-compatible#cmssw-compatibility).

```
## Add lwtnn (using instructions above)

## Add our code - the cms-ttbarAC packages!
mkdir Analysis
git clone https://github.com/cms-ttbarAC/CyMiniAna.git Analysis/CyMiniAna
```

Once everything is checked out, compile!  
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

# More Information

Please consult the wiki pages for more information.  
If you cannot find what you are looking for, please see the "Contact" section below 
or create an issue.

# Contact

If there are questions, concerns, or bugs found in this package, please contact the author.  
If you would like contribute to the project, it is recommended to follow the 
[feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow).  
See `CONTRIBUTING.md` for more information.
