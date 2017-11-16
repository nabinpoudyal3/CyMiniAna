## Some software is readily available on cvmfs / cmssw

export CERN_USER=${USER} # put your CERN username if different from enviornment name
export CYMINIANADIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export JRDATABASE=$CYMINIANADIR/../../cms-jet/JRDatabase/

echo ""
echo " * ------------------------------------------------- * "
echo " *                     CyMiniAna                     * "
echo " * ------------------------------------------------- * "
echo " * Framework to perform event selection, write-out   * "
echo " * a few histograms or efficiencies, and make plots. * "
echo " * ------------------------------------------------- * " 
 
echo ""
echo " > Setup CMS "
cmsenv

echo " > Setup Python "
export PYTHONPATH=$PYTHONPATH:$PWD/python
export PYTHONPATH=$PYTHONPATH:$PWD/examples/hepPlotter

echo " > Setup LaTeX "
export PATH=/afs/cern.ch/sw/XML/texlive/latest/bin/x86_64-linux:$PATH
export TEXINPUTS=${HOME}/texmf-CERN:

echo " > Setup LHAPDF "
export LHAPDF_DATA_PATH=/cvmfs/cms.cern.ch/lhapdf/pdfsets/6.1.h/
#cteq6l1/

# Set grid proxy if not provided
echo ""
echo "   Grid info: "
voms-proxy-info -exists 1> /dev/null 2> /dev/null
global_proxy_ok=$?
if [ $global_proxy_ok -ne 0 ]; then
    echo No valid grid proxy found. Creating new one...
    voms-proxy-init -voms cms --valid 96:00 #-out=$jobdir
if [ $? -ne 0 ]; then
    echo Failed to create grid proxy.
    fi
fi
