## Some software is readily available on cvmfs / cmssw

set sourced=($_)
setenv CYMINIANADIR `dirname $sourced[2]`
setenv CERN_USER ${USER} # put your CERN username if different from enviornment name

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
setenv PYTHONPATH ${PYTHONPATH}:${PWD}/python
setenv PYTHONPATH ${PYTHONPATH}:${PWD}/examples/hepPlotter

echo " > Setup LaTeX "
setenv PATH /afs/cern.ch/sw/XML/texlive/latest/bin/x86_64-linux:$PATH
setenv TEXINPUTS ${HOME}/texmf-CERN:


# Set grid proxy if not provided
echo ""
echo "   Grid info: "
voms-proxy-info -exists > /dev/null
setenv global_proxy_ok $?
if ( ${global_proxy_ok} != 0 ) then
    echo No valid grid proxy found. Creating new one...
    voms-proxy-init -voms cms 
    if ( $? != 0 ) then
        echo Failed to create grid proxy.
    endif
else
    echo Grid proxy already set.
endif

