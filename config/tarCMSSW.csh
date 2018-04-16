# Created: 9 April 2018
#
# Dan Marley
# daniel.edison.marley@cernSPAMNOT.ch
# Texas A&M University
# ---
# Tar the current CMSSW directory for running condor jobs


# set path for where to put the tarball
set eos_path = "/store/user/lpctop/ttbarAC/flatNtuples/CWoLa/"


echo " Make tarball ..."

tar --exclude-caches-all --exclude-vcs -zcf $CMSSW_VERSION.tgz -C $CMSSW_BASE/.. $CMSSW_VERSION \
    --exclude=$CMSSW_VERSION/src/Analysis/CyMiniAna/data \
    --exclude=$CMSSW_VERSION/src/BESTAnalysis \
    --exclude=$CMSSW_VERSION/src/ttbarAC_skim \
    --exclude=analysis_note \
    --exclude=cms-ttbarAC-AN \
    --exclude=$CMSSW_VERSION/src/PhysicsTools \
    --exclude=$CMSSW_VERSION/src/RecoEgamma \
    --exclude=$CMSSW_VERSION/src/BESTAnalysis \
    --exclude=$CMSSW_VERSION/include/slc6_amd64_gcc630/PhysicsTools \
    --exclude=$CMSSW_VERSION/tmp/slc6_amd64_gcc630/src/PhysicsTools \
    --exclude=$CMSSW_VERSION/include/slc6_amd64_gcc630/RecoEgamma \
    --exclude=$CMSSW_VERSION/tmp/slc6_amd64_gcc630/src/RecoEgamma \
    --exclude=$CMSSW_VERSION/include/slc6_amd64_gcc630/BESTAnalysis \
    --exclude=$CMSSW_VERSION/tmp/slc6_amd64_gcc630/src/BESTAnalysis \
    --verbose


echo " Copy tarball to new location on EOS "$eos_path
xrdcp $CMSSW_VERSION.tgz root://cmseos.fnal.gov/"$eos_path"
