/*
Created:        --
Last Updated:   12 April 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Create and fill TTree.
*/
#include "Analysis/CyMiniAna/interface/metadataTree.h"


metadataTree::metadataTree(configuration &cmaConfig) : 
  miniTree(cmaConfig),
  m_config(&cmaConfig){}

metadataTree::~metadataTree() {}



void metadataTree::initialize(TTree* t, TFile& outputFile, const std::string directory, const bool recalculateMetadata) {
    /* Setup the new tree

       @param t                    TTree from the original file
       @param directory            Directory the tree may be stored under
       @param recalculateMetadata  Value used in cloning the tree 
                                   (0=clone no events, just branch names; -1=clone all data)
    */
    outputFile.cd(directory.c_str());
    m_oldTTree = t;

    if (!recalculateMetadata){
        m_ttree = m_oldTTree->CloneTree(-1);
        cma::getListOfBranches(m_oldTTree,m_listOfBranches);
    }
    else{
        m_ttree = new TTree("metadata","metadata");
        createBranches();
    }

    return;
}


void metadataTree::createBranches(){
    /* Setup branches */
    m_ttree->Branch("primaryDataset", &m_sampleName);                      // string
    m_ttree->Branch("xsection",       &m_xsection,     "xsection/F");      // float
    m_ttree->Branch("kfactor",        &m_kfactor,      "kfactor/F");       // float
    m_ttree->Branch("sumOfWeights",   &m_sumOfWeights, "sumOfWeights/F");  // float
    m_ttree->Branch("NEvents",        &m_NEvents,      "NEvents/i");       // uint

    return;
}


void metadataTree::saveMetaData(const Sample& smp){
    /* Save the metadata to the ttree! */
    cma::DEBUG("METADATA : Load the entry to be saved");

    m_sampleName   = smp.primaryDataset;
    m_sumOfWeights = smp.sumOfWeights;
    m_xsection = smp.XSection;
    m_kfactor  = smp.KFactor;
    m_NEvents  = smp.NEvents;

    cma::DEBUG("METADATA : Fill the tree");
    m_ttree->Fill();

    return;
}

// THE END
