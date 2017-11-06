/*
Created:        --
Last Updated:   28 August   2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Create and fill TTree.
*/
#include "cms-ttbarAC/CyMiniAna/interface/miniTree.h"


miniTree::miniTree(configuration &cmaConfig) : 
  m_config(&cmaConfig){}

miniTree::~miniTree() {}



void miniTree::initialize(TTree* t, TFile& outputFile) {
    /* 
       Setup the new tree 
       't' represents the TTree from the original file
    */
    outputFile.cd();
    m_oldTTree = t;

    m_ttree = m_oldTTree->CloneTree(0);  // clone the tree (branches) but copy no data

    /*** setup new branches here ***/
    m_ttree->Branch( "DNN", &m_dnn );

    /*** disable branches here ***/
    // m_ttree->SetBranchStatus("", 0);

    return;
} // end initialize



void miniTree::saveEvent(Event& event) {
    /* Save the event to the ttree! */
    cma::DEBUG("MINITREE : Load the entry to be saved");
    m_oldTTree->GetEntry( event.entry() );  // make sure the original values are loaded for this event
                                            // otherwise only the branches accessed in Event are copied (!?)

    // set the new values
    m_dnn = event.DNN();

    cma::DEBUG("MINITREE : Fill the tree");
    m_ttree->Fill();

    return;
}


void miniTree::finalize(){
    /* Finalize the class */
}

// THE END
