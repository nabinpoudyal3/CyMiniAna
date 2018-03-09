/*
Created:        --
Last Updated:   28 August   2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Create and fill TTree.
*/
#include "Analysis/CyMiniAna/interface/miniTree.h"


miniTree::miniTree(configuration &cmaConfig) : 
  m_config(&cmaConfig){
    m_selections = m_config->selections();
  }

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
    // values based on the selection(s)
    m_passSelection.resize( m_selections.size() );
    unsigned int ss(0);
    for (const auto& sel : m_selections){
        m_passSelection.at(ss) = 0;
        m_ttree->Branch( sel.c_str(), &m_passSelection.at(ss), (sel+"/i").c_str() ); // unsigned int 0,1
        ss++;
    }

    if ( m_config->getDNN() )
        m_ttree->Branch( "DNN", &m_dnn, "DNN/F" );

    /*** disable branches here ***/
    // m_ttree->SetBranchStatus("", 0);

    return;
} // end initialize



void miniTree::saveEvent(Event& event, const std::vector<unsigned int>& evtsel_decisions) {
    /* Save the event to the ttree! */
    cma::DEBUG("MINITREE : Load the entry to be saved");
    m_oldTTree->GetEntry( event.entry() );  // make sure the original values are loaded for this event
                                            // otherwise only the branches accessed in Event are copied (!?)

    // set the new values
//    std::vector<Ljet> ljets = event.ljets();
//    for (const auto& ljet : ljets)
//        m_dnn = ljet.dnn;

    // set all decisions to false if they aren't passed here
    unsigned int n_sels = m_selections.size();
    bool generateDecisions = (evtsel_decisions.size()<1);

    for (unsigned int idx=0; idx<n_sels; idx++)
        m_passSelection.at(idx) = (generateDecisions) ? 0 : evtsel_decisions.at(idx); // set to 0 by default

    cma::DEBUG("MINITREE : Fill the tree");
    m_ttree->Fill();

    return;
}


void miniTree::finalize(){
    /* Finalize the class */
}

// THE END
