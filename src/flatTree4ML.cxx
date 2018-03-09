/*
Created:        --
Last Updated:   20 February 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Create and fill TTree for ML
*/
#include "Analysis/CyMiniAna/interface/flatTree4ML.h"


flatTree4ML::flatTree4ML(configuration &cmaConfig) : 
  m_config(&cmaConfig),
  m_nDeepAK8(16){}

flatTree4ML::~flatTree4ML() {}


void flatTree4ML::initialize(TFile& outputFile) {
    /*
       Setup the new tree 
       Contains features for the NN
       --  No vector<T> stored in tree: completely flat!
    */
    outputFile.cd();                                     // move to output file
    m_ttree        = new TTree("features", "features");  // Tree contains features for the NN
    m_metadataTree = new TTree("metadata","metadata");   // Tree contains metadata

    /**** Setup new branches here ****/
    // Weights
    m_ttree->Branch( "xsection", &m_xsection, "xsection/F" );
    m_ttree->Branch( "kfactor",  &m_kfactor,  "kfactor/F" );
    m_ttree->Branch( "weight",   &m_weight,   "weight/F" );
    m_ttree->Branch( "sumOfWeights",   &m_sumOfWeights,   "sumOfWeights/F" );
    m_ttree->Branch( "nominal_weight", &m_nominal_weight, "nominal_weight/F" );

    // Features
    m_ttree->Branch( "target", &m_target, "target/I" );  // target value (.e.g, 0 or 1)

    // AK8
    m_ttree->Branch( "ljet_charge",  &m_ljet_charge,  "ljet_charge/F" );
    m_ttree->Branch( "ljet_subjet0_bdisc", &m_ljet_subjet0_bdisc,   "ljet_subjet0_bdisc/F" );
    m_ttree->Branch( "ljet_subjet0_charge", &m_ljet_subjet0_charge, "ljet_subjet0_charge/F");
    m_ttree->Branch( "ljet_subjet1_bdisc", &m_ljet_subjet1_bdisc,   "ljet_subjet1_bdisc/F" );
    m_ttree->Branch( "ljet_subjet1_charge", &m_ljet_subjet1_charge, "ljet_subjet1_charge/F");


    /**** Metadata ****/
    // which sample has which target value
    // many ROOT files will be merged together to do the training
    m_metadataTree->Branch( "name",    &m_name );
    m_metadataTree->Branch( "target",  &m_target_value,  "target/I" );
    m_metadataTree->Branch( "nEvents", &m_nEvents,       "nEvents/I" );

    return;
} // end initialize



void flatTree4ML::saveEvent(const std::map<std::string,double> features) {
    /* Save the ML features to the ttree! */
    cma::DEBUG("FLATTREE4ML : Save event ");

    m_weight   = features.at("weight");
    m_kfactor  = features.at("kfactor");
    m_xsection = features.at("xsection");
    m_sumOfWeights   = features.at("sumOfWeights");
    m_nominal_weight = features.at("nominal_weight");

    m_target = features.at("target");

    m_ljet_charge = features.at("ljet_charge");
    m_ljet_subjet0_bdisc  = features.at("ljet_subjet0_bdisc");
    m_ljet_subjet0_charge = features.at("ljet_subjet0_charge");
    m_ljet_subjet1_bdisc  = features.at("ljet_subjet1_bdisc");
    m_ljet_subjet1_charge = features.at("ljet_subjet1_charge");

    /**** Fill the tree ****/
    cma::DEBUG("FLATTREE4ML : Fill the tree");
    m_ttree->Fill();

    return;
}


void flatTree4ML::finalize(){
    /* Finalize the class -- fill in the metadata (only need to do this once!) */
    m_name    = m_config->primaryDataset();
    m_nEvents = m_config->NTotalEvents();
    m_target_value = (m_config->isQCD()) ? 0 : -1;    // multiple classes for signal, choose '-1'

    cma::DEBUG("FLATTREE4ML : Fill the metadata tree");
    m_metadataTree->Fill();
}

// THE END