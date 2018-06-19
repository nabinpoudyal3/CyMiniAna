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
    m_ttree->Branch( "ljet_charge",  &m_ljet_charge,  "ljet_charge/F" );
    m_ttree->Branch( "ljet_subjet0_bdisc",  &m_ljet_subjet0_bdisc,  "ljet_subjet0_bdisc/F" );
    m_ttree->Branch( "ljet_subjet0_charge", &m_ljet_subjet0_charge, "ljet_subjet0_charge/F" );
    m_ttree->Branch( "ljet_subjet0_mass",   &m_ljet_subjet0_mass,   "ljet_subjet0_mass/F" );
    m_ttree->Branch( "ljet_subjet0_mrel",   &m_ljet_subjet0_mrel,   "ljet_subjet0_mrel/F" );
    m_ttree->Branch( "ljet_subjet0_ptrel",  &m_ljet_subjet0_ptrel,  "ljet_subjet0_ptrel/F" );
    m_ttree->Branch( "ljet_subjet0_tau1",   &m_ljet_subjet0_tau1,   "ljet_subjet0_tau1/F" );
    m_ttree->Branch( "ljet_subjet0_tau2",   &m_ljet_subjet0_tau2,   "ljet_subjet0_tau2/F" );
    m_ttree->Branch( "ljet_subjet0_tau3",   &m_ljet_subjet0_tau3,   "ljet_subjet0_tau3/F" );
    m_ttree->Branch( "ljet_subjet0_tau21",  &m_ljet_subjet0_tau21,  "ljet_subjet0_tau21/F" );
    m_ttree->Branch( "ljet_subjet0_tau32",  &m_ljet_subjet0_tau32,  "ljet_subjet0_tau32/F" );
    m_ttree->Branch( "ljet_subjet1_bdisc",  &m_ljet_subjet1_bdisc,  "ljet_subjet1_bdisc/F" );
    m_ttree->Branch( "ljet_subjet1_charge", &m_ljet_subjet1_charge, "ljet_subjet1_charge/F" );
    m_ttree->Branch( "ljet_subjet1_mass",   &m_ljet_subjet1_mass,   "ljet_subjet1_mass/F" );
    m_ttree->Branch( "ljet_subjet1_mrel",   &m_ljet_subjet1_mrel,   "ljet_subjet1_mrel/F" );
    m_ttree->Branch( "ljet_subjet1_ptrel",  &m_ljet_subjet1_ptrel,  "ljet_subjet1_ptrel/F" );
    m_ttree->Branch( "ljet_subjet1_tau1",   &m_ljet_subjet1_tau1,   "ljet_subjet1_tau1/F" );
    m_ttree->Branch( "ljet_subjet1_tau2",   &m_ljet_subjet1_tau2,   "ljet_subjet1_tau2/F" );
    m_ttree->Branch( "ljet_subjet1_tau3",   &m_ljet_subjet1_tau3,   "ljet_subjet1_tau3/F" );
    m_ttree->Branch( "ljet_subjet1_tau21",  &m_ljet_subjet1_tau21,  "ljet_subjet1_tau21/F" );
    m_ttree->Branch( "ljet_subjet1_tau32",  &m_ljet_subjet1_tau32,  "ljet_subjet1_tau32/F" );

    // AK8
    m_ttree->Branch( "ljet_BEST_t", &m_ljet_BEST_t, "ljet_BEST_t/F" );
    m_ttree->Branch( "ljet_BEST_w", &m_ljet_BEST_w, "ljet_BEST_w/F" );
    m_ttree->Branch( "ljet_BEST_z", &m_ljet_BEST_z, "ljet_BEST_z/F" );
    m_ttree->Branch( "ljet_BEST_h", &m_ljet_BEST_h, "ljet_BEST_h/F" );
    m_ttree->Branch( "ljet_BEST_j", &m_ljet_BEST_j, "ljet_BEST_j/F" );
    m_ttree->Branch( "ljet_SDmass", &m_ljet_SDmass, "ljet_SDmass/F" );
    m_ttree->Branch( "ljet_tau1",   &m_ljet_tau1,   "ljet_tau1/F" );
    m_ttree->Branch( "ljet_tau2",   &m_ljet_tau2,   "ljet_tau2/F" );
    m_ttree->Branch( "ljet_tau3",   &m_ljet_tau3,   "ljet_tau3/F" );
    m_ttree->Branch( "ljet_tau21",  &m_ljet_tau21,  "ljet_tau21/F" );
    m_ttree->Branch( "ljet_tau32",  &m_ljet_tau32,  "ljet_tau32/F" );
    m_ttree->Branch( "ljet_isHadTop",&m_ljet_isHadTop, "ljet_isHadTop/i" );
    m_ttree->Branch( "ljet_contain", &m_ljet_contain,  "ljet_contain/I" );

    /**** Metadata ****/
    // which sample has which target value
    // many ROOT files will be merged together to do the training
    m_metadataTree->Branch( "name",    &m_name );
    m_metadataTree->Branch( "target",  &m_target_value,  "target/I" );    // useful if processing targets stored in different samples
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
    m_ljet_subjet0_mass   = features.at("ljet_subjet0_mass");
    m_ljet_subjet0_mrel   = features.at("ljet_subjet0_mrel");
    m_ljet_subjet0_ptrel  = features.at("ljet_subjet0_ptrel");
    m_ljet_subjet0_tau1   = features.at("ljet_subjet0_tau1");
    m_ljet_subjet0_tau2   = features.at("ljet_subjet0_tau2");
    m_ljet_subjet0_tau3   = features.at("ljet_subjet0_tau3");
    m_ljet_subjet0_tau21  = features.at("ljet_subjet0_tau21");
    m_ljet_subjet0_tau32  = features.at("ljet_subjet0_tau32");

    m_ljet_subjet1_bdisc  = features.at("ljet_subjet1_bdisc");
    m_ljet_subjet1_charge = features.at("ljet_subjet1_charge");
    m_ljet_subjet1_mass   = features.at("ljet_subjet1_mass");
    m_ljet_subjet1_mrel   = features.at("ljet_subjet1_mrel");
    m_ljet_subjet1_ptrel  = features.at("ljet_subjet1_ptrel");
    m_ljet_subjet1_tau1   = features.at("ljet_subjet1_tau1");
    m_ljet_subjet1_tau2   = features.at("ljet_subjet1_tau2");
    m_ljet_subjet1_tau3   = features.at("ljet_subjet1_tau3");
    m_ljet_subjet1_tau21  = features.at("ljet_subjet1_tau21");
    m_ljet_subjet1_tau32  = features.at("ljet_subjet1_tau32");

    m_ljet_BEST_t = features.at("ljet_BEST_t");
    m_ljet_BEST_w = features.at("ljet_BEST_w");
    m_ljet_BEST_z = features.at("ljet_BEST_z");
    m_ljet_BEST_h = features.at("ljet_BEST_h");
    m_ljet_BEST_j = features.at("ljet_BEST_j");
    m_ljet_SDmass = features.at("ljet_SDmass");
    m_ljet_tau1   = features.at("ljet_tau1");
    m_ljet_tau2   = features.at("ljet_tau2");
    m_ljet_tau3   = features.at("ljet_tau3");
    m_ljet_tau21  = features.at("ljet_tau21");
    m_ljet_tau32  = features.at("ljet_tau32");
    m_ljet_isHadTop = static_cast<unsigned int>(features.at("ljet_isHadTop"));
    m_ljet_contain  = static_cast<int>(features.at("ljet_contain"));

    /**** Fill the tree ****/
    cma::DEBUG("FLATTREE4ML : had top "+std::to_string(features.at("ljet_isHadTop")));
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
