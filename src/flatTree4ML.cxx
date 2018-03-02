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
    m_ljet_deepAK8.resize(m_nDeepAK8);
    for (unsigned int i=0;i<m_nDeepAK8;i++){
        m_ljet_deepAK8[i] = 0;
        std::string idx   = std::to_string(i);
        m_ttree->Branch( ("ljet_deepAK8_"+idx).c_str(), &m_ljet_deepAK8.at(i), ("ljet_deepAK8_"+idx+"/F").c_str() );
    }

/*  ORIGINAL -- pre-DEEPAK8
    m_ttree->Branch( "ljet_SDmass", &m_ljet_SDmass, "ljet_SDmass/F" );
    m_ttree->Branch( "ljet_tau1",   &m_ljet_tau1,   "ljet_tau1/F" );
    m_ttree->Branch( "ljet_tau2",   &m_ljet_tau2,   "ljet_tau2/F" );
    m_ttree->Branch( "ljet_tau3",   &m_ljet_tau3,   "ljet_tau3/F" );
    m_ttree->Branch( "ljet_tau21",  &m_ljet_tau21,  "ljet_tau21/F" );
    m_ttree->Branch( "ljet_tau32",  &m_ljet_tau32,  "ljet_tau32/F" );
    m_ttree->Branch( "ljet_subjet0_bdisc", &m_ljet_subjet0_bdisc, "ljet_subjet0_bdisc/F" );
    m_ttree->Branch( "ljet_subjet0_pTrel", &m_ljet_subjet0_pTrel, "ljet_subjet0_pTrel/F" );
    //ljet_subjet0_charge
    m_ttree->Branch( "ljet_subjet1_bdisc", &m_ljet_subjet1_bdisc, "ljet_subjet1_bdisc/F" );
    m_ttree->Branch( "ljet_subjet1_pTrel", &m_ljet_subjet1_pTrel, "ljet_subjet1_pTrel/F" );
    //ljet_subjet1_charge
*/
    // AK4
    m_ttree->Branch( "jet_bdisc",  &m_jet_bdisc,  "jet_bdisc/F" );
    m_ttree->Branch( "jet_charge", &m_jet_charge, "jet_charge/F" );

    // AK8 + AK4 system
    m_ttree->Branch( "ljet_jet_m",      &m_ljet_jet_m,      "ljet_jet_m/F" );
    m_ttree->Branch( "ljet_jet_deltaR", &m_ljet_jet_deltaR, "ljet_jet_deltaR/F" );


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
    m_sumOfWeights = features.at("sumOfWeights");
    m_nominal_weight = features.at("nominal_weight");

    m_target = features.at("target");

    for (unsigned int i=0;i<m_nDeepAK8;i++){
        m_ljet_deepAK8.at(i) = features.at("ljet_deepAK8_"+std::to_string(i));
    }
/*
    m_ljet_SDmass = features.at("ljet_SDmass");
    m_ljet_tau1  = features.at("ljet_tau1");
    m_ljet_tau2  = features.at("ljet_tau2");
    m_ljet_tau3  = features.at("ljet_tau3");
    m_ljet_tau21 = features.at("ljet_tau21");
    m_ljet_tau32 = features.at("ljet_tau32");
    m_ljet_subjet0_bdisc = features.at("ljet_subjet0_bdisc");
    m_ljet_subjet0_pTrel = features.at("ljet_subjet0_pTrel");
    m_ljet_subjet1_bdisc = features.at("ljet_subjet1_bdisc");
    m_ljet_subjet1_pTrel = features.at("ljet_subjet1_pTrel");
*/

    m_jet_bdisc  = features.at("jet_bdisc");
    m_jet_charge = features.at("jet_charge");

    m_ljet_jet_m      = features.at("ljet_jet_m");
    m_ljet_jet_deltaR = features.at("ljet_jet_deltaR");


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

