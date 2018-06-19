/*
Created:        --
Last Updated:    2 March 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Tool for performing deep learning tasks
- Inference: LWTNN
- Training: save features to flat ntuple; train in python environment

-- Setup as of 2 March: Top/Antitop tagging (CHEETAH)
*/
#include "Analysis/CyMiniAna/interface/deepLearning.h"


deepLearning::deepLearning( configuration& cmaConfig ) :
  m_config(&cmaConfig),
  m_lwnn(nullptr){
    m_features.clear();

    // Setup lwtnn
    m_dnnKey = m_config->dnnKey();
    if (m_config->DNNinference()){
      std::ifstream input_cfg = cma::open_file( m_config->dnnFile() );
      lwt::JSONConfig cfg     = lwt::parse_json( input_cfg );
      m_lwnn   = new lwt::LightweightNeuralNetwork(cfg.inputs, cfg.layers, cfg.outputs);
    }
  }

deepLearning::~deepLearning() {
    delete m_lwnn;
}


void deepLearning::training(std::vector<Ljet>& ljets){
    /* Prepare inputs for training -- saving to ROOT file in flatTree4ML */
    for (auto& ljet : ljets)
        training(ljet);

    return;
}

void deepLearning::training(Ljet& ljet){
    /* Prepare inputs for training -- saving to ROOT file in flatTree4ML */
    loadFeatures(ljet);
    ljet.features = m_features;
    return;
}

void deepLearning::inference(std::vector<Ljet>& ljets){
    /* Obtain results from LWTNN */
    for (auto& ljet : ljets)
        inference(ljet);
    return;
}

void deepLearning::inference(Ljet& ljet){
    /* Calculate DNN prediction */
    loadFeatures(ljet);
    m_predictions = m_lwnn->compute(m_features);
    ljet.dnn = m_predictions;
    m_DNN    = m_predictions.at(m_dnnKey);        // set default value

    return;
}


void deepLearning::loadFeatures(const Ljet& ljet){
    /* Calculate DNN features */
    m_features.clear();

    // feature calculations
    m_features["target"] = ljet.target;

    m_features["ljet_subjet0_bdisc"]  = ljet.subjet0_bdisc;
    m_features["ljet_subjet0_charge"] = ljet.subjet0_charge;
    m_features["ljet_subjet0_mass"]   = ljet.subjet0_mass;
    m_features["ljet_subjet0_mrel"]   = ljet.subjet0_mass/ljet.p4.M();
    m_features["ljet_subjet0_ptrel"]  = ljet.subjet0_pt/ljet.p4.Pt();
    m_features["ljet_subjet0_tau1"]   = ljet.subjet0_tau1;
    m_features["ljet_subjet0_tau2"]   = ljet.subjet0_tau2;
    m_features["ljet_subjet0_tau3"]   = ljet.subjet0_tau3;
    m_features["ljet_subjet0_tau21"]  = ljet.subjet0_tau2 / ljet.subjet0_tau1;
    m_features["ljet_subjet0_tau32"]  = ljet.subjet0_tau3 / ljet.subjet0_tau2;
    m_features["ljet_subjet1_bdisc"]  = ljet.subjet1_bdisc;
    m_features["ljet_subjet1_charge"] = ljet.subjet1_charge;
    m_features["ljet_subjet1_mass"]   = ljet.subjet1_mass;
    m_features["ljet_subjet1_mrel"]   = ljet.subjet1_mass/ljet.p4.M();
    m_features["ljet_subjet1_ptrel"]  = ljet.subjet1_pt/ljet.p4.Pt();
    m_features["ljet_subjet1_tau1"]   = ljet.subjet1_tau1;
    m_features["ljet_subjet1_tau2"]   = ljet.subjet1_tau2;
    m_features["ljet_subjet1_tau3"]   = ljet.subjet1_tau3;
    m_features["ljet_subjet1_tau21"]  = ljet.subjet1_tau2 / ljet.subjet1_tau1;
    m_features["ljet_subjet1_tau32"]  = ljet.subjet1_tau3 / ljet.subjet1_tau2;

    m_features["ljet_charge"] = ljet.charge;

    m_features["weight"] = 1.;  // 1/ljet.p4.Pt() or something

    cma::DEBUG("EVENT : Set DNN input values ");

    return;
}

double deepLearning::prediction(const std::string& key) const{
    /* Just return the prediction (after execute!) */
    return m_predictions.at(key);
}

// THE END //
