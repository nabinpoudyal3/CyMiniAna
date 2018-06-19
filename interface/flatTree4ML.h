#ifndef FLATTREE_H_
#define FLATTREE_H_

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <memory>
#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/physicsObjects.h"
#include "Analysis/CyMiniAna/interface/configuration.h"


class flatTree4ML {
  public:
    // Default - so root can load based on a name;
    flatTree4ML(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~flatTree4ML();

    // Run once at the start of the job;
    virtual void initialize(TFile& outputFile);

    // Run for every event (in every systematic) that needs saving;
    virtual void saveEvent(const std::map<std::string,double> features);

    // Clear stuff;
    virtual void finalize();


  protected:

    TTree * m_ttree;
    TTree * m_metadataTree;
    configuration * m_config;

    /**** Training branches ****/
    // weights for inputs
    float m_xsection;
    float m_kfactor;
    float m_sumOfWeights;
    float m_weight;
    float m_nominal_weight;

    // Deep learning features
    unsigned int m_target;

    unsigned int m_nDeepAK8;
    std::vector<float> m_ljet_deepAK8;

    float m_ljet_BEST_t;
    float m_ljet_BEST_w;
    float m_ljet_BEST_z;
    float m_ljet_BEST_h;
    float m_ljet_BEST_j;

    float m_ljet_SDmass;
    float m_ljet_tau1;
    float m_ljet_tau2;
    float m_ljet_tau3;
    float m_ljet_tau21;
    float m_ljet_tau32;
    unsigned int m_ljet_isHadTop;
    int m_ljet_contain;

    float m_ljet_charge;
    float m_ljet_subjet0_bdisc;
    float m_ljet_subjet0_charge;
    float m_ljet_subjet0_mass;
    float m_ljet_subjet0_mrel;
    float m_ljet_subjet0_ptrel;
    float m_ljet_subjet0_tau1;
    float m_ljet_subjet0_tau2;
    float m_ljet_subjet0_tau3;
    float m_ljet_subjet0_tau21;
    float m_ljet_subjet0_tau32;
    float m_ljet_subjet1_bdisc;
    float m_ljet_subjet1_charge;
    float m_ljet_subjet1_mass;
    float m_ljet_subjet1_mrel;
    float m_ljet_subjet1_ptrel;
    float m_ljet_subjet1_tau1;
    float m_ljet_subjet1_tau2;
    float m_ljet_subjet1_tau3;
    float m_ljet_subjet1_tau21;
    float m_ljet_subjet1_tau32;

    /**** Metadata ****/
    // which sample has which target value
    // many ROOT files will be merged together to do the training!
    std::string m_name;
    unsigned int m_target_value;
    unsigned int m_nEvents;
};

#endif

