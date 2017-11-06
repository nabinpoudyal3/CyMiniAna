#ifndef EVENT_H_
#define EVENT_H_

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TParameter.h"
#include "TEnv.h"
#include "TF1.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cms-ttbarAC/CyMiniAna/interface/physicsObjects.h"
#include "cms-ttbarAC/CyMiniAna/interface/configuration.h"
#include "cms-ttbarAC/CyMiniAna/interface/dileptonTtbarReco.h"
#include "cms-ttbarAC/lwtnn/interface/LightweightNeuralNetwork.hh"
#include "cms-ttbarAC/lwtnn/interface/parse_json.hh"


// Event Class
class Event {
  public:
    // Constructor
    Event( TTreeReader &myReader, configuration &cmaConfig);
    Event( const Event &obj);

    // Destructor
    virtual ~Event();

    // create hash tables in truth/reco TTree to match truth <-> reco events
    // uses configuration option matchTruthToReco to match truth to reco (reco loop)
    // OR match reco to truth (truth loop, for acceptance studies)
    void matchTruthWithReco();
    // check during looping over truth events, if reco event match is found
    bool isValidRecoEntry();

    // Execute the event (load information and setup objects)
    virtual void execute(Long64_t entry);
    virtual void updateEntry(Long64_t entry);

    // Setup physics information
    void initialize_leptons();
    void initialize_neutrinos();
    void initialize_jets();
    void initialize_ljets();
    void initialize_eventWeights();
    void initialize_weights();
    void initialize_kinematics();
    void initialize_truth();

    virtual double getSystEventWeight(const std::string &syst, const int weightIndex=-1);

    // Clear stuff;
    virtual void finalize();
    virtual void clear();

    // Get physics information
    std::vector<Lepton> leptons();
    std::vector<Neutrino> neutrinos();
    std::vector<Ljet> ljets();
    std::vector<Jet>  jets();

    // Get truth physics information 
    std::vector<Lepton> truth_leptons();
    std::vector<Neutrino> truth_neutrinos();
    std::vector<Ljet> truth_ljets();
    std::vector<Jet>  truth_jets();

    virtual float HT();
    virtual float ST();
    virtual std::vector<int> btag_jets(const std::string &wkpt);
    virtual std::vector<int> btag_jets(); // using configured b-tag WP

//    virtual unsigned long long eventNumber();
//    virtual unsigned int runNumber();
    virtual float eventNumber();
    virtual float runNumber();
    virtual unsigned int mcChannelNumber();
    virtual float mu();
    virtual int lumiblock();
    virtual float xsection();
    virtual float kfactor();
    virtual float sumOfWeights();
    virtual float met( const std::string& met_name );

    void getDNNInputs();      // return the DNN inputs to the user
    void getDNN();            // get the DNN output
    void getHME();
    double DNN();
    double HME();
    void buildTtbar();
    void getDilepton();
    std::map<std::string,Top> ttbar();
    void truth();
    long long entry();

    // Get weights
    virtual float nominal_weight();
    float weight_mc();
    float truth_weight_mc();
    float weight_jvt();
    float weight_pileup();
    float weight_lept_eff();
    float weight_btag();
    float weight_btag(const std::string &wkpt);

    // Get weight systematics
    virtual std::map<std::string,float > weightSystematicsFloats();
    virtual std::map<std::string,std::vector<float> > weightSystematicsVectorFloats();

    virtual std::vector<std::string> listOfWeightSystematics();
    virtual std::string treeName();

  protected:

    // general information
    configuration *m_config;
    TTreeReader &m_ttree;
    TTreeReader m_truth_tree;
    std::string m_treeName;
    bool m_grid;
    bool m_isMC;
    long long m_entry;
    long long m_truth_entry;
    bool m_getDNN;
    bool m_getHME;
    bool m_buildNeutrinos;

    bool m_ee;
    bool m_mumu;
    bool m_emu;

    // event weight information
    double m_nominal_weight;
    double m_xsection;
    double m_kfactor;
    double m_sumOfWeights;
    double m_LUMI;
    std::map<int, float> m_mapXSection; // map DSID to XSection
    std::map<int, float> m_mapKFactor;  // map DSID to KFactor
    std::map<int, float> m_mapAMI;      // map DSID to sum of weights

    // physics object information
    std::vector<Lepton> m_leptons;
    std::vector<Neutrino> m_neutrinos;
    std::vector<Ljet> m_ljets;
    std::vector<Jet>  m_jets;

    // truth physics object information
    std::vector<Lepton> m_truth_leptons;
    std::vector<Neutrino> m_truth_neutrinos;
    std::vector<Ljet> m_truth_ljets;
    std::vector<Jet>  m_truth_jets;


    // b-tagged calo jets with various WP
    std::map<std::string, std::vector<int> > m_btag_jets;
    std::vector<int> m_btag_jets_default;
    float m_cMVAv2L;
    float m_cMVAv2M;
    float m_cMVAv2T;

    float m_HT;
    float m_ST;
    float m_metmet;
    float m_metphi;

    lwt::LightweightNeuralNetwork* m_lwnn;
    std::map<std::string, double> m_dnnInputs;   // values for inputs to the DNN
    std::string m_dnnKey;
    double m_DNN;   // DNN score
    double m_HME;   // HME score

    dileptonTtbarReco* m_dileptonTtbar;
    std::map<std::string,Top> m_ttbar;
    DileptonReco m_dilepton;

    // nominal b-tagging weight maps
    std::map<std::string, float> m_weight_trackjet_btag;
    std::map<std::string, float> m_weight_btag;
    float m_weight_btag_default;
    // Maps to keep track of weight systematics
    std::map<std::string,TTreeReaderValue<float> * > m_weightSystematicsFloats;
    std::map<std::string,TTreeReaderValue<std::vector<float>> * > m_weightSystematicsVectorFloats;
    std::vector<std::string> m_listOfWeightSystematics;

    // TTree variables [all possible ones]

    // *************
    // the following are from root files accessed 
    //    on 28 August 2017 
    //    at /fdata/hepx/store/user/tahuang/HHNTuples
    TTreeReaderValue<float> * m_llmetjj_DPhi_ll_jj;
    TTreeReaderValue<float> * m_llmetjj_minDR_l_j;
    TTreeReaderValue<float> * m_llmetjj_MTformula;
    TTreeReaderValue<float> * m_llmetjj_MT2;
    TTreeReaderValue<float> * m_llmetjj_M;
    TTreeReaderValue<float> * m_lljj_M;

    TTreeReaderValue<float> * m_lep1_pt;
    TTreeReaderValue<float> * m_lep2_pt;
    TTreeReaderValue<float> * m_lep1_eta;
    TTreeReaderValue<float> * m_lep1_phi;
    TTreeReaderValue<float> * m_lep1_Iso;
    TTreeReaderValue<float> * m_lep2_eta;
    TTreeReaderValue<float> * m_lep2_phi;
    TTreeReaderValue<float> * m_lep2_Iso;
    TTreeReaderValue<float> * m_ll_M;
    TTreeReaderValue<float> * m_ll_pt;
    TTreeReaderValue<float> * m_ll_DR_l_l;
    TTreeReaderValue<float> * m_ll_DPhi_l_l;
    TTreeReaderValue<float> * m_ll_DEta_l_l;
    TTreeReaderValue<float> * m_llidiso;
    TTreeReaderValue<float> * m_mumuidiso;
    TTreeReaderValue<float> * m_elelidiso;

    TTreeReaderValue<float> * m_isElEl;
    TTreeReaderValue<float> * m_isMuMu;
    TTreeReaderValue<float> * m_isElMu;
    TTreeReaderValue<float> * m_isMuEl;

    TTreeReaderValue<float> * m_jet1_pt;
    TTreeReaderValue<float> * m_jet2_pt;
    TTreeReaderValue<float> * m_jet1_eta;
    TTreeReaderValue<float> * m_jet1_phi;
    TTreeReaderValue<float> * m_jet2_eta;
    TTreeReaderValue<float> * m_jet2_phi;
    TTreeReaderValue<float> * m_jet1_cMVAv2;
    TTreeReaderValue<float> * m_jet2_cMVAv2;
    TTreeReaderValue<float> * m_nJetsL;
    TTreeReaderValue<float> * m_jjbtag_heavy;
    TTreeReaderValue<float> * m_jjbtag_light;
    TTreeReaderValue<float> * m_jj_DR_j_j;
    TTreeReaderValue<float> * m_jj_pt;
    TTreeReaderValue<float> * m_jj_M;
    TTreeReaderValue<float> * m_met_pt;
    TTreeReaderValue<float> * m_ht;
    TTreeReaderValue<float> * m_hme_h2mass_reco;
    TTreeReaderValue<float> * m_dnn_score;

    TTreeReaderValue<float> * m_cosThetaStar;
    TTreeReaderValue<bool> * m_isSF;

    TTreeReaderValue<float> * m_trigeff;
    TTreeReaderValue<float> * m_pu;
    TTreeReaderValue<float> * m_sample_weight;
    TTreeReaderValue<float> * m_event_weight;
    TTreeReaderValue<float> * m_event_pu_weight;
    TTreeReaderValue<float> * m_event_number;
    TTreeReaderValue<float> * m_event_run;
    TTreeReaderValue<float> * m_total_weight;
    // *************


    TTreeReaderValue<float> * m_weight_mc;
    TTreeReaderValue<float> * m_weight_jvt;
    TTreeReaderValue<float> * m_weight_pileup;
    TTreeReaderValue<float> * m_weight_lept_eff;
    TTreeReaderValue<float> * m_weight_btag_85;
    TTreeReaderValue<float> * m_weight_btag_77;
    TTreeReaderValue<float> * m_weight_btag_70;
    TTreeReaderValue<float> * m_weight_btag_60;
    TTreeReaderValue<float> * m_weight_trackjet_btag_70;
    TTreeReaderValue<float> * m_weight_trackjet_btag_77;
    TTreeReaderValue<float> * m_weight_jvt_UP;
    TTreeReaderValue<float> * m_weight_jvt_DOWN;    
    TTreeReaderValue<float> * m_weight_pileup_UP;
    TTreeReaderValue<float> * m_weight_pileup_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_Trigger_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_Trigger_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_Reco_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_Reco_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_ID_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_ID_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_Isol_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_EL_SF_Isol_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Trigger_STAT_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Trigger_STAT_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Trigger_SYST_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Trigger_SYST_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_ID_STAT_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_ID_STAT_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_ID_SYST_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_ID_SYST_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Isol_STAT_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Isol_STAT_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Isol_SYST_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_Isol_SYST_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_TTVA_STAT_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_TTVA_STAT_DOWN;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_TTVA_SYST_UP;
    TTreeReaderValue<float> * m_weight_leptonSF_MU_SF_TTVA_SYST_DOWN;

    TTreeReaderValue<std::vector<float>> * m_weight_bTagSF_77_eigenvars_B_up;
    TTreeReaderValue<std::vector<float>> * m_weight_bTagSF_77_eigenvars_C_up;
    TTreeReaderValue<std::vector<float>> * m_weight_bTagSF_77_eigenvars_Light_up;
    TTreeReaderValue<std::vector<float>> * m_weight_bTagSF_77_eigenvars_B_down;
    TTreeReaderValue<std::vector<float>> * m_weight_bTagSF_77_eigenvars_C_down;
    TTreeReaderValue<std::vector<float>> * m_weight_bTagSF_77_eigenvars_Light_down;
    TTreeReaderValue<float> * m_weight_bTagSF_77_extrapolation_up;
    TTreeReaderValue<float> * m_weight_bTagSF_77_extrapolation_down;
    TTreeReaderValue<float> * m_weight_bTagSF_77_extrapolation_from_charm_up;
    TTreeReaderValue<float> * m_weight_bTagSF_77_extrapolation_from_charm_down;

    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_70_extrapolation_up;
    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_70_extrapolation_down;
    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_70_extrapolation_from_charm_up;
    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_70_extrapolation_from_charm_down;
    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_77_extrapolation_up;
    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_77_extrapolation_down;
    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_77_extrapolation_from_charm_up;
    TTreeReaderValue<float> * m_weight_trackjet_bTagSF_77_extrapolation_from_charm_down;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_70_eigenvars_B_up;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_70_eigenvars_C_up;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_70_eigenvars_Light_up;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_70_eigenvars_B_down;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_70_eigenvars_C_down;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_70_eigenvars_Light_down;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_77_eigenvars_B_up;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_77_eigenvars_C_up;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_77_eigenvars_Light_up;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_77_eigenvars_B_down;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_77_eigenvars_C_down;
    TTreeReaderValue<std::vector<float>> * m_weight_trackjet_bTagSF_77_eigenvars_Light_down;

    // event info 
    TTreeReaderValue<float> * m_eventNumber;
    TTreeReaderValue<float> * m_runNumber;
//    TTreeReaderValue<unsigned long long> * m_eventNumber;
//    TTreeReaderValue<unsigned int> * m_runNumber;
    TTreeReaderValue<unsigned int> * m_mcChannelNumber;
    TTreeReaderValue<float> * m_mu;
    TTreeReaderValue<int> * m_lumiblock;
    TTreeReaderValue<float> * m_treeXSection;
    TTreeReaderValue<float> * m_treeKFactor;
    TTreeReaderValue<float> * m_treeSumOfWeights;
    TTreeReaderValue<float> * m_met_met;
    TTreeReaderValue<float> * m_met_phi;

    // lepton info
    TTreeReaderValue<std::vector<float>> * m_el_pt;
    TTreeReaderValue<std::vector<float>> * m_el_eta;
    TTreeReaderValue<std::vector<float>> * m_el_phi;
    TTreeReaderValue<std::vector<float>> * m_el_e;
    TTreeReaderValue<std::vector<float>> * m_el_charge;
    TTreeReaderValue<std::vector<float>> * m_mu_pt;
    TTreeReaderValue<std::vector<float>> * m_mu_eta;
    TTreeReaderValue<std::vector<float>> * m_mu_phi;
    TTreeReaderValue<std::vector<float>> * m_mu_e;
    TTreeReaderValue<std::vector<float>> * m_mu_charge;

    // large-R jet info
    TTreeReaderValue<std::vector<float>> * m_ljet_pt;
    TTreeReaderValue<std::vector<float>> * m_ljet_eta;
    TTreeReaderValue<std::vector<float>> * m_ljet_phi;
    TTreeReaderValue<std::vector<float>> * m_ljet_e;
    TTreeReaderValue<std::vector<float>> * m_ljet_d23;
    TTreeReaderValue<std::vector<float>> * m_ljet_tau1_wta;
    TTreeReaderValue<std::vector<float>> * m_ljet_tau2_wta;
    TTreeReaderValue<std::vector<float>> * m_ljet_tau3_wta;
    TTreeReaderValue<std::vector<float>> * m_ljet_tau21_wta;
    TTreeReaderValue<std::vector<float>> * m_ljet_tau32_wta;
    TTreeReaderValue<std::vector<int>> * m_ljet_isGood;
    TTreeReaderValue<std::vector<float>> * m_ljet_charge;

    // truth large-R jet info
    TTreeReaderValue<std::vector<float>> * m_truth_ljet_pt;
    TTreeReaderValue<std::vector<float>> * m_truth_ljet_eta;
    TTreeReaderValue<std::vector<float>> * m_truth_ljet_phi;
    TTreeReaderValue<std::vector<float>> * m_truth_ljet_e;
    TTreeReaderValue<std::vector<float>> * m_truth_ljet_split23;
    TTreeReaderValue<std::vector<float>> * m_truth_ljet_Qw;
    TTreeReaderValue<std::vector<float>> * m_truth_ljet_tau32_wta;

    // jet info
    TTreeReaderValue<std::vector<float>> * m_jet_pt;
    TTreeReaderValue<std::vector<float>> * m_jet_eta;
    TTreeReaderValue<std::vector<float>> * m_jet_phi;
    TTreeReaderValue<std::vector<float>> * m_jet_e;
    TTreeReaderValue<std::vector<float>> * m_jet_mv2c10;
    TTreeReaderValue<std::vector<float>> * m_jet_mv2c20;
    TTreeReaderValue<std::vector<char>> * m_jet_isbtagged_77;
    TTreeReaderValue<std::vector<char>> * m_jet_isbtagged_70;
    TTreeReaderValue<std::vector<float>> * m_jet_jvt;
    TTreeReaderValue<std::vector<int>> * m_jet_true_flavor;

    // truth jet info
    TTreeReaderValue<std::vector<float>> * m_truth_jet_pt;
    TTreeReaderValue<std::vector<float>> * m_truth_jet_eta;
    TTreeReaderValue<std::vector<float>> * m_truth_jet_phi;
    TTreeReaderValue<std::vector<float>> * m_truth_jet_e;
    // flavor? other information?

    // truth information
    TTreeReaderValue<unsigned long long> * m_truthEventNumber;
    TTreeReaderValue<unsigned int> * m_truthRunNumber;
    TTreeReaderValue<float> * m_truth_weight_mc;
    TTreeReaderValue<float> * MC_b_from_t_pt;
    TTreeReaderValue<float> * MC_b_from_t_eta;
    TTreeReaderValue<float> * MC_b_from_t_phi;
    TTreeReaderValue<float> * MC_b_from_t_m;
    TTreeReaderValue<float> * MC_W_from_t_pt;
    TTreeReaderValue<float> * MC_W_from_t_eta;
    TTreeReaderValue<float> * MC_W_from_t_phi;
    TTreeReaderValue<float> * MC_W_from_t_m;
    TTreeReaderValue<float> * MC_b_from_tbar_pt;
    TTreeReaderValue<float> * MC_b_from_tbar_eta;
    TTreeReaderValue<float> * MC_b_from_tbar_phi;
    TTreeReaderValue<float> * MC_b_from_tbar_m;
    TTreeReaderValue<float> * MC_W_from_tbar_pt;
    TTreeReaderValue<float> * MC_W_from_tbar_eta;
    TTreeReaderValue<float> * MC_W_from_tbar_phi;
    TTreeReaderValue<float> * MC_W_from_tbar_m;
};
#endif

