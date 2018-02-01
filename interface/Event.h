#ifndef EVENT_H
#define EVENT_H

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

#include "Analysis/CyMiniAna/interface/physicsObjects.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/dileptonTtbarReco.h"
#include "lwtnn/lwtnn/interface/LightweightNeuralNetwork.hh"
#include "lwtnn/lwtnn/interface/parse_json.hh"


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

    virtual float met( const std::string& met_name );
    virtual float HT();
    virtual float ST();
    virtual void getBtaggedJets( Jet& jet );
    virtual std::vector<int> btag_jets(const std::string &wkpt);
    virtual std::vector<int> btag_jets(); // using configured b-tag WP

//    virtual unsigned long long eventNumber();
//    virtual unsigned int runNumber();
    long long entry() { return m_entry; }
    virtual float eventNumber();
    virtual float runNumber();
    virtual float mu();
    virtual int lumiblock();
    virtual std::string treeName();
    virtual float xsection();
    virtual float kfactor();
    virtual float sumOfWeights();
    void truth();

    // Neural network information
    void getDNNInputs();      // return the DNN inputs to the user
    void getDNN();            // get the DNN output
    double DNN();

    // kinematic reconstruction
    void buildTtbar();
    void getDilepton();
    std::map<std::string,LepTop> ttbar();

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

  protected:

    // general information
    configuration *m_config;
    TTreeReader &m_ttree;
    TTreeReader m_truth_tree;
    std::string m_treeName;
    std::string m_fileName;
    bool m_grid;
    bool m_isMC;
    long long m_entry;
    long long m_truth_entry;

    // 0/1/2-lepton analyses
    bool m_isZeroLeptonAnalysis;
    bool m_isOneLeptonAnalysis;
    bool m_isTwoLeptonAnalysis;

    // neural network & kinematic reconstruction
    bool m_kinematicReco;
    bool m_buildNeutrinos;
    bool m_getDNN;
    lwt::LightweightNeuralNetwork* m_lwnn;
    std::map<std::string, double> m_dnnInputs;   // values for inputs to the DNN
    std::string m_dnnKey;                        // string to access DNN from lwtnn
    double m_DNN;                                // DNN score

    // dilepton
    bool m_ee;
    bool m_mumu;
    bool m_emu;
    dileptonTtbarReco* m_dileptonTtbar;
    std::map<std::string,LepTop> m_ttbar;
    DileptonReco m_dilepton;

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

    // kinematics
    float m_HT;
    float m_ST;
    float m_metmet;
    float m_metphi;

    // nominal b-tagging weight maps
    std::map<std::string, float> m_weight_btag;
    float m_weight_btag_default;
    // Maps to keep track of weight systematics
    std::map<std::string,TTreeReaderValue<float> * > m_weightSystematicsFloats;
    std::map<std::string,TTreeReaderArray<float> * > m_weightSystematicsVectorFloats;
    std::vector<std::string> m_listOfWeightSystematics;




    // TTree variables [all possible ones]
    // ************
    // the following are from root files accessed 
    //  on 6 November 2017 (N. Eminizer)
    //  /store/user/eminizer/B2GTTrees_Jun19/
    unsigned int m_array_size = 500;

    TTreeReaderValue<float> * m_weight_mc;          // evt_Gen_Weight
    TTreeReaderValue<float> * m_weight_pileup;
    TTreeReaderValue<float> * m_weight_lept_eff;
    TTreeReaderValue<float> * m_weight_pileup_UP;
    TTreeReaderValue<float> * m_weight_pileup_DOWN;

    // Event info 
    TTreeReaderValue<int> * m_eventNumber;          // evt_EventNumber
    TTreeReaderValue<int> * m_runNumber;            // evt_RunNumber
    TTreeReaderValue<float> * m_mu;
    TTreeReaderValue<double> * m_rho;               // evt_rho
    TTreeReaderValue<int> * m_lumiblock;            // evt_LumiBlock
    TTreeReaderValue<float> * m_treeXSection;       // evt_XSec
    TTreeReaderValue<float> * m_treeKFactor;
    TTreeReaderValue<float> * m_treeSumOfWeights;
    TTreeReaderArray<int> * m_NGoodVtx;             // evt_NGoodVtx
    TTreeReaderArray<int> * m_LHAPDF_ID;            // evt_LHA_PDF_ID
    TTreeReaderArray<int> * m_NIsoTrk;              // evt_NIsoTrk
    TTreeReaderArray<int> * m_pu_NtrueInt;          // pu_NtrueInt

    TTreeReaderValue<unsigned int> * m_scale_size;  // scale_size
    TTreeReaderArray<float> * m_scale_Weights;      // scale_Weights
    TTreeReaderValue<unsigned int> * m_pdf_size;    // pdf_size
    TTreeReaderArray<float> * m_weights_pdf;        // pdf_Weights
    TTreeReaderValue<unsigned int> * m_alphas_size; // alphas_size
    TTreeReaderArray<float> * m_weights_alphas;     // alphas_Weights

    TTreeReaderArray<int> * m_Flag_BadPFMuonFilter;                    // Flag_BadPFMuonFilter
    TTreeReaderArray<int> * m_Flag_BadChargedCandidateFilter;          // Flag_BadChargedCandidateFilter
    TTreeReaderArray<int> * m_Flag_HBHENoiseFilter;                    // Flag_HBHENoiseFilter
    TTreeReaderArray<int> * m_Flag_HBHENoiseIsoFilter;                 // Flag_HBHENoiseIsoFilter
    TTreeReaderArray<int> * m_Flag_CSCTightHaloFilter;                 // Flag_CSCTightHaloFilter
    TTreeReaderArray<int> * m_Flag_CSCTightHaloTrkMuUnvetoFilter;      // Flag_CSCTightHaloTrkMuUnvetoFilter
    TTreeReaderArray<int> * m_Flag_CSCTightHalo2015Filter;             // Flag_CSCTightHalo2015Filter
    TTreeReaderArray<int> * m_Flag_globalTightHalo2016Filter;          // Flag_globalTightHalo2016Filter
    TTreeReaderArray<int> * m_Flag_globalSuperTightHalo2016Filter;     // Flag_globalSuperTightHalo2016Filter
    TTreeReaderArray<int> * m_Flag_HcalStripHaloFilter;                // Flag_HcalStripHaloFilter
    TTreeReaderArray<int> * m_Flag_hcalLaserEventFilter;               // Flag_hcalLaserEventFilter
    TTreeReaderArray<int> * m_Flag_EcalDeadCellTriggerPrimitiveFilter; // Flag_EcalDeadCellTriggerPrimitiveFilter
    TTreeReaderArray<int> * m_Flag_EcalDeadCellBoundaryEnergyFilter;   // Flag_EcalDeadCellBoundaryEnergyFilter
    TTreeReaderArray<int> * m_Flag_goodVertices;                       // Flag_goodVertices
    TTreeReaderArray<int> * m_Flag_eeBadScFilter;                      // Flag_eeBadScFilter
    TTreeReaderArray<int> * m_Flag_ecalLaserCorrFilter;                // Flag_ecalLaserCorrFilter
    TTreeReaderArray<int> * m_Flag_trkPOGFilters;                      // Flag_trkPOGFilters
    TTreeReaderArray<int> * m_Flag_chargedHadronTrackResolutionFilter; // Flag_chargedHadronTrackResolutionFilter
    TTreeReaderArray<int> * m_Flag_muonBadTrackFilter;                 // Flag_muonBadTrackFilter
    TTreeReaderArray<int> * m_Flag_trkPOG_manystripclus53X;            // Flag_trkPOG_manystripclus53X
    TTreeReaderArray<int> * m_Flag_trkPOG_toomanystripclus53X;         // Flag_trkPOG_toomanystripclus53X
    TTreeReaderArray<int> * m_Flag_trkPOG_logErrorTooManyClusters;     // Flag_trkPOG_logErrorTooManyClusters
    TTreeReaderArray<int> * m_Flag_METFilters;                         // Flag_METFilters
    TTreeReaderArray<int> * m_Flag_badMuons;                           // Flag_badMuons
    TTreeReaderArray<int> * m_Flag_duplicateMuons;                     // Flag_duplicateMuons
    TTreeReaderArray<int> * m_Flag_noBadMuons;                         // Flag_noBadMuons



    // MET
    TTreeReaderValue<unsigned int> * m_met_size;  // met_size
    TTreeReaderArray<float> * m_met_met;          // met_Pt
    TTreeReaderArray<float> * m_met_phi;          // met_Phi
    TTreeReaderArray<float> * m_met_met_uncor;    // met_uncorPt
    TTreeReaderArray<float> * m_met_phi_uncor;    // met_uncorPhi

    TTreeReaderValue<unsigned int> * m_met_muCleanOnly_size;  // met_MuCleanOnly_size
    TTreeReaderArray<float> * m_met_muCleanOnly_met;          // met_MuCleanOnly_Pt
    TTreeReaderArray<float> * m_met_muCleanOnly_phi;          // met_MuCleanOnly_Phi
    TTreeReaderArray<float> * m_met_muCleanOnly_met_uncor;    // met_MuCleanOnly_uncorPt
    TTreeReaderArray<float> * m_met_muCleanOnly_phi_uncor;    // met_MuCleanOnly_uncorPhi

    TTreeReaderValue<unsigned int> * m_met_syst_size;  // metsyst_size
    TTreeReaderArray<float> * m_met_syst_met;          // metsyst_Pt
    TTreeReaderArray<float> * m_met_syst_phi;          // metsyst_Phi
    TTreeReaderArray<float> * m_met_syst_muCleanOnly_met;    // metsyst_MuCleanOnly_Pt
    TTreeReaderArray<float> * m_met_syst_muCleanOnly_phi;    // metsyst_MuCleanOnly_Phi

    // Leptons
    TTreeReaderValue<unsigned int> * m_el_size;    // el_size
    TTreeReaderArray<float> * m_el_pt;             // el_Pt
    TTreeReaderArray<float> * m_el_eta;            // el_Eta
    TTreeReaderArray<float> * m_el_phi;            // el_Phi
    TTreeReaderArray<float> * m_el_e;              // el_E
    TTreeReaderArray<float> * m_el_charge;         // el_Charge
    TTreeReaderArray<float> * m_el_key;            // el_Key
    TTreeReaderArray<float> * m_el_iso03;          // el_Iso03
    TTreeReaderArray<float> * m_el_iso03db;        // el_Iso03db
    TTreeReaderArray<float> * m_el_miniIso;        // el_MiniIso
    TTreeReaderArray<float> * m_el_SCEta;          // el_SCEta
    TTreeReaderArray<float> * m_el_SCPhi;          // el_SCPhi
    TTreeReaderArray<float> * m_el_vidVeto;        // el_vidVeto
    TTreeReaderArray<float> * m_el_vidLoose;       // el_vidLoose
    TTreeReaderArray<float> * m_el_vidMedium;      // el_vidMedium
    TTreeReaderArray<float> * m_el_vidTight;       // el_vidTight
    TTreeReaderArray<float> * m_el_vidHEEP;        // el_vidHEEP
    TTreeReaderArray<float> * m_el_vidVetonoiso;   // el_vidVetonoiso
    TTreeReaderArray<float> * m_el_vidLoosenoiso;  // el_vidLoosenoiso
    TTreeReaderArray<float> * m_el_vidMediumnoiso; // el_vidMediumnoiso
    TTreeReaderArray<float> * m_el_vidTightnoiso;  // el_vidTightnoiso
    TTreeReaderArray<float> * m_el_vidHEEPnoiso;   // el_vidHEEPnoiso
    TTreeReaderArray<float> * m_el_vidMvaGPvalue;  // el_vidMvaGPvalue
    TTreeReaderArray<float> * m_el_vidMvaGPcateg;  // el_vidMvaGPcateg
    TTreeReaderArray<float> * m_el_vidMvaHZZvalue; // el_vidMvaHZZvalue
    TTreeReaderArray<float> * m_el_vidMvaHZZcateg; // el_vidMvaHZZcateg
    TTreeReaderArray<int> * m_el_veto_NoIsoID;     // el_IDVeto_NoIso
    TTreeReaderArray<int> * m_el_loose_NoIsoID;    // el_IDLoose_NoIso
    TTreeReaderArray<int> * m_el_medium_NoIsoID;   // el_IDMedium_NoIso
    TTreeReaderArray<int> * m_el_tight_NoIsoID;    // el_IDTight_NoIso
    TTreeReaderArray<int> * m_el_isoVeto;          // el_IsoVeto
    TTreeReaderArray<int> * m_el_isoLoose;         // el_IsoLoose
    TTreeReaderArray<int> * m_el_isoMedium;        // el_IsoMedium
    TTreeReaderArray<int> * m_el_isoTight;         // el_IsoTight
    TTreeReaderArray<int> * m_el_vetoID;           // el_IDVeto
    TTreeReaderArray<int> * m_el_looseID;          // el_IDLoose
    TTreeReaderArray<int> * m_el_mediumID;         // el_IDMedium
    TTreeReaderArray<int> * m_el_tightID;          // el_IDTight


    TTreeReaderValue<int> * m_mu_size;         // mu_size
    TTreeReaderArray<float> * m_mu_pt;         // mu_Pt
    TTreeReaderArray<float> * m_mu_eta;        // mu_Eta
    TTreeReaderArray<float> * m_mu_phi;        // mu_Phi
    TTreeReaderArray<float> * m_mu_e;          // mu_E
    TTreeReaderArray<float> * m_mu_charge;     // mu_Charge
    TTreeReaderArray<float> * m_mu_key;        // mu_Key
    TTreeReaderArray<float> * m_mu_iso04;      // mu_Iso04
    TTreeReaderArray<float> * m_mu_miniIso;    // mu_MiniIso
    TTreeReaderArray<float> * m_mu_soft;       // mu_IsSoftMuon
    TTreeReaderArray<float> * m_mu_loose;      // mu_IsLooseMuon
    TTreeReaderArray<float> * m_mu_medium;     // mu_IsMediumMuon
    TTreeReaderArray<float> * m_mu_medium2016; // mu_IsMediumMuon2016
    TTreeReaderArray<float> * m_mu_tight;      // mu_IsTightMuon
    TTreeReaderArray<float> * m_mu_hightPt;    // mu_IsHighPtMuon

    // large-R jet info
    TTreeReaderValue<float> * m_dnn_score;

    TTreeReaderValue<unsigned int> * m_ljet_size; // jetAK8CHS_size
    TTreeReaderArray<float> * m_ljet_pt;          // jetAK8CHS_Pt
    TTreeReaderArray<float> * m_ljet_eta;         // jetAK8CHS_Eta
    TTreeReaderArray<float> * m_ljet_phi;         // jetAK8CHS_Phi
    TTreeReaderArray<float> * m_ljet_e;           // jetAK8CHS_E
    TTreeReaderArray<float> * m_ljet_tau1_CHS;    // jetAK8CHS_tau1CHS
    TTreeReaderArray<float> * m_ljet_tau2_CHS;    // jetAK8CHS_tau1CHS
    TTreeReaderArray<float> * m_ljet_tau3_CHS;    // jetAK8CHS_tau1CHS
    TTreeReaderArray<float> * m_ljet_charge;      // jetAK8CHS_Charge
    TTreeReaderArray<float> * m_ljet_softDropMass_CHS;  // jetAK8CHS_softDropMassCHS

    TTreeReaderArray<float> * m_ljet_CSVv2;                   // jetAK8CHS_CSVv2
    TTreeReaderArray<float> * m_ljet_CMVAv2;                  // jetAK8CHS_CMVAv2
    TTreeReaderArray<float> * m_ljet_CvsL;                    // jetAK8CHS_CvsL
    TTreeReaderArray<float> * m_ljet_CvsB;                    // jetAK8CHS_CvsB
    TTreeReaderArray<float> * m_ljet_PartonFlavour;           // jetAK8CHS_PartonFlavour
    TTreeReaderArray<float> * m_ljet_HadronFlavour;           // jetAK8CHS_HadronFlavour
    TTreeReaderArray<float> * m_ljet_neutralMultiplicity;     // jetAK8CHS_neutralMultiplicity
    TTreeReaderArray<float> * m_ljet_neutralHadronEnergyFrac; // jetAK8CHS_neutralHadronEnergyFrac
    TTreeReaderArray<float> * m_ljet_neutralEmEnergyFrac;     // jetAK8CHS_neutralEmEnergyFrac
    TTreeReaderArray<float> * m_ljet_chargedHadronEnergyFrac; // jetAK8CHS_chargedHadronEnergyFrac
    TTreeReaderArray<float> * m_ljet_chargedEmEnergyFrac;     // jetAK8CHS_chargedEmEnergyFrac
    TTreeReaderArray<float> * m_ljet_chargedMultiplicity;     // jetAK8CHS_chargedMultiplicity
    TTreeReaderArray<float> * m_ljet_jecFactor0;              // jetAK8CHS_jecFactor0
    TTreeReaderArray<float> * m_ljet_jetArea;                 // jetAK8CHS_jetArea
    TTreeReaderArray<float> * m_ljet_jecUncertainty;          // jetAK8CHS_jecUncertainty
    TTreeReaderArray<float> * m_ljet_PtResolution;            // jetAK8CHS_PtResolution
    TTreeReaderArray<float> * m_ljet_JERSF;                   // jetAK8CHS_JERSF
    TTreeReaderArray<float> * m_ljet_JERSFUp;                 // jetAK8CHS_JERSFUp
    TTreeReaderArray<float> * m_ljet_JERSFDown;               // jetAK8CHS_JERSFDown
    TTreeReaderArray<float> * m_ljet_SmearedPt;               // jetAK8CHS_SmearedPt
    TTreeReaderArray<float> * m_ljet_vSubjetIndex0;           // jetAK8CHS_vSubjetIndex0
    TTreeReaderArray<float> * m_ljet_vSubjetIndex1;           // jetAK8CHS_vSubjetIndex1
    std::vector<std::vector<int> > * m_ljet_keys;             // jetAK8CHS_Keys

    TTreeReaderValue<unsigned int> * m_ljet_subjet_size;             // subjetAK8CHS_size
    TTreeReaderArray<float> * m_ljet_subjet_pt;                      // subjetAK8CHS_Pt
    TTreeReaderArray<float> * m_ljet_subjet_eta;                     // subjetAK8CHS_Eta
    TTreeReaderArray<float> * m_ljet_subjet_phi;                     // subjetAK8CHS_Phi
    TTreeReaderArray<float> * m_ljet_subjet_e;                       // subjetAK8CHS_E
    TTreeReaderArray<float> * m_ljet_subjet_charge;                  // subjetAK8CHS_Charge
    TTreeReaderArray<float> * m_ljet_subjet_CSVv2;                   // subjetAK8CHS_CSVv2
    TTreeReaderArray<float> * m_ljet_subjet_CMVAv2;                  // subjetAK8CHS_CMVAv2
    TTreeReaderArray<float> * m_ljet_subjet_CvsL;                    // subjetAK8CHS_CvsL
    TTreeReaderArray<float> * m_ljet_subjet_CvsB;                    // subjetAK8CHS_CvsB
    TTreeReaderArray<float> * m_ljet_subjet_partonFlavour;           // subjetAK8CHS_PartonFlavour
    TTreeReaderArray<float> * m_ljet_subjet_hadronFlavour;           // subjetAK8CHS_HadronFlavour
    TTreeReaderArray<float> * m_ljet_subjet_neutralMultiplicity;     // subjetAK8CHS_neutralMultiplicity
    TTreeReaderArray<float> * m_ljet_subjet_neutralHadronEnergyFrac; // subjetAK8CHS_neutralHadronEnergyFrac
    TTreeReaderArray<float> * m_ljet_subjet_neutralEmEnergyFrac;     // subjetAK8CHS_neutralEmEnergyFrac
    TTreeReaderArray<float> * m_ljet_subjet_chargedHadronEnergyFrac; // subjetAK8CHS_chargedHadronEnergyFrac
    TTreeReaderArray<float> * m_ljet_subjet_chargedEmEnergyFrac;     // subjetAK8CHS_chargedEmEnergyFrac
    TTreeReaderArray<float> * m_ljet_subjet_chargedMultiplicity;     // subjetAK8CHS_chargedMultiplicity
    TTreeReaderArray<float> * m_ljet_subjet_jecFactor0;              // subjetAK8CHS_jecFactor0
    TTreeReaderArray<float> * m_ljet_subjet_jetArea;                 // subjetAK8CHS_jetArea
    std::vector<std::vector<int> > * m_ljet_subjet_keys;             // subjetAK8CHS_Keys

    // truth large-R jet info
    TTreeReaderArray<float> * m_truth_ljet_pt;            // jetAK8CHS_GenJetPt
    TTreeReaderArray<float> * m_truth_ljet_eta;           // jetAK8CHS_GenJetEta
    TTreeReaderArray<float> * m_truth_ljet_phi;           // jetAK8CHS_GenJetPhi
    TTreeReaderArray<float> * m_truth_ljet_e;             // jetAK8CHS_GenJetE
    TTreeReaderArray<float> * m_truth_ljet_charge;        // jetAK8CHS_GenJetCharge
    TTreeReaderArray<float> * m_truth_ljet_subjet_pt;     // subjetAK8CHS_GenJetPt
    TTreeReaderArray<float> * m_truth_ljet_subjet_eta;    // subjetAK8CHS_GenJetEta
    TTreeReaderArray<float> * m_truth_ljet_subjet_phi;    // subjetAK8CHS_GenJetPhi
    TTreeReaderArray<float> * m_truth_ljet_subjet_e;      // subjetAK8CHS_GenJetE
    TTreeReaderArray<float> * m_truth_ljet_subjet_charge; // subjetAK8CHS_GenJetCharge
    TTreeReaderArray<float> * m_truth_ljet_Qw;
    TTreeReaderArray<float> * m_truth_ljet_tau32_wta;
    TTreeReaderArray<float> * m_truth_ljet_split23;

    // Jet info
    TTreeReaderValue<unsigned int> * m_jet_size;             // jetAK4CHS_size
    TTreeReaderArray<float> * m_jet_pt;                      // jetAK4CHS_Pt
    TTreeReaderArray<float> * m_jet_eta;                     // jetAK4CHS_Eta
    TTreeReaderArray<float> * m_jet_phi;                     // jetAK4CHS_Phi
    TTreeReaderArray<float> * m_jet_e;                       // jetAK4CHS_E
    TTreeReaderArray<float> * m_jet_charge;                  // jetAK4CHS_Charge
    TTreeReaderArray<float> * m_jet_CSVv2;                   // jetAK4CHS_CSVv2
    TTreeReaderArray<float> * m_jet_CMVAv2;                  // jetAK4CHS_CMVAv2
    TTreeReaderArray<float> * m_jet_CvsL;                    // jetAK4CHS_CvsL
    TTreeReaderArray<float> * m_jet_CvsB;                    // jetAK4CHS_CvsB
    TTreeReaderArray<float> * m_jet_partonFlavour;           // jetAK4CHS_PartonFlavour
    TTreeReaderArray<float> * m_jet_hadronFlavour;           // jetAK4CHS_HadronFlavour
    TTreeReaderArray<float> * m_jet_neutralMultiplicity;     // jetAK4CHS_neutralMultiplicity
    TTreeReaderArray<float> * m_jet_neutralHadronEnergyFrac; // jetAK4CHS_neutralHadronEnergyFrac
    TTreeReaderArray<float> * m_jet_neutralEmEnergyFrac;     // jetAK4CHS_neutralEmEnergyFrac
    TTreeReaderArray<float> * m_jet_chargedHadronEnergyFrac; // jetAK4CHS_chargedHadronEnergyFrac
    TTreeReaderArray<float> * m_jet_chargedEmEnergyFrac;     // jetAK4CHS_chargedEmEnergyFrac
    TTreeReaderArray<float> * m_jet_chargedMultiplicity;     // jetAK4CHS_chargedMultiplicity
    TTreeReaderArray<float> * m_jet_jecFactor0;              // jetAK4CHS_jecFactor0
    TTreeReaderArray<float> * m_jet_jetArea;                 // jetAK4CHS_jetArea
    TTreeReaderArray<float> * m_jet_jecUncertainty;          // jetAK4CHS_jecUncertainty
    TTreeReaderArray<float> * m_jet_ptResolution;            // jetAK4CHS_PtResolution
    TTreeReaderArray<float> * m_jet_JERSF;                   // jetAK4CHS_JERSF
    TTreeReaderArray<float> * m_jet_JERSFUp;                 // jetAK4CHS_JERSFUp
    TTreeReaderArray<float> * m_jet_JERSFDown;               // jetAK4CHS_JERSFDown
    TTreeReaderArray<float> * m_jet_smearedPt;               // jetAK4CHS_SmearedPt
    std::vector<std::vector<int> > * m_jet_keys;             // jetAK4CHS_Keys



    // Reconstructed neutrinos
    TTreeReaderArray<float> * m_nu_pt;
    TTreeReaderArray<float> * m_nu_eta;
    TTreeReaderArray<float> * m_nu_phi;

    // Reconstructed ttbar
    TTreeReaderValue<float> * m_top_pt;
    TTreeReaderValue<float> * m_top_eta;
    TTreeReaderValue<float> * m_top_phi;
    TTreeReaderValue<float> * m_top_e;
    TTreeReaderValue<int> * m_lepton_top_index;
    TTreeReaderValue<int> * m_jet_top_index;
    TTreeReaderValue<int> * m_nu_top_index;
    TTreeReaderValue<float> * m_antitop_pt;
    TTreeReaderValue<float> * m_antitop_eta;
    TTreeReaderValue<float> * m_antitop_phi;
    TTreeReaderValue<float> * m_antitop_e;
    TTreeReaderValue<int> * m_lepton_antitop_index;
    TTreeReaderValue<int> * m_jet_antitop_index;
    TTreeReaderValue<int> * m_nu_antitop_index;

    TTreeReaderValue<float> * m_dileptonTtbarWeight;
    TTreeReaderValue<float> * m_semileptonTtbarWeight;
    TTreeReaderValue<float> * m_allhadTtbarWeight;


    // Truth jet info
    TTreeReaderArray<float> * m_truth_jet_pt;           // jetAK4CHS_GenJetPt
    TTreeReaderArray<float> * m_truth_jet_eta;          // jetAK4CHS_GenJetEta
    TTreeReaderArray<float> * m_truth_jet_phi;          // jetAK4CHS_GenJetPhi
    TTreeReaderArray<float> * m_truth_jet_e;            // jetAK4CHS_GenJetCharge

    // General truth information
    TTreeReaderValue<unsigned long long> * m_truthEventNumber;
    TTreeReaderValue<unsigned int> * m_truthRunNumber;
    TTreeReaderValue<float> * m_truth_weight_mc;

    TTreeReaderValue<float> * m_mc_ht;                // evt_Gen_Ht

    TTreeReaderArray<float> * m_MC_part1_factor;      // MC_part1_factor
    TTreeReaderArray<float> * m_MC_part1_ID;          // MC_part1_ID
    TTreeReaderArray<float> * m_MC_part2_factor;      // MC_part2_factor
    TTreeReaderArray<float> * m_MC_part2_ID;          // MC_part2_ID
    TTreeReaderArray<float> * m_MC_t_pt;              // MC_t_pt
    TTreeReaderArray<float> * m_MC_t_eta;             // MC_t_eta
    TTreeReaderArray<float> * m_MC_t_phi;             // MC_t_phi
    TTreeReaderArray<float> * m_MC_t_e;               // MC_t_E
    TTreeReaderArray<float> * m_MC_tbar_pt;           // MC_tbar_pt
    TTreeReaderArray<float> * m_MC_tbar_eta;          // MC_tbar_eta
    TTreeReaderArray<float> * m_MC_tbar_phi;          // MC_tbar_phi
    TTreeReaderArray<float> * m_MC_tbar_e;            // MC_tbar_E
    TTreeReaderArray<float> * m_MC_lep_pt;            // MC_lep_pt
    TTreeReaderArray<float> * m_MC_lep_eta;           // MC_lep_eta
    TTreeReaderArray<float> * m_MC_lep_phi;           // MC_lep_phi
    TTreeReaderArray<float> * m_MC_lep_e;             // MC_lep_E
    TTreeReaderArray<float> * m_MC_lep_ID;            // MC_lep_ID
    TTreeReaderArray<float> * m_MC_nu_pt;             // MC_nu_pt
    TTreeReaderArray<float> * m_MC_nu_eta;            // MC_nu_eta
    TTreeReaderArray<float> * m_MC_nu_phi;            // MC_nu_phi
    TTreeReaderArray<float> * m_MC_nu_e;              // MC_nu_E
    TTreeReaderArray<float> * m_MC_lepb_pt;           // MC_lepb_pt
    TTreeReaderArray<float> * m_MC_lepb_eta;          // MC_lepb_eta
    TTreeReaderArray<float> * m_MC_lepb_phi;          // MC_lepb_phi
    TTreeReaderArray<float> * m_MC_lepb_e;            // MC_lepb_E
    TTreeReaderArray<float> * m_MC_hadW_pt;           // MC_hadW_pt
    TTreeReaderArray<float> * m_MC_hadW_eta;          // MC_hadW_eta
    TTreeReaderArray<float> * m_MC_hadW_phi;          // MC_hadW_phi
    TTreeReaderArray<float> * m_MC_hadW_e;            // MC_hadW_E
    TTreeReaderArray<float> * m_MC_hadb_pt;           // MC_hadb_pt
    TTreeReaderArray<float> * m_MC_hadb_eta;          // MC_hadb_eta
    TTreeReaderArray<float> * m_MC_hadb_phi;          // MC_hadb_phi
    TTreeReaderArray<float> * m_MC_hadb_e;            // MC_hadb_E
    TTreeReaderArray<float> * m_MC_cstar;             // MC_cstar
    TTreeReaderArray<float> * m_MC_x_F;               // MC_x_F
    TTreeReaderArray<float> * m_MC_Mtt;               // MC_Mtt

    // HLT 
    TTreeReaderArray<int> * m_HLT_Ele45_WPLoose_Gsf;          // HLT_Ele45_WPLoose_Gsf
    TTreeReaderArray<int> * m_HLT_Ele45_WPLoose_Gsf_prescale; // HLT_Ele45_WPLoose_Gsf_prescale

    TTreeReaderArray<int> * m_HLT_Mu50;            // HLT_Mu50
    TTreeReaderArray<int> * m_HLT_Mu50_prescale;   // HLT_Mu50_prescale
    TTreeReaderArray<int> * m_HLT_TkMu50;          // HLT_TkMu50
    TTreeReaderArray<int> * m_HLT_TkMu50_prescale; // HLT_TkMu50_prescale

//   TTreeReaderArray<int> * m_HLT_Mu17;          // HLT_Mu17
//   TTreeReaderArray<int> * m_HLT_Mu17_prescale;          // HLT_Mu17_prescale
//   TTreeReaderArray<int> * m_HLT_Mu20;          // HLT_Mu20
//   TTreeReaderArray<int> * m_HLT_Mu20_prescale;          // HLT_Mu20_prescale
//   TTreeReaderArray<int> * m_HLT_Mu27;          // HLT_Mu27
//   TTreeReaderArray<int> * m_HLT_Mu27_prescale;          // HLT_Mu27_prescale
//   TTreeReaderArray<int> * m_HLT_Mu55;          // HLT_Mu55
//   TTreeReaderArray<int> * m_HLT_Mu55_prescale;          // HLT_Mu55_prescale
//   TTreeReaderArray<int> * m_HLT_TkMu17;          // HLT_TkMu17
//   TTreeReaderArray<int> * m_HLT_TkMu17_prescale;          // HLT_TkMu17_prescale
//   TTreeReaderArray<int> * m_HLT_TkMu20;          // HLT_TkMu20
//   TTreeReaderArray<int> * m_HLT_TkMu20_prescale;          // HLT_TkMu20_prescale
//   TTreeReaderArray<int> * m_HLT_TkMu27;          // HLT_TkMu27
//   TTreeReaderArray<int> * m_HLT_TkMu27_prescale;          // HLT_TkMu27_prescale
//   TTreeReaderArray<int> * m_HLT_IsoMu18;          // HLT_IsoMu18
//   TTreeReaderArray<int> * m_HLT_IsoMu18_prescale;          // HLT_IsoMu18_prescale
//   TTreeReaderArray<int> * m_HLT_IsoMu20;          // HLT_IsoMu20
//   TTreeReaderArray<int> * m_HLT_IsoMu20_prescale;          // HLT_IsoMu20_prescale
//   TTreeReaderArray<int> * m_HLT_IsoMu22;          // HLT_IsoMu22
//   TTreeReaderArray<int> * m_HLT_IsoMu22_prescale;          // HLT_IsoMu22_prescale
//   TTreeReaderArray<int> * m_HLT_IsoMu24;          // HLT_IsoMu24
//   TTreeReaderArray<int> * m_HLT_IsoMu24_prescale;          // HLT_IsoMu24_prescale
//   TTreeReaderArray<int> * m_HLT_IsoMu27;          // HLT_IsoMu27
//   TTreeReaderArray<int> * m_HLT_IsoMu27_prescale;          // HLT_IsoMu27_prescale
//   TTreeReaderArray<int> * m_HLT_IsoTkMu18;          // HLT_IsoTkMu18
//   TTreeReaderArray<int> * m_HLT_IsoTkMu18_prescale;          // HLT_IsoTkMu18_prescale
//   TTreeReaderArray<int> * m_HLT_IsoTkMu20;          // HLT_IsoTkMu20
//   TTreeReaderArray<int> * m_HLT_IsoTkMu20_prescale;          // HLT_IsoTkMu20_prescale
//   TTreeReaderArray<int> * m_HLT_IsoTkMu22;          // HLT_IsoTkMu22
//   TTreeReaderArray<int> * m_HLT_IsoTkMu22_prescale;          // HLT_IsoTkMu22_prescale
//   TTreeReaderArray<int> * m_HLT_IsoTkMu24;          // HLT_IsoTkMu24
//   TTreeReaderArray<int> * m_HLT_IsoTkMu24_prescale;          // HLT_IsoTkMu24_prescale
//   TTreeReaderArray<int> * m_HLT_IsoTkMu27;          // HLT_IsoTkMu27
//   TTreeReaderArray<int> * m_HLT_IsoTkMu27_prescale;          // HLT_IsoTkMu27_prescale
//   TTreeReaderArray<int> * m_HLT_Ele17_CaloIdL_GsfTrkIdVL;          // HLT_Ele17_CaloIdL_GsfTrkIdVL
//   TTreeReaderArray<int> * m_HLT_Ele17_CaloIdL_GsfTrkIdVL_prescale;          // HLT_Ele17_CaloIdL_GsfTrkIdVL_prescale
//   TTreeReaderArray<int> * m_HLT_Ele22_eta2p1_WPLoose_Gsf;          // HLT_Ele22_eta2p1_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele22_eta2p1_WPLoose_Gsf_prescale;          // HLT_Ele22_eta2p1_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele23_WPLoose_Gsf;          // HLT_Ele23_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele23_WPLoose_Gsf_prescale;          // HLT_Ele23_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele24_eta2p1_WPLoose_Gsf;          // HLT_Ele24_eta2p1_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele24_eta2p1_WPLoose_Gsf_prescale;          // HLT_Ele24_eta2p1_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele25_WPTight_Gsf;          // HLT_Ele25_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele25_WPTight_Gsf_prescale;          // HLT_Ele25_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele25_eta2p1_WPLoose_Gsf;          // HLT_Ele25_eta2p1_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele25_eta2p1_WPLoose_Gsf_prescale;          // HLT_Ele25_eta2p1_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele25_eta2p1_WPTight_Gsf;          // HLT_Ele25_eta2p1_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele25_eta2p1_WPTight_Gsf_prescale;          // HLT_Ele25_eta2p1_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele27_WPLoose_Gsf;          // HLT_Ele27_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele27_WPLoose_Gsf_prescale;          // HLT_Ele27_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele27_WPTight_Gsf;          // HLT_Ele27_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele27_WPTight_Gsf_prescale;          // HLT_Ele27_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele27_eta2p1_WPLoose_Gsf;          // HLT_Ele27_eta2p1_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele27_eta2p1_WPLoose_Gsf_prescale;          // HLT_Ele27_eta2p1_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele27_eta2p1_WPTight_Gsf;          // HLT_Ele27_eta2p1_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele27_eta2p1_WPTight_Gsf_prescale;          // HLT_Ele27_eta2p1_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele30_WPTight_Gsf;          // HLT_Ele30_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele30_WPTight_Gsf_prescale;          // HLT_Ele30_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele30_eta2p1_WPLoose_Gsf;          // HLT_Ele30_eta2p1_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele30_eta2p1_WPLoose_Gsf_prescale;          // HLT_Ele30_eta2p1_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele30_eta2p1_WPTight_Gsf;          // HLT_Ele30_eta2p1_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele30_eta2p1_WPTight_Gsf_prescale;          // HLT_Ele30_eta2p1_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele32_WPTight_Gsf;          // HLT_Ele32_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele32_WPTight_Gsf_prescale;          // HLT_Ele32_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele32_eta2p1_WPLoose_Gsf;          // HLT_Ele32_eta2p1_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele32_eta2p1_WPLoose_Gsf_prescale;          // HLT_Ele32_eta2p1_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele32_eta2p1_WPTight_Gsf;          // HLT_Ele32_eta2p1_WPTight_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele32_eta2p1_WPTight_Gsf_prescale;          // HLT_Ele32_eta2p1_WPTight_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele35_WPLoose_Gsf;          // HLT_Ele35_WPLoose_Gsf
//   TTreeReaderArray<int> * m_HLT_Ele35_WPLoose_Gsf_prescale;          // HLT_Ele35_WPLoose_Gsf_prescale
//   TTreeReaderArray<int> * m_HLT_Ele105_CaloIdVT_GsfTrkIdT;          // HLT_Ele105_CaloIdVT_GsfTrkIdT
//   TTreeReaderArray<int> * m_HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale;          // HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale
//   TTreeReaderArray<int> * m_HLT_Ele115_CaloIdVT_GsfTrkIdT;          // HLT_Ele115_CaloIdVT_GsfTrkIdT
//   TTreeReaderArray<int> * m_HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale;          // HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale
//   TTreeReaderArray<int> * m_HLT_Ele145_CaloIdVT_GsfTrkIdT;          // HLT_Ele145_CaloIdVT_GsfTrkIdT
//   TTreeReaderArray<int> * m_HLT_Ele145_CaloIdVT_GsfTrkIdT_prescale;          // HLT_Ele145_CaloIdVT_GsfTrkIdT_prescale
//   TTreeReaderArray<int> * m_HLT_Ele200_CaloIdVT_GsfTrkIdT;          // HLT_Ele200_CaloIdVT_GsfTrkIdT
//   TTreeReaderArray<int> * m_HLT_Ele200_CaloIdVT_GsfTrkIdT_prescale;          // HLT_Ele200_CaloIdVT_GsfTrkIdT_prescale
//   TTreeReaderArray<int> * m_HLT_Ele250_CaloIdVT_GsfTrkIdT;          // HLT_Ele250_CaloIdVT_GsfTrkIdT
//   TTreeReaderArray<int> * m_HLT_Ele250_CaloIdVT_GsfTrkIdT_prescale;          // HLT_Ele250_CaloIdVT_GsfTrkIdT_prescale
//   TTreeReaderArray<int> * m_HLT_Ele300_CaloIdVT_GsfTrkIdT;          // HLT_Ele300_CaloIdVT_GsfTrkIdT
//   TTreeReaderArray<int> * m_HLT_Ele300_CaloIdVT_GsfTrkIdT_prescale;          // HLT_Ele300_CaloIdVT_GsfTrkIdT_prescale
//   TTreeReaderArray<int> * m_HLT_Mu30_eta2p1_PFJet150_PFJet50;          // HLT_Mu30_eta2p1_PFJet150_PFJet50
//   TTreeReaderArray<int> * m_HLT_Mu30_eta2p1_PFJet150_PFJet50_prescale;          // HLT_Mu30_eta2p1_PFJet150_PFJet50_prescale
//   TTreeReaderArray<int> * m_HLT_Mu40_eta2p1_PFJet200_PFJet50;          // HLT_Mu40_eta2p1_PFJet200_PFJet50
//   TTreeReaderArray<int> * m_HLT_Mu40_eta2p1_PFJet200_PFJet50_prescale;          // HLT_Mu40_eta2p1_PFJet200_PFJet50_prescale
//   TTreeReaderArray<int> * m_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50;          // HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50
//   TTreeReaderArray<int> * m_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale;          // HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale
//   TTreeReaderArray<int> * m_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;          // HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50
//   TTreeReaderArray<int> * m_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale;          // HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale
};

#endif
