/*
Created:        --
Last Updated:   22 August 2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Event class
Contains all the objects (& structs) with event information
-- check header for functions that return member variables
*/
#include "Analysis/CyMiniAna/interface/Event.h"

// constructor
Event::Event( TTreeReader &myReader, configuration &cmaConfig ) :
  m_config(&cmaConfig),
  m_ttree(myReader),
  m_treeName("SetMe"),
  m_fileName("SetMe"){
    m_isMC     = m_config->isMC();
    m_useTruth = m_config->useTruth();
    m_grid     = m_config->isGridFile();             // file directly from EDM->FlatNtuple step
    m_treeName = m_ttree.GetTree()->GetName();       // for systematics
    m_fileName = m_config->filename();               // for accessing file metadata

    // Neural network and kinematic reconstruction information
    m_isZeroLeptonAnalysis = m_config->isZeroLeptonAnalysis();
    m_isOneLeptonAnalysis  = m_config->isOneLeptonAnalysis();
    m_isTwoLeptonAnalysis  = m_config->isTwoLeptonAnalysis();
    m_mapOfContainment     = m_config->mapOfPartonContainment(); // containment map for truth tops matched to jets

    m_useJets       = m_config->useJets();
    m_useLargeRJets = m_config->useLargeRJets();
    m_useLeptons    = m_config->useLeptons();
    m_useNeutrinos  = m_config->useNeutrinos();
    m_neutrinoReco  = m_config->neutrinoReco();            // reconstruct neutrino
    m_DNNinference  = m_config->DNNinference();            // use DNN to predict values
    m_DNNtraining   = m_config->DNNtraining();             // load DNN features (save/use later)
    m_getDNN        = (m_DNNinference || m_DNNtraining);   // CWoLa
    m_useDNN        = m_config->useDNN();                  // Access CWoLa from TTree
    m_kinematicReco = m_config->kinematicReco();           // build the ttbar system

    // b-tagging working points
    m_CSVv2L = m_config->CSVv2L();
    m_CSVv2M = m_config->CSVv2M();
    m_CSVv2T = m_config->CSVv2T();


    //** Access branches from Tree **//
    m_eventNumber  = new TTreeReaderValue<unsigned long long>(m_ttree,"eventNumber");
    m_runNumber    = new TTreeReaderValue<unsigned int>(m_ttree,"runNumber");
    m_lumiblock    = new TTreeReaderValue<unsigned int>(m_ttree,"lumiblock");

    m_npv = new TTreeReaderValue<unsigned int>(m_ttree,"npv");
    m_rho = new TTreeReaderValue<float>(m_ttree,"rho");
    m_true_pileup = new TTreeReaderValue<unsigned int>(m_ttree,"true_pileup");

    /** Triggers **/
    m_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50 = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50");
    m_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165");
    m_HLT_Ele115_CaloIdVT_GsfTrkIdT    = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_Ele115_CaloIdVT_GsfTrkIdT");
    m_HLT_Mu40_Eta2P1_PFJet200_PFJet50 = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_Mu40_Eta2P1_PFJet200_PFJet50");
    m_HLT_Mu50    = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_Mu50");
    m_HLT_TkMu50  = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_TkMu50");
    m_HLT_PFHT800 = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_PFHT800");
    m_HLT_PFHT900 = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_PFHT900");
    m_HLT_AK8PFJet450 = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_AK8PFJet450");
    m_HLT_PFHT700TrimMass50  = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_PFHT700TrimMass50");
    m_HLT_PFJet360TrimMass30 = new TTreeReaderValue<unsigned int>(m_ttree,"HLT_PFJet360TrimMass30");
    //m_HLT_Ele45_WPLoose_Gsf = new TTreeReaderValue<int>(m_ttree,"HLT_Ele45_WPLoose_Gsf");
    //m_HLT_TkMu50 = new TTreeReaderValue<int>(m_ttree,"HLT_TkMu50");

    /** Filters **/
    m_Flag_goodVertices  = new TTreeReaderValue<unsigned int>(m_ttree,"Flag_goodVertices");
    m_Flag_eeBadScFilter = new TTreeReaderValue<unsigned int>(m_ttree,"Flag_eeBadScFilter");
    m_Flag_HBHENoiseFilter    = new TTreeReaderValue<unsigned int>(m_ttree,"Flag_HBHENoiseFilter");
    m_Flag_HBHENoiseIsoFilter = new TTreeReaderValue<unsigned int>(m_ttree,"Flag_HBHENoiseIsoFilter");
    m_Flag_globalTightHalo2016Filter = new TTreeReaderValue<unsigned int>(m_ttree,"Flag_globalTightHalo2016Filter");
    m_Flag_EcalDeadCellTriggerPrimitiveFilter = new TTreeReaderValue<unsigned int>(m_ttree,"Flag_EcalDeadCellTriggerPrimitiveFilter");

    /** JETS **/
    if (m_useJets){
      // small-R jet information
      m_jet_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4pt");
      m_jet_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4eta");
      m_jet_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4phi");
      m_jet_m   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4mass");
      m_jet_bdisc    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4bDisc");
      m_jet_deepCSV  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4deepCSV");
      m_jet_area     = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4area");
      m_jet_uncorrPt = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4uncorrPt");
      m_jet_uncorrE  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4uncorrE");
      m_jet_jerSF    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4jerSF");
      m_jet_jerSF_UP = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4jerSF_UP");
      m_jet_jerSF_DOWN = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4jerSF_DOWN");
    }

    if (m_useLargeRJets){
      // large-R Jet information
      m_ljet_pt     = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8pt");
      m_ljet_eta    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8eta");
      m_ljet_phi    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8phi");
      m_ljet_m      = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8mass");
      m_ljet_SDmass = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8SDmass");
      m_ljet_tau1   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8tau1");
      m_ljet_tau2   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8tau2");
      m_ljet_tau3   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8tau3");
      m_ljet_area   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8area");
      m_ljet_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8charge");
      m_ljet_chargeSD = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8chargeSD");
      m_ljet_charge3  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8charge3");
      m_ljet_charge10 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8charge10");

      m_ljet_subjet0_bdisc    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0bDisc");
      m_ljet_subjet0_deepCSV  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0deepCSV");
      m_ljet_subjet0_charge   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0charge");
      m_ljet_subjet0_charge3  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0charge3");
      m_ljet_subjet0_charge10 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0charge10");
      m_ljet_subjet0_pt   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0pt");
      m_ljet_subjet0_mass = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0mass");
      m_ljet_subjet0_tau1 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0tau1");
      m_ljet_subjet0_tau2 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0tau2");
      m_ljet_subjet0_tau3 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet0tau3");

      m_ljet_subjet1_bdisc    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1bDisc");
      m_ljet_subjet1_deepCSV  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1deepCSV");
      m_ljet_subjet1_charge   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1charge");
      m_ljet_subjet1_charge3  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1charge3");
      m_ljet_subjet1_charge10 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1charge10");
      m_ljet_subjet1_pt   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1pt");
      m_ljet_subjet1_mass = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1mass");
      m_ljet_subjet1_tau1 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1tau1");
      m_ljet_subjet1_tau2 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1tau2");
      m_ljet_subjet1_tau3 = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8subjet1tau3");

      m_ljet_BEST_class = new TTreeReaderValue<std::vector<int>>(m_ttree,"AK8BEST_class");
      m_ljet_BEST_t = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8BEST_t");
      m_ljet_BEST_w = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8BEST_w");
      m_ljet_BEST_z = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8BEST_z");
      m_ljet_BEST_h = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8BEST_h");
      m_ljet_BEST_j = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8BEST_j");

      m_ljet_uncorrPt = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8uncorrPt");
      m_ljet_uncorrE  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8uncorrE");
      m_ljet_jerSF    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8jerSF");
      m_ljet_jerSF_UP = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8jerSF_UP");
      m_ljet_jerSF_DOWN = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8jerSF_DOWN");
    }


    /** LEPTONS **/
    if (m_useLeptons){
      m_el_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELpt");
      m_el_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELeta");
      m_el_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELphi");
      m_el_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELenergy");
      m_el_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELcharge");
      //m_el_iso = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELiso");
      m_el_id_loose  = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"ELlooseID");
      m_el_id_medium = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"ELmediumID");
      m_el_id_tight  = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"ELtightID");
      m_el_id_loose_noIso  = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"ELlooseIDnoIso");
      m_el_id_medium_noIso = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"ELmediumIDnoIso");
      m_el_id_tight_noIso  = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"ELtightIDnoIso");

      m_mu_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUpt");
      m_mu_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUeta");
      m_mu_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUphi");
      m_mu_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUenergy");
      m_mu_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUcharge");
      m_mu_iso = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUcorrIso");
      m_mu_id_loose  = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"MUlooseID");
      m_mu_id_medium = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"MUmediumID");
      m_mu_id_tight  = new TTreeReaderValue<std::vector<unsigned int>>(m_ttree,"MUtightID");
    }

    m_met_met  = new TTreeReaderValue<float>(m_ttree,"METpt");
    m_met_phi  = new TTreeReaderValue<float>(m_ttree,"METphi");

    m_HTAK8    = new TTreeReaderValue<float>(m_ttree,"HTak8");
    m_HTAK4    = new TTreeReaderValue<float>(m_ttree,"HTak4");

    // set some event weights and access necessary branches
    m_xsection       = 1.0;
    m_kfactor        = 1.0;
    m_sumOfWeights   = 1.0;
    m_LUMI           = m_config->LUMI();

    Sample ss = m_config->sample();

    // MC information
    if (m_isMC){
      //m_weight_mc    = 1;//new TTreeReaderValue<float>(m_ttree,"evt_Gen_Weight");
      m_xsection     = ss.XSection;
      m_kfactor      = ss.KFactor;        // most likely =1
      m_sumOfWeights = ss.sumOfWeights;

      if (m_config->isTtbar()){
        m_mc_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"GENpt");
        m_mc_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"GENeta");
        m_mc_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"GENphi");
        m_mc_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"GENenergy");
        m_mc_pdgId  = new TTreeReaderValue<std::vector<int>>(m_ttree,"GENid");
        m_mc_status = new TTreeReaderValue<std::vector<int>>(m_ttree,"GENstatus");
        m_mc_parent_idx = new TTreeReaderValue<std::vector<int>>(m_ttree,"GENparent_idx");
        m_mc_child0_idx = new TTreeReaderValue<std::vector<int>>(m_ttree,"GENchild0_idx");
        m_mc_child1_idx = new TTreeReaderValue<std::vector<int>>(m_ttree,"GENchild1_idx");
        m_mc_isHadTop = new TTreeReaderValue<std::vector<int>>(m_ttree,"GENisHadTop");
      }
/*
      m_mc_ht = new TTreeReaderValue<float>(m_ttree,"evt_Gen_Ht");
      m_truth_jet_pt  = new TTreeReaderValue<float>(m_ttree,"jetAK4CHS_GenJetPt");
      m_truth_jet_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK4CHS_GenJetEta");
      m_truth_jet_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK4CHS_GenJetPhi");
      m_truth_jet_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK4CHS_GenJetCharge");
      m_truth_ljet_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK8CHS_GenJetPt");
      m_truth_ljet_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK8CHS_GenJetEta");
      m_truth_ljet_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK8CHS_GenJetPhi");
      m_truth_ljet_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK8CHS_GenJetE");
      m_truth_ljet_charge     = new TTreeReaderValue<std::vector<float>>(m_ttree,"jetAK8CHS_GenJetCharge");
      m_truth_ljet_subjet_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"subjetAK8CHS_GenJetCharge");
*/
    } // end isMC


    // Truth matching tool
    m_truthMatchingTool = new truthMatching(cmaConfig);
    m_truthMatchingTool->initialize();

    // DNN material
    m_deepLearningTool = new deepLearning(cmaConfig);
    if (!m_getDNN && m_useDNN)
        m_dnn_score = new TTreeReaderValue<float>(m_ttree,"ljet_CWoLa");


    // Kinematic reconstruction algorithms
    m_ttbarRecoTool    = new ttbarReco(cmaConfig);
    m_neutrinoRecoTool = new neutrinoReco(cmaConfig);

    if (!m_neutrinoReco && m_useNeutrinos){
        // Neutrinos aren't stored in the baseline ntuples, requires 'kinematicReco' to create
        m_nu_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree, "nu_pt");
        m_nu_eta = new TTreeReaderValue<std::vector<float>>(m_ttree, "nu_eta");
        m_nu_phi = new TTreeReaderValue<std::vector<float>>(m_ttree, "nu_phi");
    }

    if (!m_kinematicReco){
        // Load ttbar variables from file
        m_leptop_jet  = new TTreeReaderValue<int>(m_ttree, "leptop_jet");
        m_hadtop_ljet = new TTreeReaderValue<int>(m_ttree, "hadtop_ljet");
    }
} // end constructor


Event::~Event() {}



void Event::initialize_eventWeights(){
    /* Create vectors of the systematics that are weights for the nominal events
       Must be called from the constructor for the access to TTreeReaderValues to work!
    */
    std::map<std::string,unsigned int> mapWeightSystematics = m_config->mapOfWeightVectorSystematics();

    m_listOfWeightSystematics = m_config->listOfWeightSystematics();

    m_weightSystematicsFloats.clear();
    m_weightSystematicsVectorFloats.clear();

    // systematics from the nominal tree that are floats
    for (const auto& nom_syst : m_listOfWeightSystematics){
        if (!m_useLeptons && nom_syst.find("leptonSF")!=std::string::npos)
            continue;
        m_weightSystematicsFloats[nom_syst] = new TTreeReaderValue<float>(m_ttree,nom_syst.c_str());
    }

    // systematics from the nominal tree that are vectors
    for (const auto& syst : mapWeightSystematics)
        m_weightSystematicsVectorFloats[syst.first] = new TTreeReaderValue<float>(m_ttree,syst.first.c_str());

    return;
}



void Event::updateEntry(Long64_t entry){
    /* Update the entry -> update all TTree variables */
    cma::DEBUG("EVENT : Update Entry "+std::to_string(entry) );
    m_entry = entry;

    // make sure the entry exists
    if(isValidRecoEntry())
        m_ttree.SetEntry(m_entry);
    else
        cma::ERROR("EVENT : Invalid Reco entry "+std::to_string(m_entry)+"!");

    return;
}


void Event::clear(){
    /* Clear many of the vectors/maps for each event -- SAFETY PRECAUTION */
    m_truth_ljets.clear();
    m_truth_jets.clear();
    m_truth_leptons.clear();
    m_truth_neutrinos.clear();

    m_jets.clear();
    m_ljets.clear();
    m_leptons.clear();
    m_neutrinos.clear();

    m_btag_jets.clear();
    m_btag_jets_default.clear();
    m_weight_btag_default = 1.0;
    m_nominal_weight = 1.0;

    m_HT = 0;
    m_ST = 0;

    return;
}


void Event::execute(Long64_t entry){
    /* Get the values from the event */
    cma::DEBUG("EVENT : Execute event " );

    // Load data from root tree for this event
    updateEntry(entry);

    // Reset many event-level values
    clear();

    // Get the event weights (for cutflow & histograms)
    initialize_weights();
    cma::DEBUG("EVENT : Setup weights ");

    // Filters
    initialize_filters();

    // Triggers
    initialize_triggers();

    // Truth Information
    if (m_useTruth){
        initialize_truth();
        cma::DEBUG("EVENT : Setup truth information ");
    }

    // Jets
    if (m_useJets){
        initialize_jets();
        cma::DEBUG("EVENT : Setup small-R jets ");
    }

    // Leptons
    if (m_useLeptons){
        initialize_leptons();
        cma::DEBUG("EVENT : Setup leptons ");
    }

    // Large-R Jets
    if (m_useLargeRJets){
        initialize_ljets();
        cma::DEBUG("EVENT : Setup large-R jets ");
    }

    // Get some kinematic variables (MET, HT, ST)
    initialize_kinematics();
    cma::DEBUG("EVENT : Setup kinematic variables ");


    // Neutrinos
    if (m_useNeutrinos){
        // relies on kinematic reconstruction, unless the information is saved in root file
        initialize_neutrinos();
        cma::DEBUG("EVENT : Setup neutrinos ");
    }

    // Kinematic reconstruction (if they values aren't in the root file)
    m_ttbar0L = {};
    m_ttbar1L = {};
    m_ttbar2L = {};
    if (m_kinematicReco) ttbarReconstruction();
    else{
        Jet dummy_jet;
        Ljet dummy_ljet;
        Neutrino dummy_nu;
        Lepton dummy_lep;

        if (m_isOneLeptonAnalysis){
            int ljetidx = **m_hadtop_ljet;
            int jetidx  = **m_leptop_jet;
            m_ttbar1L.jet  = (jetidx>=0)  ? m_jets.at(jetidx)   : dummy_jet;
            m_ttbar1L.ljet = (ljetidx>=0) ? m_ljets.at(ljetidx) : dummy_ljet;

            if (m_muons.size()>0)
                m_ttbar1L.lepton = m_muons.at(0);
            else if (m_electrons.size()>0)
                m_ttbar1L.lepton = m_electrons.at(0);
            else if (m_muons.size()+m_electrons.size() == 0)
                m_ttbar1L.lepton = dummy_lep;

            m_ttbar1L.neutrino = (m_neutrinos.size()>0) ? m_neutrinos.at(0) : dummy_nu;
        }
    }

    cma::DEBUG("EVENT : Setup Event ");

    return;
}


void Event::ttbarReconstruction(){
    /* Reconstruct ttbar system -- after event selection! */
    m_ttbar0L = {};
    m_ttbar1L = {};
    m_ttbar2L = {};

    if (m_isZeroLeptonAnalysis){
        m_ttbarRecoTool->execute(m_ljets);
        m_ttbar0L = m_ttbarRecoTool->ttbar0L();
    }
    if (m_isOneLeptonAnalysis){
        m_ttbarRecoTool->execute(m_leptons,m_neutrinos,m_jets,m_ljets);
        m_ttbar1L = m_ttbarRecoTool->ttbar1L();
    }
    if (m_isTwoLeptonAnalysis){
        m_ttbarRecoTool->execute(m_electrons,m_muons,m_jets);
        m_ttbar2L = m_ttbarRecoTool->ttbar2L();
    }

    return;
}


void Event::initialize_filters(){
    /* Setup the filters */
    m_filters.clear();

    m_filters["goodVertices"] = **m_Flag_goodVertices;
    m_filters["eeBadScFilter"] = **m_Flag_eeBadScFilter;
    m_filters["HBHENoiseFilter"] = **m_Flag_HBHENoiseFilter;
    m_filters["HBHENoiseIsoFilter"] = **m_Flag_HBHENoiseIsoFilter;
    m_filters["globalTightHalo2016Filter"] = **m_Flag_globalTightHalo2016Filter;
    m_filters["EcalDeadCellTriggerPrimitiveFilter"] = **m_Flag_EcalDeadCellTriggerPrimitiveFilter;

    return;
}


void Event::initialize_triggers(){
    /* Setup triggers */
    m_triggers.clear();

    m_triggers["HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50"] = **m_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
    m_triggers["HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"] = **m_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
    m_triggers["HLT_Ele115_CaloIdVT_GsfTrkIdT"]    = **m_HLT_Ele115_CaloIdVT_GsfTrkIdT;
    m_triggers["HLT_Mu40_Eta2P1_PFJet200_PFJet50"] = **m_HLT_Mu40_Eta2P1_PFJet200_PFJet50;
    m_triggers["HLT_Mu50"]    = **m_HLT_Mu50;
    m_triggers["HLT_TkMu50"]  = **m_HLT_TkMu50;
    m_triggers["HLT_PFHT800"] = **m_HLT_PFHT800;
    m_triggers["HLT_PFHT900"] = **m_HLT_PFHT900;
    m_triggers["HLT_AK8PFJet450"] = **m_HLT_AK8PFJet450;
    m_triggers["HLT_PFHT700TrimMass50"]  = **m_HLT_PFHT700TrimMass50;
    m_triggers["HLT_PFJet360TrimMass30"] = **m_HLT_PFJet360TrimMass30;

    return;
}


void Event::initialize_truth(){
    /* Setup struct of truth information */
    m_truth_partons.clear();
    m_truth_tops.clear();

    if (!m_config->isTtbar()) return;   // don't need this for MC other than ttbar

    // only care about this for ttbar
    unsigned int nPartons( (*m_mc_pt)->size() );
    cma::DEBUG("EVENT : N Partons = "+std::to_string(nPartons));

    // Collect truth top information into one value
    unsigned int t_idx(0);  // keeping track of tops in m_truth_tops

    // loop over truth partons
    unsigned int p_idx(0);
    for (unsigned int i=0; i<nPartons; i++){

        Parton parton;
        parton.p4.SetPtEtaPhiE((*m_mc_pt)->at(i),(*m_mc_eta)->at(i),(*m_mc_phi)->at(i),(*m_mc_e)->at(i));

        int status = (*m_mc_status)->at(i);
        int pdgId  = (*m_mc_pdgId)->at(i);
        unsigned int abs_pdgId = std::abs(pdgId);

        parton.pdgId  = pdgId;
        parton.status = status;

        // simple booleans for type
        parton.isTop = ( abs_pdgId==6 );
        parton.isW   = ( abs_pdgId==24 );
        parton.isLepton = ( abs_pdgId>=11 && abs_pdgId<=16 );
        parton.isQuark  = ( abs_pdgId<7 );

        if (parton.isLepton){
            parton.isTau  = ( abs_pdgId==15 ) ? 1 : 0;
            parton.isMuon = ( abs_pdgId==13 ) ? 1 : 0;
            parton.isElectron = ( abs_pdgId==11 ) ? 1 : 0;
            parton.isNeutrino = ( abs_pdgId==12 || abs_pdgId==14 || abs_pdgId==16 ) ? 1 : 0;
        }
        else if (parton.isQuark){
            parton.isLight  = ( abs_pdgId<5 ) ? 1 : 0;
            parton.isBottom = ( abs_pdgId==5 ) ? 1 : 0;
        }

        parton.index      = p_idx;                    // index in vector of truth_partons
        parton.top_index  = -1;                       // index in truth_tops vector
        parton.containment = 0;                       // value for containment calculation

        parton.parent_idx = (*m_mc_parent_idx)->at(i);
        parton.child0_idx = (*m_mc_child0_idx)->at(i);
        parton.child1_idx = (*m_mc_child1_idx)->at(i);

        // skip replicated top/W in truth record
        if (parton.isTop && parton.status<60) continue;
        if (parton.isW && (parton.child0_idx<0 || parton.child1_idx<0)) continue;

        // build truth top structs
        // in truth parton record, the top should arrive before its children
        TruthTop top;

        if (parton.isTop){
            cma::DEBUG("EVENT : is top ");
            top.Wdecays.clear();    // for storing W daughters
            top.daughters.clear();  // for storing non-W/bottom daughters

            top.Top       = parton.index;
            top.isTop     = (pdgId>0);
            top.isAntiTop = (pdgId<0);
            top.isHadronic = (*m_mc_isHadTop)->at(p_idx);
            top.isLeptonic = !(*m_mc_isHadTop)->at(p_idx);
            parton.top_index   = t_idx;
            parton.containment = m_mapOfContainment.at("FULL");   // only considering truth tops right now, not the decay products
            if (parton.pdgId<0) parton.containment *= -1;         // negative value for anti-tops

            m_truth_tops.push_back(top);   // store tops now, add information from children in future iterations
            t_idx++;
        }
        else if (!parton.isTop && parton.parent_idx>0) {
            int parent_pdgid = (*m_mc_pdgId)->at(parton.parent_idx);
            cma::DEBUG("EVENT : it's not a top, it's a "+std::to_string(pdgId)+"; parent idx = "+std::to_string(parton.parent_idx)+"; parent pdgid = "+std::to_string(parent_pdgid));

            // check if W is decaying to itself
            if (std::abs(parent_pdgid) == 24 && parent_pdgid == parton.pdgId) {// look at grandparent
                int gparent_idx = (*m_mc_parent_idx)->at(parton.parent_idx);
                parent_pdgid = (*m_mc_pdgId)->at(gparent_idx);
            }
            else if (parent_pdgid==parton.pdgId) continue;    // other particles self-decaying, just skip

            // get the parent from the list of partons
            Parton parent;
            int top_index(-1);
            for (const auto& t : m_truth_partons){
                if (t.pdgId==parent_pdgid) {
                    parent    = t;
                    top_index = t.top_index;
                    break;
                }
            }
            if (top_index<0) continue;    // weird element in truth record, just skip it
            parton.top_index = top_index;
            cma::DEBUG("EVENT : Top index = "+std::to_string(top_index));

            // Parent is Top (W or b)
            if (parent.isTop){
                top = m_truth_tops.at(parent.top_index);
                if (parton.isW) top.W = parton.index;
                else if (parton.isBottom) {
                    top.bottom = parton.index;
                    parton.containment = m_mapOfContainment.at("BONLY");
                    if (top.isAntiTop) parton.containment*=-1;
                }
                else top.daughters.push_back( parton.index );        // non-W/bottom daughter
                m_truth_tops[parent.top_index] = top;                // update entry
            }
            // Parent is W
            else if (parent.isW){
                top = m_truth_tops.at(top_index);
                top.Wdecays.push_back(parton.index);
                top.isHadronic = (parton.isQuark);
                top.isLeptonic = (parton.isLepton);

                parton.containment = m_mapOfContainment.at("QONLY");
                if (top.isAntiTop) parton.containment*=-1;

                m_truth_tops[top_index] = top;      // update entry
            }
        } // end else if not top

        // store for later access
        m_truth_partons.push_back( parton );
        p_idx++;
    } // end loop over truth partons

    m_truthMatchingTool->setTruthPartons(m_truth_partons);
    m_truthMatchingTool->setTruthTops(m_truth_tops);

    return;
}


void Event::initialize_jets(){
    /* Setup struct of jets (small-r) and relevant information 
     * b-tagging: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
        CSVv2L 0.5426
        CSVv2M 0.8484
        CSVv2T 0.9535
     */
    unsigned int nJets = (*m_jet_pt)->size();
    m_jets.clear();
    m_jets_iso.clear();  // jet collection for lepton 2D isolation

    for (const auto& btagWP : m_config->btagWkpts() ){
        m_btag_jets[btagWP].clear();
    }

    unsigned int idx(0);
    unsigned int idx_iso(0);
    for (unsigned int i=0; i<nJets; i++){
        Jet jet;
        jet.p4.SetPtEtaPhiM( (*m_jet_pt)->at(i),(*m_jet_eta)->at(i),(*m_jet_phi)->at(i),(*m_jet_m)->at(i));

        bool isGoodIso( jet.p4.Pt()>15 && std::abs(jet.p4.Eta())<2.4);
        bool isGood(jet.p4.Pt()>50 && std::abs(jet.p4.Eta())<2.4);

        if (!isGood && !isGoodIso) continue;

        jet.isGood = isGood;

        jet.bdisc    = (*m_jet_bdisc)->at(i);
        jet.deepCSV  = (*m_jet_deepCSV)->at(i);
        jet.area     = (*m_jet_area)->at(i);
        jet.uncorrE  = (*m_jet_uncorrE)->at(i);
        jet.uncorrPt = (*m_jet_uncorrPt)->at(i);
        jet.jerSF    = (*m_jet_jerSF)->at(i);
        jet.jerSF_UP = (*m_jet_jerSF_UP)->at(i);
        jet.jerSF_DOWN = (*m_jet_jerSF_DOWN)->at(i);

        jet.index  = idx;

        if (isGood){
            m_jets.push_back(jet);
            getBtaggedJets(jet);          // only care about b-tagging for 'real' AK4
            idx++;
        }
        if (isGoodIso){
            m_jets_iso.push_back(jet);    // used for 2D isolation
            idx_iso++;
        }
    }

    m_btag_jets_default = m_btag_jets.at(m_config->jet_btagWkpt());

    return;
}


void Event::initialize_ljets(){
    /* Setup struct of large-R jets and relevant information 
      0 :: Top      (lepton Q < 0)
      1 :: Anti-top (lepton Q > 0)
    */
    unsigned int nLjets = (*m_ljet_pt)->size();
    m_ljets.clear();

    // Define CWoLa classification based on lepton charge (only single lepton events)
    int target(-1);
    if (m_config->isOneLeptonAnalysis() && m_leptons.size()>0){
        int charge = m_leptons.at(0).charge;
        target = (charge>0) ? 1:0;
    }

    unsigned int idx(0);
    for (unsigned int i=0; i<nLjets; i++){
        Ljet ljet;
        ljet.p4.SetPtEtaPhiM( (*m_ljet_pt)->at(i),(*m_ljet_eta)->at(i),(*m_ljet_phi)->at(i),(*m_ljet_m)->at(i));
        ljet.softDropMass = (*m_ljet_SDmass)->at(i);

        ljet.tau1   = (*m_ljet_tau1)->at(i);
        ljet.tau2   = (*m_ljet_tau2)->at(i);
        ljet.tau3   = (*m_ljet_tau3)->at(i);
        ljet.tau21  = ljet.tau2 / ljet.tau1;
        ljet.tau32  = ljet.tau3 / ljet.tau2;
        //bool toptag = (ljet.softDropMass>105. && ljet.softDropMass<210 && ljet.tau32<0.65);  // apply in eventSelection

        float subjet0_bdisc = (*m_ljet_subjet0_bdisc)->at(i);  // want the subjets to have "real" b-disc values
        float subjet1_bdisc = (*m_ljet_subjet1_bdisc)->at(i);

        // check if the AK8 is 'good'
        bool isGood(ljet.p4.Pt()>400. && fabs(ljet.p4.Eta())<2.4 && subjet0_bdisc>=0 && subjet1_bdisc>=0);
        ljet.isGood = isGood;

        if (!isGood) continue;

        ljet.BEST_t = (*m_ljet_BEST_t)->at(i);
        ljet.BEST_w = (*m_ljet_BEST_w)->at(i);
        ljet.BEST_z = (*m_ljet_BEST_z)->at(i);
        ljet.BEST_h = (*m_ljet_BEST_h)->at(i);
        ljet.BEST_j = (*m_ljet_BEST_j)->at(i);
        ljet.BEST_class = (*m_ljet_BEST_class)->at(i);

        ljet.subjet0_bdisc  = subjet0_bdisc;   // (*m_ljet_subjet0_bdisc)->at(i);
        ljet.subjet0_charge = (*m_ljet_subjet0_charge)->at(i);
        ljet.subjet0_mass   = (*m_ljet_subjet0_mass)->at(i);
        ljet.subjet0_pt     = (*m_ljet_subjet0_pt)->at(i);
        ljet.subjet0_tau1   = (*m_ljet_subjet0_tau1)->at(i);
        ljet.subjet0_tau2   = (*m_ljet_subjet0_tau2)->at(i);
        ljet.subjet0_tau3   = (*m_ljet_subjet0_tau3)->at(i);

        ljet.subjet1_bdisc  = subjet1_bdisc;   // (*m_ljet_subjet1_bdisc)->at(i);
        ljet.subjet1_charge = (*m_ljet_subjet1_charge)->at(i);
        ljet.subjet1_mass   = (*m_ljet_subjet1_mass)->at(i);
        ljet.subjet1_pt     = (*m_ljet_subjet1_pt)->at(i);
        ljet.subjet1_tau1   = (*m_ljet_subjet1_tau1)->at(i);
        ljet.subjet1_tau2   = (*m_ljet_subjet1_tau2)->at(i);
        ljet.subjet1_tau3   = (*m_ljet_subjet1_tau3)->at(i);

        ljet.charge = (*m_ljet_charge)->at(i);
        ljet.target = target;
        ljet.index  = idx;

        ljet.area     = (*m_ljet_area)->at(i);
        ljet.uncorrE  = (*m_ljet_uncorrE)->at(i);
        ljet.uncorrPt = (*m_ljet_uncorrPt)->at(i);

        ljet.jerSF    = 1.0; //(*m_ljet_jerSF)->at(i);
        ljet.jerSF_UP = 1.0; //(*m_ljet_jerSF_UP)->at(i);
        ljet.jerSF_DOWN = 1.0; //(*m_ljet_jerSF_DOWN)->at(i);

        // Truth-matching to jet
        ljet.truth_partons.clear();
        if (m_useTruth && m_config->isTtbar()) {
            cma::DEBUG("EVENT : Truth match AK8");          // match subjets (and then the AK8 jet) to truth tops

            m_truthMatchingTool->matchJetToTruthTop(ljet);  // match to partons

            cma::DEBUG("EVENT : ++ Ljet had top = "+std::to_string(ljet.isHadTop)+" for truth top "+std::to_string(ljet.matchId));
        } // end truth matching ljet to partons

        m_ljets.push_back(ljet);
        idx++;
    }

    deepLearningPrediction();   // store features in map (easily access later)

    return;
}


void Event::initialize_leptons(){
    /* Setup struct of lepton and relevant information */
    m_leptons.clear();
    m_electrons.clear();  // not using right now
    m_muons.clear();      // not using right now

    // Muons
    unsigned int nMuons = (*m_mu_pt)->size();

    for (unsigned int i=0; i<nMuons; i++){
        Lepton mu;
        mu.p4.SetPtEtaPhiE( (*m_mu_pt)->at(i),(*m_mu_eta)->at(i),(*m_mu_phi)->at(i),(*m_mu_e)->at(i));
        bool isMedium   = (*m_mu_id_medium)->at(i);
        bool isTight    = (*m_mu_id_tight)->at(i);

        bool iso = customIsolation(mu);    // 2D isolation cut between leptons & AK4 (need AK4 initialized first!)

        bool isGood(mu.p4.Pt()>50 && std::abs(mu.p4.Eta())<2.4 && isMedium && iso);
        mu.isGood = isGood;

        if (!isGood) continue;

        mu.charge = (*m_mu_charge)->at(i);
        mu.loose  = (*m_mu_id_loose)->at(i);
        mu.medium = isMedium; 
        mu.tight  = isTight; 
        mu.iso    = iso;       // use 2D isolation instead -- (*m_mu_iso)->at(i);

        mu.isMuon = true;
        mu.isElectron = false;

        m_leptons.push_back(mu);
    }

    // Electrons
    unsigned int nElectrons = (*m_el_pt)->size();
    for (unsigned int i=0; i<nElectrons; i++){
        Lepton el;
        el.p4.SetPtEtaPhiE( (*m_el_pt)->at(i),(*m_el_eta)->at(i),(*m_el_phi)->at(i),(*m_el_e)->at(i));
        bool isTightNoIso  = (*m_el_id_tight_noIso)->at(i);
        bool isMediumNoIso = (*m_el_id_medium_noIso)->at(i);

        bool iso = customIsolation(el);    // 2D isolation cut between leptons & AK4 (need AK4 initialized first!)

        bool isGood(el.p4.Pt()>50 && std::abs(el.p4.Eta())<2.4 && isMediumNoIso && iso);
        el.isGood = isGood;

        if (!isGood) continue;

        el.charge = (*m_el_charge)->at(i);
        el.loose  = (*m_el_id_loose)->at(i);
        el.medium = (*m_el_id_medium)->at(i);
        el.tight  = (*m_el_id_tight)->at(i);
        el.loose_noIso  = (*m_el_id_loose_noIso)->at(i);
        el.medium_noIso = isMediumNoIso;
        el.tight_noIso  = isTightNoIso;
        el.iso = iso;      // 2D isolation with ID+noIso

        el.isMuon = false;
        el.isElectron = true;

        m_leptons.push_back(el);
    }

    cma::DEBUG("EVENT : Found "+std::to_string(m_leptons.size())+" leptons!");

    return;
}


void Event::initialize_neutrinos(){
    /* Build the neutrinos */
    m_neutrinos.clear();

    if (!m_isOneLeptonAnalysis && !m_isTwoLeptonAnalysis){
        cma::ERROR("EVENT : 'initialize_neutrinos()' called for analysis without neutrinos!");
        return;
    }

    Neutrino nu1;
    nu1.p4.SetPtEtaPhiM( m_met.p4.Pt(), 0, m_met.p4.Phi(), 0);   // "dummy"; pz=0
    Neutrino nu2;
    nu2.p4.SetPtEtaPhiM( m_met.p4.Pt(), 0, m_met.p4.Phi(), 0);   // "dummy"; pz=0

    int nlep = m_leptons.size();
    if (nlep<1){
        // not enough leptons to do reconstruction, so just create dummy value
        m_neutrinos.push_back(nu1);
        return;
    }

    m_neutrinoRecoTool->setObjects(m_leptons.at(0),m_met);
    if (m_neutrinoReco){
        // reconstruct neutrinos!
        if (m_isOneLeptonAnalysis){
            nu1 = m_neutrinoRecoTool->execute();        // tool assumes 1-lepton final state
            m_neutrinos.push_back(nu1);
        }
        else if (m_isTwoLeptonAnalysis){
            // part of ttbar reconstruction instead?
            m_neutrinos.push_back(nu1);
            m_neutrinos.push_back(nu2);
        }
    }
    else{
        // Assign neutrinos from root file
        unsigned int nNus = (*m_nu_pt)->size();

        for (unsigned int i=0; i<nNus; i++){
            Neutrino nu;
            nu.p4.SetPtEtaPhiM( (*m_nu_pt)->at(i),(*m_nu_eta)->at(i),(*m_nu_phi)->at(i),0.0);
            m_neutrinos.push_back(nu);
        }
    } // end use neutrinos but don't reconstruct them

    return;
}


void Event::initialize_weights(){
    /* Event weights */
    m_nominal_weight = 1.0;

    m_weight_btag.clear();
    if (m_isMC){
        m_nominal_weight  = 1.0; //(**m_weight_pileup) * (**m_weight_mc);
        m_nominal_weight *= (m_xsection) * (m_kfactor) * m_LUMI / (m_sumOfWeights);
/*      // event weights
        m_weight_btag["70"] = (**m_weight_btag_70);
        m_weight_btag["77"] = (**m_weight_btag_77);
        m_weight_btag_default = m_weight_btag[m_config->jet_btagWkpt()];
        m_nominal_weight *= m_weight_btag_default;
*/
    }

    return;
}


void Event::initialize_kinematics(){
    /* Kinematic variables (HT, ST, MET) */
    m_HT_ak4 = **m_HTAK4;
    m_HT_ak8 = **m_HTAK8;

    m_HT = 0.0;
    m_ST = 0.0;

    // Get hadronic transverse energy
    if (m_useJets){
        // include small-R jet pT
        for (auto &small_jet : m_jets ){
            m_HT += small_jet.p4.Pt();
        }
    }
    else{
        // include large-R jet pT
        for (auto &large_jet : m_ljets){
            m_HT += large_jet.p4.Pt();
        }
    }

    // set MET
    m_met.p4.SetPtEtaPhiM(**m_met_met,0.,**m_met_phi,0.);
    cma::DEBUG("EVENT : MET = "+std::to_string(m_met.p4.Pt()));

    // Get MET and lepton transverse energy
    m_ST += m_HT;
    m_ST += m_met.p4.Pt();

    if (m_useLeptons){
        for (const auto& lep : m_leptons)
            m_ST += lep.p4.Pt(); 
    }

    // transverse mass of the W (only relevant for 1-lepton)
    float mtw(0.0);

    if (m_leptons.size()>0){
        Lepton lep = m_leptons.at(0);
        float dphi = m_met.p4.Phi() - lep.p4.Phi();
        mtw = sqrt( 2 * lep.p4.Pt() * m_met.p4.Pt() * (1-cos(dphi)) );
    }
    m_met.mtw = mtw;

    return;
}



bool Event::customIsolation( Lepton& lep ){
    /* 2D isolation cut for leptons 
       - Check that the lepton and nearest AK4 jet satisfies
         DeltaR() < 0.4 || pTrel>25
    */
    bool pass(false);
    //int min_index(-1);                    // index of AK4 closest to lep
    float drmin(100.0);                   // min distance between lep and AK4s
    float ptrel(0.0);                     // pTrel between lepton and AK4s

    if (m_jets_iso.size()<1) return false;    // no AK4 -- event will fail anyway

    for (const auto& jet : m_jets_iso){
        float dr = lep.p4.DeltaR( jet.p4 );
        if (dr < drmin) {
            drmin = dr;
            ptrel = cma::ptrel( lep.p4,jet.p4 );
            //min_index = jet.index;
        }
    }

    lep.drmin = drmin;
    lep.ptrel = ptrel;

    if (drmin > 0.4 || ptrel > 30) pass = true;

    return pass;
}


/*** GETTER FUNCTIONS ***/
void Event::getBtaggedJets( Jet& jet ){
    /* Determine the b-tagging */
    jet.isbtagged["L"] = false;
    jet.isbtagged["M"] = false;
    jet.isbtagged["T"] = false;

    if (jet.bdisc > m_CSVv2L){
        jet.isbtagged["L"] = true;
        m_btag_jets["L"].push_back(jet.index);  // 0 = index of this jet
        if (jet.bdisc > m_CSVv2M){
            jet.isbtagged["M"] = true;
            m_btag_jets["M"].push_back(jet.index);
            if (jet.bdisc > m_CSVv2T){
                jet.isbtagged["T"] = true;
                m_btag_jets["T"].push_back(jet.index);
            }
        }
    }

    return;
}


double Event::getSystEventWeight( const std::string &syst, const int weightIndex ){
    /* Calculate the event weight given some systematic
       -- only call for nominal events and systematic weights
       -- for non-nominal tree systematics, use the nominal event weight

       @param syst          Name of systematic (nominal or some weight systematic)
       @param weightIndex   Index of btagging SF; default to -1
    */
    double syst_event_weight(1.0);

    if (syst.compare("nominal")==0){
        // nominal event weight
        syst_event_weight  = m_nominal_weight;
    }
/*    else if (syst.find("leptonSF")!=std::string::npos){
        // leptonSF event weight
        syst_event_weight  = (**m_weight_pileup) * (**m_weight_mc);
        syst_event_weight *= m_weight_btag_default;
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);

        syst_event_weight *= **m_weightSystematicsFloats.at(syst);
    }
    else if (syst.find("bTagSF")!=std::string::npos){
        // bTagSF event weight -- check indices for eigenvector systematics
        syst_event_weight  = (**m_weight_pileup) * (**m_weight_mc);
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);
    }
    else{
        // safety to catch something weird -- just return 1.0
        cma::WARNING("EVENT : Passed systematic variation, "+syst+", to Event::getSystEventWeight() ");
        cma::WARNING("EVENT : that is inconsistent with the CyMiniAna options of ");
        cma::WARNING("EVENT :     nominal, jvt, pileup, leptonSF, and bTagSF. ");
        cma::WARNING("EVENT : Returning a weight of 1.0. ");
        syst_event_weight = 1.0;
    }
*/
    return syst_event_weight;
}




/*** RETURN PHYSICS INFORMATION ***/
std::vector<int> Event::btag_jets(const std::string &wkpt) const{
    /* Small-R Jet b-tagging */
    std::string tmp_wkpt(wkpt);
    if(m_btag_jets.find(wkpt) == m_btag_jets.end()){
        cma::WARNING("EVENT : B-tagging working point "+wkpt+" does not exist.");
        cma::WARNING("EVENT : Return vector of b-tagged jets for default working point "+m_config->jet_btagWkpt());
        tmp_wkpt = m_config->jet_btagWkpt();
    }
    return m_btag_jets.at(tmp_wkpt);
}

void Event::deepLearningPrediction(){
    /* Deep learning for large-R jets -- CWoLa */
    if (m_DNNinference){
        cma::DEBUG("EVENT : Calculate DNN ");
        m_deepLearningTool->inference(m_ljets);     // decorate the ljet with DNN values
    }
    else if (m_DNNtraining){
        for (auto& ljet : m_ljets){
            m_deepLearningTool->training(ljet);
            ljet.features = m_deepLearningTool->features();  // store the features on the ljet to make easily accessible later
        }
    }

    return;
}


/*** RETURN WEIGHTS ***/
float Event::weight_mc(){
    return 1.0; //**m_weight_mc;
}
float Event::weight_pileup(){
    return 1.0; //**m_weight_pileup;
}

float Event::weight_btag(const std::string &wkpt){
    std::string tmp_wkpt(wkpt);
    if(m_weight_btag.find(wkpt) == m_weight_btag.end()){
        cma::WARNING("EVENT : B-tagging working point "+wkpt+" does not exist");
        cma::WARNING("EVENT : Return calo-jet b-tag SF for default working point "+m_config->jet_btagWkpt());
        tmp_wkpt = m_config->jet_btagWkpt();
    }
    return m_weight_btag[tmp_wkpt];
}

float Event::weight_btag(){
    /* Default b-tag weight */
    return m_weight_btag_default;
}

// Get weight systematics
std::map<std::string,float> Event::weightSystematicsFloats(){
    /* systematics floats */
    std::map<std::string,float> tmp_weightSystematicsFloats;
    for (const auto& wsf : m_weightSystematicsFloats)
        tmp_weightSystematicsFloats[wsf.first] = **wsf.second;

    return tmp_weightSystematicsFloats;
}

std::map<std::string,std::vector<float> > Event::weightSystematicsVectorFloats(){
    /* weight systematics stored as vectors */
    std::map<std::string,std::vector<float> > tmp_weightSystematicsVectorFloats;
    return tmp_weightSystematicsVectorFloats;
}

std::vector<std::string> Event::listOfWeightSystematics(){
    /* list of weight systematics */
    return m_listOfWeightSystematics;
}



/*** DELETE VARIABLES ***/
void Event::finalize(){
    // delete variables
    cma::DEBUG("EVENT : Finalize() ");
    delete m_eventNumber;
    delete m_runNumber;
    delete m_lumiblock;

    if (m_useJets){
      delete m_jet_pt;
      delete m_jet_eta;
      delete m_jet_phi;
      delete m_jet_m;
      delete m_jet_bdisc;
      delete m_jet_deepCSV;
      delete m_jet_area;
      delete m_jet_uncorrPt;
      delete m_jet_uncorrE;
    }

    if (m_useLargeRJets){
      delete m_ljet_pt;
      delete m_ljet_eta;
      delete m_ljet_phi;
      delete m_ljet_m;
      delete m_ljet_tau1;
      delete m_ljet_tau2;
      delete m_ljet_tau3;
      delete m_ljet_charge;
      delete m_ljet_SDmass;
      delete m_ljet_subjet0_charge;
      delete m_ljet_subjet0_bdisc;
      delete m_ljet_subjet0_deepCSV;
      delete m_ljet_subjet0_pt;
      delete m_ljet_subjet0_mass;
      delete m_ljet_subjet1_charge;
      delete m_ljet_subjet1_bdisc;
      delete m_ljet_subjet1_deepCSV;
      delete m_ljet_subjet1_pt;
      delete m_ljet_subjet1_mass;
      delete m_ljet_area;
      delete m_ljet_uncorrPt;
      delete m_ljet_uncorrE;
    }

    if (m_useLeptons){
      delete m_el_pt;
      delete m_el_eta;
      delete m_el_phi;
      delete m_el_e;
      delete m_el_charge;
      //delete m_el_iso;
      delete m_el_id_loose;
      delete m_el_id_medium;
      delete m_el_id_tight;
      delete m_el_id_loose_noIso;
      delete m_el_id_medium_noIso;
      delete m_el_id_tight_noIso;

      delete m_mu_pt;
      delete m_mu_eta;
      delete m_mu_phi;
      delete m_mu_e;
      delete m_mu_charge;
      delete m_mu_iso;
      delete m_mu_id_loose;
      delete m_mu_id_medium;
      delete m_mu_id_tight;
    }

    delete m_met_met;
    delete m_met_phi;
    delete m_HTAK8;
    delete m_HTAK4;

    delete m_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
    delete m_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
    delete m_HLT_Ele115_CaloIdVT_GsfTrkIdT;
    delete m_HLT_Mu40_Eta2P1_PFJet200_PFJet50;
    delete m_HLT_Mu50;
    delete m_HLT_PFHT800;
    delete m_HLT_PFHT900;
    delete m_HLT_AK8PFJet450;
    delete m_HLT_PFHT700TrimMass50;
    delete m_HLT_PFJet360TrimMass30;

    delete m_Flag_goodVertices;
    delete m_Flag_eeBadScFilter;
    delete m_Flag_HBHENoiseFilter;
    delete m_Flag_HBHENoiseIsoFilter;
    delete m_Flag_globalTightHalo2016Filter;
    delete m_Flag_EcalDeadCellTriggerPrimitiveFilter;

    if (m_isMC){
      if (m_config->isTtbar()){
        delete m_mc_pt;
        delete m_mc_eta;
        delete m_mc_phi;
        delete m_mc_e;
        delete m_mc_pdgId;
        delete m_mc_status;
        delete m_mc_isHadTop;
      }
/*
        delete m_weight_mc;
        delete m_weight_pileup;
        delete m_weight_pileup_UP;
        delete m_weight_pileup_DOWN;

        delete m_mc_ht;

        delete m_truth_jet_pt;
        delete m_truth_jet_eta;
        delete m_truth_jet_phi;
        delete m_truth_jet_e;

        delete m_truth_ljet_pt;
        delete m_truth_ljet_eta;
        delete m_truth_ljet_phi;
        delete m_truth_ljet_m;
        delete m_truth_ljet_charge;
        delete m_truth_ljet_subjet0_charge;
        delete m_truth_ljet_subjet0_bdisc;
        delete m_truth_ljet_subjet1_charge;
        delete m_truth_ljet_subjet1_bdisc;
*/
    } // end isMC

    return;
}

// THE END

