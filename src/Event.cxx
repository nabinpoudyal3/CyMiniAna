/*
Created:        --
Last Updated:   22 August 2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109

-----

Event class
Contains all the objects (& structs) with event information

*/
#include "diHiggs/CyMiniAna/interface/Event.h"

// constructor
Event::Event( TTreeReader &myReader, configuration &cmaConfig ) :
  m_config(&cmaConfig),
  m_ttree(myReader),
  m_treeName("SetMe"),
  m_lwnn(nullptr),
  m_DNN(0.0),
  m_HME(0.0),
  m_amwt(nullptr){
    m_isMC     = m_config->isMC();
    m_grid     = m_config->isGridFile();             // file directly from original analysis team
    m_treeName = m_ttree.GetTree()->GetName();       // for systematics
    m_getDNN   = m_config->getDNN();                 // build DNN
    m_getHME   = m_config->getHME();                 // build HME
    m_buildNeutrinos = m_config->buildNeutrinos();   // build the neutrinos using AMWT

    m_cMVAv2L = m_config->cMVAv2L();
    m_cMVAv2M = m_config->cMVAv2M();
    m_cMVAv2T = m_config->cMVAv2T();

/*  // not available:
    m_mcChannelNumber = new TTreeReaderValue<unsigned int>(m_ttree,"mcChannelNumber");
    m_mu        = new TTreeReaderValue<float>(m_ttree,"mu");
    m_lumiblock = new TTreeReaderValue<int>(m_ttree,"lumiblock");
*/
    // now access branches we want to save or use
    // based names on files accessed from
    //     /fdata/hepx/store/user/tahuang/HHNTuples (brazos)
    // on 28 August 2017
    m_eventNumber  = new TTreeReaderValue<float>(m_ttree,"event_number");     // unsigned long long
    m_runNumber    = new TTreeReaderValue<float>(m_ttree,"event_run");        // unsigned int
    m_trigeff      = new TTreeReaderValue<float>(m_ttree,"trigeff");
    m_pu           = new TTreeReaderValue<float>(m_ttree,"pu");
    m_cosThetaStar = new TTreeReaderValue<float>(m_ttree,"cosThetaStar");
    m_isSF         = new TTreeReaderValue<bool>(m_ttree,"isSF");

    if (m_config->useLargeRJets()){
      // not available in this analysis (for now)
      m_ljet_pt     = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_pt" );
      m_ljet_eta    = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_eta");
      m_ljet_phi    = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_phi" );
      m_ljet_e      = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_e" );
      m_ljet_isGood = new TTreeReaderValue<std::vector<int>>(m_ttree,"ljet_isGood");
      m_ljet_d23       = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_d23" );
      m_ljet_tau1_wta  = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_tau1_wta" );
      m_ljet_tau2_wta  = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_tau2_wta" );
      m_ljet_tau3_wta  = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_tau3_wta" );
      m_ljet_tau21_wta = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_tau21_wta" );
      m_ljet_tau32_wta = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_tau32_wta" );
      m_ljet_charge    = new TTreeReaderValue<std::vector<float>>(m_ttree,"ljet_charge");
    }

    if (m_config->useJets()){
        m_jet1_pt      = new TTreeReaderValue<float>(m_ttree,"jet1_pt");
        m_jet1_eta     = new TTreeReaderValue<float>(m_ttree,"jet1_eta");
        m_jet1_phi     = new TTreeReaderValue<float>(m_ttree,"jet1_phi");
        m_jet1_cMVAv2  = new TTreeReaderValue<float>(m_ttree,"jet1_cMVAv2");
        m_jet2_pt      = new TTreeReaderValue<float>(m_ttree,"jet2_pt");
        m_jet2_eta     = new TTreeReaderValue<float>(m_ttree,"jet2_eta");
        m_jet2_phi     = new TTreeReaderValue<float>(m_ttree,"jet2_phi");
        m_jet2_cMVAv2  = new TTreeReaderValue<float>(m_ttree,"jet2_cMVAv2");

        m_nJetsL       = new TTreeReaderValue<float>(m_ttree,"nJetsL");
        m_jjbtag_heavy = new TTreeReaderValue<float>(m_ttree,"jjbtag_heavy");
        m_jjbtag_light = new TTreeReaderValue<float>(m_ttree,"jjbtag_light");
        m_jj_DR_j_j    = new TTreeReaderValue<float>(m_ttree,"jj_DR_j_j");
        m_jj_pt        = new TTreeReaderValue<float>(m_ttree,"jj_pt");
        m_jj_M         = new TTreeReaderValue<float>(m_ttree,"jj_M");
/*
      m_jet_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"jet_pt" );
      m_jet_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"jet_eta" );
      m_jet_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"jet_phi" );
      m_jet_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"jet_e" );
      m_jet_mv2c10 = new TTreeReaderValue<std::vector<float>>(m_ttree,"jet_mv2c10" );
      if (m_isMC)
        m_jet_true_flavor = new TTreeReaderValue<std::vector<int>>(m_ttree,"jet_truthflav" );
*/
    }


    if (m_config->useLeptons()){
      // lepton information in the grid
      m_lep1_pt  = new TTreeReaderValue<float>(m_ttree,"lep1_pt");
      m_lep1_eta = new TTreeReaderValue<float>(m_ttree,"lep1_eta");
      m_lep1_phi = new TTreeReaderValue<float>(m_ttree,"lep1_phi");
      m_lep1_Iso = new TTreeReaderValue<float>(m_ttree,"lep1_Iso");

      m_lep2_pt  = new TTreeReaderValue<float>(m_ttree,"lep2_pt");
      m_lep2_eta = new TTreeReaderValue<float>(m_ttree,"lep2_eta");
      m_lep2_phi = new TTreeReaderValue<float>(m_ttree,"lep2_phi");
      m_lep2_Iso = new TTreeReaderValue<float>(m_ttree,"lep2_Iso");

      m_ll_M        = new TTreeReaderValue<float>(m_ttree,"ll_M");
      m_ll_pt       = new TTreeReaderValue<float>(m_ttree,"ll_pt");
      m_ll_DR_l_l   = new TTreeReaderValue<float>(m_ttree,"ll_DR_l_l");
      m_ll_DPhi_l_l = new TTreeReaderValue<float>(m_ttree,"ll_DPhi_l_l");
      m_ll_DEta_l_l = new TTreeReaderValue<float>(m_ttree,"ll_DEta_l_l");
      m_llidiso     = new TTreeReaderValue<float>(m_ttree,"llidiso");
      m_mumuidiso   = new TTreeReaderValue<float>(m_ttree,"mumuidiso");
      m_elelidiso   = new TTreeReaderValue<float>(m_ttree,"elelidiso");

      m_isElEl   = new TTreeReaderValue<float>(m_ttree,"isElEl");
      m_isMuMu   = new TTreeReaderValue<float>(m_ttree,"isMuMu");
      m_isElMu   = new TTreeReaderValue<float>(m_ttree,"isElMu");
      m_isMuEl   = new TTreeReaderValue<float>(m_ttree,"isMuEl");


      if (m_config->useJets()){
        m_llmetjj_DPhi_ll_jj = new TTreeReaderValue<float>(m_ttree,"llmetjj_DPhi_ll_jj");
        m_llmetjj_minDR_l_j  = new TTreeReaderValue<float>(m_ttree,"llmetjj_minDR_l_j");
        m_llmetjj_MTformula  = new TTreeReaderValue<float>(m_ttree,"llmetjj_MTformula");
        m_llmetjj_MT2 = new TTreeReaderValue<float>(m_ttree,"llmetjj_MT2");
        m_llmetjj_M   = new TTreeReaderValue<float>(m_ttree,"llmetjj_M");
        m_lljj_M      = new TTreeReaderValue<float>(m_ttree,"lljj_M");
      }
/*
      m_el_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,   "el_pt");
      m_el_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,   "el_eta");
      m_el_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,   "el_phi");
      m_el_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,   "el_e");
      m_el_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"el_charge");
      m_mu_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,   "mu_pt");
      m_mu_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,   "mu_eta");
      m_mu_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,   "mu_phi");
      m_mu_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,   "mu_e");
      m_mu_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_charge");
*/
    }

    m_met_met = new TTreeReaderValue<float>(m_ttree,"met_pt");
    m_met_phi = new TTreeReaderValue<float>(m_ttree,"met_phi");
    m_ht      = new TTreeReaderValue<float>(m_ttree,"ht");

    // set some global weights and access necessary branches
    m_xsection       = 1.0;
    m_kfactor        = 1.0;
    m_sumOfWeights   = 1.0;
    m_LUMI           = m_config->LUMI();
    if (m_isMC){
      m_total_weight     = new TTreeReaderValue<float>(m_ttree,"total_weight");
      m_event_weight     = new TTreeReaderValue<float>(m_ttree,"event_weight");
      m_weight_mc        = new TTreeReaderValue<float>(m_ttree,"sample_weight");       //m_sample_weight
      m_weight_pileup    = new TTreeReaderValue<float>(m_ttree,"event_pu_weight");     //m_event_pu_weight

      // Get weights (global to all events; available as branches in Tao's samples)
      if (m_grid){  // samples from original analysis
        TParameter<double>* xsection = (TParameter<double>*)m_ttree.GetTree()->GetCurrentFile()->Get("cross_section");
        TParameter<double>* sum_of_weights = (TParameter<double>*)m_ttree.GetTree()->GetCurrentFile()->Get("event_weight_sum");
        m_xsection     = xsection->GetVal();       //m_config->XSectionMap( *(*m_mcChannelNumber) );
        m_kfactor      = 1.0;                      //m_config->KFactorMap( *(*m_mcChannelNumber) );
        m_sumOfWeights = sum_of_weights->GetVal(); //m_config->sumWeightsMap( *(*m_mcChannelNumber) );
      }
      else{
          m_treeXSection     = new TTreeReaderValue<float>(m_ttree,"cross_section");    // Tao's samples
          m_treeSumOfWeights = new TTreeReaderValue<float>(m_ttree,"event_weight_sum"); // Tao's samples
      }

      if (m_config->useTruth()){
        // not in this analysis, yet
        TTree* truth_tree = (TTree*)m_ttree.GetTree()->GetCurrentFile()->Get("truth");
        m_truth_tree.SetTree(truth_tree);
        m_truthEventNumber = new TTreeReaderValue<unsigned long long>(m_truth_tree,  "eventNumber");
        m_truthRunNumber   = new TTreeReaderValue<unsigned int>(m_truth_tree,"runNumber");
        m_truth_weight_mc  = new TTreeReaderValue<float>(m_truth_tree,"weight_mc");

        m_truth_ljet_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"truth_ljet_pt");
        m_truth_ljet_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"truth_ljet_eta");
        m_truth_ljet_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"truth_ljet_phi");
        m_truth_ljet_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"truth_ljet_e");
        m_truth_ljet_Qw  = new TTreeReaderValue<std::vector<float>>(m_ttree,"truth_ljet_Qw");
        m_truth_ljet_tau32_wta = new TTreeReaderValue<std::vector<float>>(m_ttree,"truth_ljet_tau32_wta");
        m_truth_ljet_split23   = new TTreeReaderValue<std::vector<float>>(m_ttree,"truth_ljet_split23");
      } // end useTruth
/*
      m_weight_btag_77   = new TTreeReaderValue<float>(m_ttree,"weight_bTagSF_77");
      m_weight_btag_70   = new TTreeReaderValue<float>(m_ttree,"weight_bTagSF_70");
      if (m_config->useLeptons())
        m_weight_lept_eff  = new TTreeReaderValue<float>(m_ttree,"weight_leptonSF");
      if ( m_config->isNominalTree( m_treeName ) ){
        initialize_eventWeights();
      } // end isNominal
*/
    } // end isMC
    m_truth_entry = 0;


    // HME material (TBD)
    if (m_getHME){
        m_mmc_hme = new MMC(cmaConfig);
        m_mmc_hme->initialize();
    }
    else if (!m_getHME && !m_grid)
        m_hme_h2mass_reco = new TTreeReaderValue<float>(m_ttree,"hme_h2mass_reco");
    else
        m_HME = 0.0;

    // DNN material
    bool useDNN(false);
    if (!m_getDNN && useDNN)  // always false for now
        m_dnn_score = new TTreeReaderValue<float>(m_ttree,"dnn_score");

    std::ifstream input_cfg = cma::open_file( m_config->dnnFile() );
    lwt::JSONConfig cfg     = lwt::parse_json( input_cfg );
    m_lwnn   = new lwt::LightweightNeuralNetwork(cfg.inputs, cfg.layers, cfg.outputs);
    m_dnnKey = m_config->dnnKey();

    // AMWT
    m_amwt = new AMWT(cmaConfig);
    m_amwt->initialize();
} // end constructor


Event::~Event() {}



void Event::initialize_eventWeights(){
    /* Create vectors of the systematics that are weights for the nominal events */
    std::map<std::string,unsigned int> mapWeightSystematics = m_config->mapOfWeightVectorSystematics();

    m_listOfWeightSystematics = m_config->listOfWeightSystematics();

    m_weightSystematicsFloats.clear();
    m_weightSystematicsVectorFloats.clear();

    // systematics from the nominal tree that are floats
    for (const auto& nom_syst : m_listOfWeightSystematics){
        if (!m_config->useLeptons() && nom_syst.find("leptonSF")!=std::string::npos)
            continue;
        m_weightSystematicsFloats[nom_syst] = new TTreeReaderValue<float>(m_ttree,nom_syst.c_str());
    }

    // systematics from the nominal tree that are vectors
    for (const auto& syst : mapWeightSystematics)
        m_weightSystematicsVectorFloats[syst.first] = new TTreeReaderValue<std::vector<float>>(m_ttree,syst.first.c_str());

    return;
}



void Event::updateEntry(Long64_t entry){
    /* Update the entry -> update all TTree variables */
    cma::DEBUG("EVENT : Update Entry "+std::to_string(entry) );
    if(m_config->matchTruthToReco()) {
        // this option means we loop over reco events and then match truth events (rather than the opposite)
        // legacy option from ATLAS work -- not needed for now
        m_entry = entry;
    }
    else {
        m_truth_tree.SetEntry(entry);
        m_truth_entry = entry;
        // try to find corresponding reco tree entry
        m_entry = m_ttree.GetTree()->GetEntryNumberWithIndex(**m_truthRunNumber, **m_truthEventNumber);
    }
    cma::DEBUG("EVENT : Set entry for updating ");

    // make sure the entry exists
    // when looping over truth events, this condition is not always met
    if(isValidRecoEntry())
        m_ttree.SetEntry(m_entry);
    else
        cma::ERROR("EVENT : Invalid Reco entry "+std::to_string(m_entry)+"!");

    return;
}


bool Event::isValidRecoEntry(){
    return (m_entry > (long long)-1);
}

void Event::clear(){
    /* Clear many of the vectors/maps for each event -- SAFETY PRECAUTION */
    m_truth_ljets.clear();
    m_truth_jets.clear();
    m_truth_leptons.clear();
    m_truth_neutrinos.clear();
    
    m_ljets.clear();
    m_jets.clear();
    m_leptons.clear();
    m_neutrinos.clear();

    m_btag_jets.clear();
    m_btag_jets_default.clear();
    m_weight_btag_default = 1.0;

    m_dilepton = {};
    m_dnnInputs.clear();

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

    // Get the event weights (for cutflow)
    initialize_weights();
    cma::DEBUG("EVENT : Setup weights ");

    // Truth Information
    if (m_config->useTruth() && m_config->isMC()){
        initialize_truth();
        cma::DEBUG("EVENT : Setup truth information ");
    }

    // Jets
    if (m_config->useJets()){
        initialize_jets();
        cma::DEBUG("EVENT : Setup small-R jets ");
    }

    // Large-R Jets
    if (m_config->useLargeRJets()){
        initialize_ljets();
        cma::DEBUG("EVENT : Setup large-R jets ");
    }

    // Leptons
    if (m_config->useLeptons()){
        initialize_leptons();
        cma::DEBUG("EVENT : Setup leptons ");
    }

    // Get some kinematics
    initialize_kinematics();
    cma::DEBUG("EVENT : Setup kinematic variables ");

    getDilepton();
    cma::DEBUG("EVENT : Setup the dilepton struct ");

    // Neutrinos
    if (m_config->useNeutrinos()){
        // Need ALL other information from the event to do this
        initialize_neutrinos();
        cma::DEBUG("EVENT : Setup neutrinos ");
    }

    // HME && DNN
    if (m_getHME){
        cma::DEBUG("EVENT : Calculate HME ");
        getHME();
    }
    else{
        // load from ntuple
        m_HME = *(*m_hme_h2mass_reco);
    }

    if (m_getDNN){
        cma::DEBUG("EVENT : Calculate DNN ");
        getDNN();
    }
    else{
        // load from ntuple -- not available yet!
        m_DNN = 0.0;  //*(*m_dnn_score);
    }

    cma::DEBUG("EVENT : Setup Event ");

    return;
}


void Event::initialize_truth(){
    /* Setup struct of truth information */
    return;
}


void Event::initialize_jets(){
    /* Setup struct of jets (small-r) and relevant information 
     * b-tagging: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
        cMVAv2L -0.5884	
        cMVAv2M 0.4432	 	 
        cMVAv2T 0.9432
     */
    m_jets.resize( **m_nJetsL );   // (*m_jet_pt)->size());
    m_btag_jets["L"].clear();
    m_btag_jets["M"].clear();
    m_btag_jets["T"].clear();

    Jet jet1;
    jet1.p4.SetPtEtaPhiM(*(*m_jet1_pt),*(*m_jet1_eta),*(*m_jet1_phi),4.18);  // use the b-quark mass [GeV]
    jet1.cMVAv2 = *(*m_jet1_cMVAv2);

    Jet jet2;
    jet2.p4.SetPtEtaPhiM(*(*m_jet2_pt),*(*m_jet2_eta),*(*m_jet2_phi),4.18);  // use the b-quark mass [GeV]
    jet2.cMVAv2 = *(*m_jet2_cMVAv2);

    m_jets[0] = jet1;
    m_jets[1] = jet2;

    // b-tagging
    if (jet1.cMVAv2 > m_cMVAv2L) m_btag_jets["L"].push_back(0);  // 0 = index of this jet
    if (jet1.cMVAv2 > m_cMVAv2M) m_btag_jets["M"].push_back(0);
    if (jet1.cMVAv2 > m_cMVAv2T) m_btag_jets["T"].push_back(0);
    if (jet2.cMVAv2 > m_cMVAv2L) m_btag_jets["L"].push_back(1);  // 1 = index of this jet
    if (jet2.cMVAv2 > m_cMVAv2M) m_btag_jets["M"].push_back(1);
    if (jet2.cMVAv2 > m_cMVAv2T) m_btag_jets["T"].push_back(1);

    m_btag_jets_default = m_btag_jets.at(m_config->jet_btagWkpt());

    // for JER
    calculateRho();         // decorate each jet with value rho

/*
    m_btag_jets["70"].clear();
    m_btag_jets["77"].clear();

    for (unsigned int i=0,size=(*m_jet_pt)->size(); i<size; i++){
        Jet jet;
        jet.p4.SetPtEtaPhiE((*m_jet_pt)->at(i),(*m_jet_eta)->at(i),(*m_jet_phi)->at(i),(*m_jet_e)->at(i));
        jet.mv2c10 = (*m_jet_mv2c10)->at(i);
        jet.mv2c20 = (*m_jet_mv2c20)->at(i);
        jet.jvt    = (*m_jet_jvt)->at(i);
        jet.true_flavor = (m_isMC) ? (*m_jet_true_flavor)->at(i) : -1;
        jet.isbtagged["70"] = (*m_jet_isbtagged_70)->at(i);
        jet.isbtagged["77"] = (*m_jet_isbtagged_77)->at(i);

        if ( (*m_jet_isbtagged_70)->at(i)==0x1 ) m_btag_jets["70"].push_back(i);
        if ( (*m_jet_isbtagged_77)->at(i)==0x1 ) m_btag_jets["77"].push_back(i);

        m_jets[i] = jet;
    }

    // store the configured b-tag WP information for fast+convenient access
    m_btag_jets_default = m_btag_jets.at(m_config->jet_btagWkpt());
*/
    return;
}


void Event::initialize_ljets(){
    /* Setup struct of large-R jets and relevant information */
    m_ljets.resize((*m_ljet_pt)->size());

    for (unsigned int i=0,size=(*m_ljet_pt)->size(); i<size; i++){
        // pre-selection (should be done in AnalysisTop?)
        // leading pT>500 GeV; others pT>350 GeV

        Ljet ljet;
        ljet.p4.SetPtEtaPhiE( (*m_ljet_pt)->at(i),(*m_ljet_eta)->at(i),(*m_ljet_phi)->at(i),(*m_ljet_e)->at(i));
        ljet.charge    = (*m_ljet_charge)->at(i);
        ljet.Split23   = (*m_ljet_d23)->at(i);
        ljet.tau1_wta  = (*m_ljet_tau1_wta)->at(i);
        ljet.tau2_wta  = (*m_ljet_tau2_wta)->at(i);
        ljet.tau3_wta  = (*m_ljet_tau3_wta)->at(i);
        ljet.tau21_wta = (*m_ljet_tau21_wta)->at(i);
        ljet.tau32_wta = (*m_ljet_tau32_wta)->at(i);
        ljet.isGood    = (ljet.p4.Pt()>200000. && fabs(ljet.p4.Eta())<2.0) ? 1 : 0;

        m_ljets[i] = ljet;
    }

    return;
}


void Event::initialize_leptons(){
    /* Setup struct of lepton and relevant information */
    m_ee   = ( *(*m_isElEl)>0.5 );
    m_mumu = ( *(*m_isMuMu)>0.5 );
    m_emu  = ( (*(*m_isElMu)>0.5) || (*(*m_isMuEl)>0.5) );

    m_leptons.resize(2);  // di-lepton analysis
    Lepton lep_1;
    Lepton lep_2;

    if (m_ee){
        lep_1.p4.SetPtEtaPhiM( *(*m_lep1_pt),*(*m_lep1_eta),*(*m_lep1_phi),0.0);  // assume no mass
        //lep_1.charge = (*m_el_charge)->at(0);
        lep_1.isElectron = true;
        lep_1.isMuon     = false;
        lep_1.Iso        = *(*m_lep1_Iso);

        lep_2.p4.SetPtEtaPhiM( *(*m_lep2_pt),*(*m_lep2_eta),*(*m_lep2_phi),0.0);  // assume no mass
        //lep_2.charge = (*m_el_charge)->at(0);
        lep_2.isElectron = true;
        lep_2.isMuon     = false;
        lep_2.Iso        = *(*m_lep2_Iso);
    }
    else if (m_mumu){
        lep_1.p4.SetPtEtaPhiM( *(*m_lep1_pt),*(*m_lep1_eta),*(*m_lep1_phi),0.0);  // assume no mass
        //lep_1.charge = (*m_el_charge)->at(0);
        lep_1.isElectron = false;
        lep_1.isMuon     = true;
        lep_1.Iso        = *(*m_lep1_Iso);

        lep_2.p4.SetPtEtaPhiM( *(*m_lep2_pt),*(*m_lep2_eta),*(*m_lep2_phi),0.0);  // assume no mass
        //lep_2.charge = (*m_el_charge)->at(0);
        lep_2.isElectron = false;
        lep_2.isMuon     = true;
        lep_2.Iso        = *(*m_lep2_Iso);
    }
    else{
        lep_1.p4.SetPtEtaPhiM( *(*m_lep1_pt),*(*m_lep1_eta),*(*m_lep1_phi),0.0);  // assume no mass
        //lep_1.charge = (*m_el_charge)->at(0);
        lep_1.isElectron = ( *(*m_isElMu)>0.5 );
        lep_1.isMuon     = ( *(*m_isMuEl)>0.5 );
        lep_1.Iso        = *(*m_lep1_Iso);

        lep_2.p4.SetPtEtaPhiM( *(*m_lep2_pt),*(*m_lep2_eta),*(*m_lep2_phi),0.0);  // assume no mass
        //lep_2.charge = (*m_el_charge)->at(0);
        lep_2.isElectron = ( *(*m_isMuEl)>0.5 );
        lep_2.isMuon     = ( *(*m_isElMu)>0.5 );
        lep_2.Iso        = *(*m_lep2_Iso);
    }

    m_leptons[0] = lep_1;
    m_leptons[1] = lep_2;

    return;
}


void Event::initialize_neutrinos(){
    /* Build the neutrinos */
    m_neutrinos.clear();

    if (m_buildNeutrinos){
        Neutrino nu1 = m_ttbar["top"].neutrino;
        Neutrino nu2 = m_ttbar["antitop"].neutrino;

        if (nu1.p4.Pt() > nu2.p4.Pt()){
            m_neutrinos.push_back( nu1 );
            m_neutrinos.push_back( nu2 );
        }
        else{
            m_neutrinos.push_back( nu2 );
            m_neutrinos.push_back( nu1 );
        }
    }
    else{
        // not supported yet
        Neutrino nu1;
        Neutrino nu2;
        m_neutrinos.push_back(nu1);
        m_neutrinos.push_back(nu2);

        m_ttbar = {};
    }

    return;
}


void Event::calculateRho(){
    /* Calculate jet energy density.
     * Using approximation: only two jets saved and don't know areas, assume 0.4
     *   -> rho = median( {jet.p4.Pt() / area}_i for 0<i<N_jets ) 
     */
    double rho(0.0);
    double area(0.4);  // approximation, not actual value!

    // Calculate jet_pt/area for each jet and store in vector
    std::vector<double> values;
    for (const auto& jet : m_jets){
        double value = jet.p4.Pt() / area;
        values.push_back(value);
    }

    // Get rho
    rho = cma::median<double>( values );

    // Store it as atttribute for each jet
    for (auto& jet : m_jets){
        jet.rho = rho;
    }

    return;
}


void Event::getDilepton(){
    /* Organize information into struct */
    // set dilepton struct
    m_dilepton = {};             // struct of information needed to build neutrinos

    // leptons
    Lepton lepton_p;
    Lepton lepton_n;
    if (m_leptons.at(0).charge > 0){
        lepton_p = m_leptons.at(0);
        lepton_n = m_leptons.at(1);
    }
    else{
        lepton_p = m_leptons.at(1);
        lepton_n = m_leptons.at(0);
    }

    m_dilepton.lepton_pos = lepton_p;
    m_dilepton.lepton_neg = lepton_n;

    // MET
    TVector2 dilep_met;
    dilep_met.SetX( m_metmet*cos(m_metphi) );
    dilep_met.SetY( m_metmet*sin(m_metphi) );
    m_dilepton.met = dilep_met;

    // Jets (only two 'b'-jets at this time)
    std::vector<Jet> dilep_jets;
    dilep_jets.push_back(m_jets.at(0));
    dilep_jets.push_back(m_jets.at(1));
    m_dilepton.jets = dilep_jets;

    cma::DEBUG("EVENT : Dilepton");
    cma::DEBUG("EVENT : met    = "+std::to_string(m_metmet));
    cma::DEBUG("EVENT : lepton = "+std::to_string(m_leptons.at(0).p4.Pt()));
    cma::DEBUG("EVENT : jet    = "+std::to_string(m_jets.at(0).p4.Pt()));

    return;
}


void Event::buildTtbar(){
    /* Build ttbar dilepton system */
    m_ttbar = m_amwt->findMass(m_dilepton);

    return;
}


void Event::getDNNInputs(){
    /* Load DNN values and return them */
    m_dnnInputs.clear();

    cma::DEBUG("EVENT : Set DNN input values ");

    return;
}


void Event::getDNN(){
    /* Dan Guest's lightweight DNN framework */
    getDNNInputs();

    std::map<std::string,double> discriminant = m_lwnn->compute(m_dnnInputs);
    m_DNN = discriminant.at(m_dnnKey);

    return;
}


double Event::DNN(){
    /* Return the DNN value */
    return m_DNN;
}


void Event::getHME(){
    /* Algorithm for Heavy Mass Estimator from Luca & Tao (their paper) */
    m_mmc_hme->execute( m_entry,m_dilepton );
    TH1F h_hme   = (TH1F)m_mmc_hme->getMMCh2();
    m_HME = (h_hme.GetXaxis())->GetBinCenter(h_hme.GetMaximumBin());

    return;
}


double Event::HME(){
    /* Return the HME value */
    return m_HME;
}


void Event::truth(){
    /* Do something with truth information (possibly change type and return information?) */
    return;
}


void Event::initialize_weights(){
    /* Event weights */
    m_nominal_weight = 1.0;

    m_weight_btag.clear();
    if (m_isMC){
        if (!m_grid){
            m_xsection     = **m_treeXSection;
            m_sumOfWeights = **m_treeSumOfWeights;
        }
        m_nominal_weight  = (**m_weight_pileup) * (**m_weight_mc);       // * (**m_weight_jvt);
        m_nominal_weight *= (m_xsection) * (m_kfactor) * m_LUMI / (m_sumOfWeights);
/*      // weights not in CMS (so far):
        m_weight_btag["70"] = (**m_weight_btag_70);
        m_weight_btag["77"] = (**m_weight_btag_77);
        m_weight_btag_default = m_weight_btag[m_config->jet_btagWkpt()];
        m_nominal_weight *= m_weight_btag_default;
        if (m_config->useLeptons())  m_nominal_weight *= (**m_weight_lept_eff);
*/
    }

    return;
}


void Event::initialize_kinematics(){
    /* Kinematic variables (HT, ST, MET) */
    m_metmet = *(*m_met_met);
    m_metphi = *(*m_met_phi);

    m_HT = *(*m_ht);   // total transverse hadronic energy
    m_ST = 0.0;        // total transverse energy

    m_ST += m_HT;

    m_ST = m_metmet;
    if (m_config->useLeptons()){
        for (const auto& lep : m_leptons)
            m_ST += lep.p4.Pt(); 
    }
/*
    // Not necessary for this analysis, yet
    if (m_config->useJets()){
        // include small-R jet pT
        for (auto &small_jet : m_jets ){
            m_ST += small_jet.p4.Pt();
            m_HT += small_jet.p4.Pt();
        }
    }
    else{
        // include large-R jet pT
        for (auto &large_jet : m_ljets){
            m_ST += large_jet.p4.Pt();
            m_HT += large_jet.p4.Pt();
        }
    }
*/
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
    else if (syst.find("jvt")!=std::string::npos){
        // pileup event weight
        syst_event_weight  = (**m_weight_pileup) * (**m_weight_mc);
        syst_event_weight *= m_weight_btag_default;
        if (m_config->useLeptons())  syst_event_weight *= (**m_weight_lept_eff);
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);

        syst_event_weight *= **m_weightSystematicsFloats.at(syst);
    }
    else if (syst.find("pileup")!=std::string::npos){
        // pileup event weight
        syst_event_weight  = (**m_weight_mc);
        syst_event_weight *= (**m_weight_jvt) * m_weight_btag_default;
        if (m_config->useLeptons())  syst_event_weight *= (**m_weight_lept_eff);
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);

        syst_event_weight *= **m_weightSystematicsFloats.at(syst);
    }
    else if (syst.find("leptonSF")!=std::string::npos){
        // leptonSF event weight
        syst_event_weight  = (**m_weight_pileup) * (**m_weight_mc);
        syst_event_weight *= (**m_weight_jvt) * m_weight_btag_default;
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);

        syst_event_weight *= **m_weightSystematicsFloats.at(syst);
    }
    else if (syst.find("bTagSF")!=std::string::npos){
        // bTagSF event weight -- check indices for eigenvector systematics
        syst_event_weight  = (**m_weight_pileup) * (**m_weight_mc);
        syst_event_weight *= (**m_weight_jvt);  syst_event_weight *= (**m_weight_lept_eff);
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);

        if (weightIndex>=0){
            syst_event_weight *= (*m_weightSystematicsVectorFloats.at(syst))->at(weightIndex);
        }
        else{
            syst_event_weight *= **m_weightSystematicsFloats.at(syst);
        }
    }
    else{
        // safety to catch something weird -- just return 1.0
        cma::WARNING("EVENT : Passed systematic variation, "+syst+", to Event::getSystEventWeight() ");
        cma::WARNING("EVENT : that is inconsistent with the CyMiniAna options of ");
        cma::WARNING("EVENT :     nominal, jvt, pileup, leptonSF, and bTagSF. ");
        cma::WARNING("EVENT : Returning a weight of 1.0. ");
        syst_event_weight = 1.0;
    }

    return syst_event_weight;
}


std::vector<Lepton> Event::leptons(){
    // Leptons
    return m_leptons;
}

std::vector<Lepton> Event::truth_leptons(){
    // Truth leptons
    return m_leptons;
}

std::vector<Neutrino> Event::neutrinos(){
    // Neutrinos
    return m_neutrinos;
}

std::vector<Neutrino> Event::truth_neutrinos(){
    // Truth neutrinos
    return m_truth_neutrinos;
}

std::vector<Ljet> Event::ljets(){
    // Large-R jets
    return m_ljets;
}

std::vector<Ljet> Event::truth_ljets(){
    // Truth Large-R jets
    return m_truth_ljets;
}

std::vector<Jet> Event::jets(){
    // Jets
    return m_jets;
}

std::vector<Jet> Event::truth_jets(){
    // Truth jets
    return m_truth_jets;
}

std::vector<int> Event::btag_jets(const std::string &wkpt){
    /* Small-R Jet b-tagging */
    std::string tmp_wkpt(wkpt);
    if(m_btag_jets.find(wkpt) == m_btag_jets.end()){
        cma::WARNING("EVENT : B-tagging working point "+wkpt+" does not exist.");
        cma::WARNING("EVENT : Return vector of b-tagged jets for default working point "+m_config->jet_btagWkpt());
        tmp_wkpt = m_config->jet_btagWkpt();
    }
    return m_btag_jets.at(tmp_wkpt);
}

std::vector<int> Event::btag_jets(){
    /* Small-R Jet b-tagging. Configured b-tag WP. */
    return m_btag_jets_default;
}

float Event::met( const std::string& met_name ){
    // MET
    float met_value(0.0);
    if (met_name.compare("met")==0)
        met_value = m_metmet;
    else if (met_name.compare("phi")==0)
        met_value = m_metphi;
    else{
        cma::WARNING("EVENT : Request for MET variable that is neither 'met' nor 'phi'");
        cma::WARNING("EVENT : Returning 0.0");
    }

    return met_value;
}

float Event::HT(){
    // Sum of hadronic transverse energy
    return m_HT;
}

float Event::ST(){
    // Sum of all transverse energy
    return m_ST;
}

std::map<std::string,Top> Event::ttbar(){
    return m_ttbar;
}

long long Event::entry(){
    return m_entry;
}


float Event::nominal_weight(){
    return m_nominal_weight;
}
float Event::truth_weight_mc(){
    return **m_truth_weight_mc;
}
float Event::weight_mc(){
    return **m_weight_mc;
}
float Event::weight_jvt(){
    return **m_weight_jvt;
}
float Event::weight_pileup(){
    return **m_weight_pileup;
}
float Event::weight_lept_eff(){
    return **m_weight_lept_eff;
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
    for (const auto& wsf : m_weightSystematicsVectorFloats)
        tmp_weightSystematicsVectorFloats[wsf.first] = **wsf.second;

    return tmp_weightSystematicsVectorFloats;
}

std::vector<std::string> Event::listOfWeightSystematics(){
    /* list of weight systematics */
    return m_listOfWeightSystematics;
}

std::string Event::treeName(){
    return m_treeName;
}

float Event::xsection(){
    return m_xsection;
}

float Event::kfactor(){
    return m_kfactor;
}

float Event::sumOfWeights(){
    return m_sumOfWeights;
}

//unsigned long long Event::eventNumber(){
float Event::eventNumber(){
    return **m_eventNumber;
}

//unsigned int Event::runNumber(){
float Event::runNumber(){
    return **m_runNumber;
}

unsigned int Event::mcChannelNumber(){
    return **m_mcChannelNumber;
}

float Event::mu(){
    return **m_mu;
}

int Event::lumiblock(){
    return **m_lumiblock;
}


void Event::finalize(){
    // delete variables
    cma::DEBUG("EVENT : Finalize() ");
    delete m_lwnn;
    delete m_amwt;
    delete m_mmc_hme;
    delete m_eventNumber;
    delete m_runNumber;
//    delete m_mcChannelNumber;
//    delete m_mu;
//    delete m_lumiblock;
    if (m_config->useLargeRJets()){
      delete m_ljet_pt;
      delete m_ljet_eta;
      delete m_ljet_phi;
      delete m_ljet_e;
      delete m_ljet_d23;
      delete m_ljet_tau1_wta;
      delete m_ljet_tau2_wta;
      delete m_ljet_tau3_wta;
      delete m_ljet_tau21_wta;
      delete m_ljet_tau32_wta;
    }
    if (m_config->useLeptons()){
      delete m_lep1_pt;
      delete m_lep1_eta;
      delete m_lep1_phi;
      delete m_lep1_Iso;
      delete m_lep2_pt;
      delete m_lep2_eta;
      delete m_lep2_phi;
      delete m_lep2_Iso;
      delete m_ll_M;
//      delete m_ll_pt;
      delete m_ll_DR_l_l;
      delete m_ll_DPhi_l_l;
      delete m_ll_DEta_l_l;
      delete m_llidiso;
      delete m_mumuidiso;
      delete m_elelidiso;
      delete m_isElEl;
      delete m_isMuMu;
      delete m_isElMu;
      delete m_isMuEl;
      if (m_config->useJets()){
        delete m_llmetjj_DPhi_ll_jj;
        delete m_llmetjj_minDR_l_j;
        delete m_llmetjj_MTformula;
        //delete m_llmetjj_MT2;
        delete m_llmetjj_M;
        //delete m_lljj_M;
      }
    }
    if (m_config->useJets()){
      delete m_jet1_pt;
      delete m_jet1_eta;
      delete m_jet1_phi;
      delete m_jet1_cMVAv2;
      delete m_jet2_pt;
      delete m_jet2_eta;
      delete m_jet2_phi;
      delete m_jet2_cMVAv2;
      delete m_nJetsL;
      delete m_jjbtag_heavy;
      delete m_jjbtag_light;
      delete m_jj_DR_j_j;
//      delete m_jj_pt;
      delete m_jj_M;

/*
      delete m_jet_e;
      delete m_jet_mv2c10;
      delete m_jet_mv2c20;
      delete m_jet_jvt;
      delete m_jet_isbtagged_70;
      delete m_jet_isbtagged_77;
      if (m_isMC){
        delete m_jet_true_flavor;
      }
*/
    }
    delete m_met_met;
    delete m_met_phi;
    delete m_ht;
    if (!m_getHME && !m_grid)
        delete m_hme_h2mass_reco;

/*
    delete m_treeXSection;
    delete m_treeKFactor;
    delete m_treeAMI;
*/
    if (m_isMC){
      delete m_trigeff;
      delete m_pu;
      delete m_cosThetaStar;
      delete m_isSF;
      delete m_weight_mc;
      delete m_weight_pileup;
/*
      delete m_weight_jvt;
      delete m_weight_btag_77;
      delete m_weight_btag_70;
      if (m_config->useLeptons())
        delete m_weight_lept_eff;
      if ( m_config->isNominalTree( m_ttree.GetName() ) ){
        for ( auto& x: m_weightSystematicsFloats )
            delete x.second;
        m_weightSystematicsFloats.clear();
        for ( auto& x: m_weightSystematicsVectorFloats )
            delete x.second;
        m_weightSystematicsVectorFloats.clear();
      }
      if (m_config->useTruth()){
        delete m_truth_ljet_pt;
        delete m_truth_ljet_eta;
        delete m_truth_ljet_phi;
        delete m_truth_ljet_e;
        delete m_truth_ljet_Qw;
        delete m_truth_ljet_split23;
        delete m_truth_ljet_tau32_wta;
      } // end useTruth
*/
    } // end isMC

    return;
}

// THE END
