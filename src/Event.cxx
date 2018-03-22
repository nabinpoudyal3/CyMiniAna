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
  m_fileName("SetMe"),
  m_dileptonTtbar(nullptr){
    m_isMC     = m_config->isMC();
    m_grid     = m_config->isGridFile();             // file directly from EDM->FlatNtuple step
    m_treeName = m_ttree.GetTree()->GetName();       // for systematics
    m_fileName = m_config->filename();               // for accessing file metadata

    // Neural network and kinematic reconstruction information
    m_isZeroLeptonAnalysis = m_config->isZeroLeptonAnalysis();
    m_isOneLeptonAnalysis  = m_config->isOneLeptonAnalysis();
    m_isTwoLeptonAnalysis  = m_config->isTwoLeptonAnalysis();

    m_useJets       = m_config->useJets();
    m_useLargeRJets = m_config->useLargeRJets();
    m_useLeptons    = m_config->useLeptons();
    m_useNeutrinos  = m_config->useNeutrinos();
    m_neutrinoReco  = m_config->neutrinoReco();            // reconstruct neutrino
    m_DNNinference  = m_config->DNNinference();            // use DNN to predict values
    m_DNNtraining   = m_config->DNNtraining();             // load DNN features (save/use later)
    m_getDNN        = (m_DNNinference || m_DNNtraining);   // CWoLa
    m_kinematicReco = m_config->kinematicReco();           // build the ttbar system

    // b-tagging working points
    m_CSVv2L = m_config->CSVv2L();
    m_CSVv2M = m_config->CSVv2M();
    m_CSVv2T = m_config->CSVv2T();


    //** Access branches from Tree **//
    m_eventNumber  = new TTreeReaderValue<unsigned int>(m_ttree,"eventNum");
    m_runNumber    = new TTreeReaderValue<unsigned int>(m_ttree,"runNum");
    m_lumiblock    = new TTreeReaderValue<unsigned int>(m_ttree,"lumiNum");

/*
    m_HLT_Ele45_WPLoose_Gsf          = new TTreeReaderValue<int>(m_ttree,"HLT_Ele45_WPLoose_Gsf");
    m_HLT_Ele45_WPLoose_Gsf_prescale = new TTreeReaderValue<int>(m_ttree,"HLT_Ele45_WPLoose_Gsf_prescale");
    m_HLT_Mu50            = new TTreeReaderValue<int>(m_ttree,"HLT_Mu50");
    m_HLT_Mu50_prescale   = new TTreeReaderValue<int>(m_ttree,"HLT_Mu50_prescale");
    m_HLT_TkMu50          = new TTreeReaderValue<int>(m_ttree,"HLT_TkMu50");
    m_HLT_TkMu50_prescale = new TTreeReaderValue<int>(m_ttree,"HLT_TkMu50_prescale");
*/

    /** JETS **/
    if (m_useJets){
      // small-R jet information
      m_jet_pt   = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4pt");
      m_jet_eta  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4eta");
      m_jet_phi  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4phi");
      m_jet_m    = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4mass");
      m_jet_bdisc  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK4bDisc");
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
      m_ljet_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8charge");
      m_ljet_subjet0_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8ChargeSubjet1");
      m_ljet_subjet0_bdisc  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8bDiscSubjet1");
      m_ljet_subjet1_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8ChargeSubjet2");
      m_ljet_subjet1_bdisc  = new TTreeReaderValue<std::vector<float>>(m_ttree,"AK8bDiscSubjet2");

      m_ljet_BEST_class = new TTreeReaderValue<std::vector<float>>(m_ttree,"BESTclass");
      m_ljet_BEST_t = new TTreeReaderValue<std::vector<float>>(m_ttree,"BESTProb_t");
      m_ljet_BEST_w = new TTreeReaderValue<std::vector<float>>(m_ttree,"BESTProb_W");
      m_ljet_BEST_z = new TTreeReaderValue<std::vector<float>>(m_ttree,"BESTProb_Z");
      m_ljet_BEST_h = new TTreeReaderValue<std::vector<float>>(m_ttree,"BESTProb_H");
      m_ljet_BEST_j = new TTreeReaderValue<std::vector<float>>(m_ttree,"BESTProb_j");
    }


    /** LEPTONS **/
    if (m_useLeptons){
      m_el_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELpt");
      m_el_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELeta");
      m_el_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELphi");
      m_el_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELenergy");
      m_el_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELcharge");
      m_el_iso = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELiso");
      m_el_id  = new TTreeReaderValue<std::vector<float>>(m_ttree,"ELid");

      m_mu_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUpt");
      m_mu_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUeta");
      m_mu_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUphi");
      m_mu_e   = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUenergy");
      m_mu_charge = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUcharge");
      m_mu_iso = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUcorrIso");
      m_mu_id  = new TTreeReaderValue<std::vector<float>>(m_ttree,"MUlooseID");
    }

    if (!m_kinematicReco && m_useNeutrinos){
        // Neutrinos aren't stored in the baseline ntuples, requires 'kinematicReco' to create
        m_nu_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree, "nu_pt");
        m_nu_eta = new TTreeReaderValue<std::vector<float>>(m_ttree, "nu_eta");
        m_nu_phi = new TTreeReaderValue<std::vector<float>>(m_ttree, "nu_phi");
    }

    if (!m_kinematicReco && m_getDNN){
        // Load ttbar variables from file
        m_leptop_jet  = new TTreeReaderValue<int>(m_ttree, "leptop_jet");
        m_hadtop_ljet = new TTreeReaderValue<int>(m_ttree, "hadtop_ljet");
    }

    m_met_met  = new TTreeReaderValue<float>(m_ttree,"METpt");
    m_met_phi  = new TTreeReaderValue<float>(m_ttree,"METphi");

    // set some event weights and access necessary branches
    m_xsection       = 1.0;
    m_kfactor        = 1.0;
    m_sumOfWeights   = 1.0;
    m_LUMI           = m_config->LUMI();

    // MC information
    if (m_isMC){
//      m_weight_mc    = 1;//new TTreeReaderValue<float>(m_ttree,"evt_Gen_Weight");
      m_xsection     = 1;//m_config->XSectionMap( m_fileName );
      m_kfactor      = 1;//m_config->KFactorMap(  m_fileName );
      m_sumOfWeights = 1;//m_config->sumWeightsMap( m_fileName );
/*
      m_mc_ht = new TTreeReaderValue<float>(m_ttree,"evt_Gen_Ht");
      m_mc_pdgId;

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
    if (!m_getDNN && m_useDNN)  // always false for now
        m_dnn_score = new TTreeReaderValue<float>(m_ttree,"ljet_CWoLa");


    // Kinematic reconstruction algorithms
    m_ttbarRecoTool = new ttbarReco(cmaConfig);
    // m_semileptonTtbar = new semileptonTtbarReco(m_config);  // semi-leptonic ttbar kinematic reco
    // m_allhadTtbar = new allhadTtbarReco(m_config);          // all-hadronic ttbar kinematic reco
    //m_dileptonTtbar = new dileptonTtbarReco(cmaConfig, configuration::run2_13tev_2016_25ns, 2, true);
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
        if (!m_config->useLeptons() && nom_syst.find("leptonSF")!=std::string::npos)
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

    m_dilepton = {};

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

    // Truth Information
    if (m_config->useTruth() && m_config->isMC()){
        initialize_truth();
        cma::DEBUG("EVENT : Setup truth information ");
    }

    // Leptons
    if (m_config->useLeptons()){
        initialize_leptons();
        cma::DEBUG("EVENT : Setup leptons ");
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

    // Get some kinematic variables (MET, HT, ST)
    initialize_kinematics();
    cma::DEBUG("EVENT : Setup kinematic variables ");


    // Neutrinos
    if (m_config->useNeutrinos()){
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
        if (m_isOneLeptonAnalysis){
            int ljetidx = **m_hadtop_ljet;
            int jetidx  = **m_leptop_jet;
            m_ttbar1L.jet  = m_jets.at(jetidx);
            m_ttbar1L.ljet = m_ljets.at(ljetidx);
            if (m_muons.size()>0)
                m_ttbar1L.lepton = m_muons.at(0);
            else if (m_electrons.size()>0)
                m_ttbar1L.lepton = m_electrons.at(0);
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
    else if (m_isOneLeptonAnalysis){
        m_ttbarRecoTool->execute(m_electrons,m_muons,m_jets,m_ljets);
        m_ttbar1L = m_ttbarRecoTool->ttbar1L();
    }
    else if (m_isTwoLeptonAnalysis){
        m_ttbarRecoTool->execute(m_electrons,m_muons,m_jets);
        m_ttbar2L = m_ttbarRecoTool->ttbar2L();
    }

    return;
}


void Event::initialize_truth(){
    /* Setup struct of truth information */
    return;
}


void Event::initialize_jets(){
    /* Setup struct of jets (small-r) and relevant information 
     * b-tagging: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
        CSVv2L -0.5884
        CSVv2M 0.4432
        CSVv2T 0.9432
     */
    unsigned int nJets = (*m_jet_pt)->size();
    m_jets.clear();

    for (const auto& btagWP : m_config->btagWkpts() ){
        m_btag_jets[btagWP].clear();
    }

    unsigned int idx(0);
    for (unsigned int i=0; i<nJets; i++){
        Jet jet;
        jet.p4.SetPtEtaPhiM( (*m_jet_pt)->at(i),(*m_jet_eta)->at(i),(*m_jet_phi)->at(i),(*m_jet_m)->at(i));

        bool isGood(jet.p4.Pt()>50 && std::abs(jet.p4.Eta())<2.4);
        if (!isGood) continue;

        jet.bdisc  = (*m_jet_bdisc)->at(i);
        jet.index  = idx;
        jet.isGood = isGood;

        getBtaggedJets(jet);

        m_jets.push_back(jet);
        idx++;
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
    if (m_config->isOneLeptonAnalysis() && (m_electrons.size()+m_muons.size())>0){
        int charge = (m_electrons.size()>0) ? m_electrons.at(0).charge : m_muons.at(0).charge;
        target = (charge>0) ? 1:0;
    }

    unsigned int idx(0);
    for (unsigned int i=0; i<nLjets; i++){
        Ljet ljet;
        ljet.p4.SetPtEtaPhiM( (*m_ljet_pt)->at(i),(*m_ljet_eta)->at(i),(*m_ljet_phi)->at(i),(*m_ljet_m)->at(i));
        ljet.softDropMass = (*m_ljet_SDmass)->at(i);

        bool isGood(ljet.p4.Pt()>400. && fabs(ljet.p4.Eta())<2.4 && ljet.softDropMass>10.);
        if (!isGood) continue;

        ljet.charge = (*m_ljet_charge)->at(i);
        ljet.tau1   = (*m_ljet_tau1)->at(i);
        ljet.tau2   = (*m_ljet_tau2)->at(i);
        ljet.tau3   = (*m_ljet_tau3)->at(i);
        ljet.tau21  = ljet.tau2 / ljet.tau1;
        ljet.tau32  = ljet.tau3 / ljet.tau2;

        ljet.BEST_t = (*m_ljet_BEST_t)->at(i);
        ljet.BEST_w = (*m_ljet_BEST_w)->at(i);
        ljet.BEST_z = (*m_ljet_BEST_z)->at(i);
        ljet.BEST_h = (*m_ljet_BEST_h)->at(i);
        ljet.BEST_j = (*m_ljet_BEST_j)->at(i);
        ljet.BEST_class = (*m_ljet_BEST_class)->at(i);

        ljet.subjet0_bdisc  = (*m_ljet_subjet0_bdisc)->at(i);
        ljet.subjet0_charge = (*m_ljet_subjet0_charge)->at(i);
        ljet.subjet1_bdisc  = (*m_ljet_subjet1_bdisc)->at(i);
        ljet.subjet1_charge = (*m_ljet_subjet1_charge)->at(i);

        ljet.target = target;
        ljet.isGood = isGood;
        ljet.index  = idx;

        if (isGood && ljet.p4.Pt()<400) std::cout << " LARGE R JET SNUCK THROUGH " << std::endl;

        m_ljets.push_back(ljet);
        idx++;
    }

    deepLearningPrediction();   // store features in map (easily access later)

    return;
}


void Event::initialize_leptons(){
    /* Setup struct of lepton and relevant information */
    m_ee   = false;
    m_mumu = false;
    m_emu  = false;

    m_leptons.clear();
    m_electrons.clear();
    m_muons.clear();

    // Muons
    unsigned int nMuons = (*m_mu_pt)->size();
    m_muons.clear();

    for (unsigned int i=0; i<nMuons; i++){
        Muon mu;
        mu.p4.SetPtEtaPhiE( (*m_mu_pt)->at(i),(*m_mu_eta)->at(i),(*m_mu_phi)->at(i),(*m_mu_e)->at(i));

        bool isGood(mu.p4.Pt()>50 && std::abs(mu.p4.Eta())<2.1);
        if (!isGood) continue;

        mu.charge = (*m_mu_charge)->at(i);
        mu.id  = (*m_mu_id)->at(i);
        mu.iso = (*m_mu_iso)->at(i);
        mu.isGood = isGood;

        m_muons.push_back(mu);
    }

    // Electrons
    unsigned int nElectrons = (*m_el_pt)->size();
    m_electrons.clear();
    for (unsigned int i=0; i<nElectrons; i++){
        Electron el;
        el.p4.SetPtEtaPhiE( (*m_el_pt)->at(i),(*m_el_eta)->at(i),(*m_el_phi)->at(i),(*m_el_e)->at(i));

        bool isGood(el.p4.Pt()>50 && std::abs(el.p4.Eta())<2.1);
        if (!isGood) continue;

        el.charge = (*m_el_charge)->at(i);
        el.id  = (*m_el_id)->at(i);
        el.iso = (*m_el_iso)->at(i);
        el.isGood = isGood;

        m_electrons.push_back(el);
    }

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
    Neutrino nu2;

    if (m_neutrinoReco){}  // reconstruct neutrinos!
        if (m_isOneLeptonAnalysis){
        }
        else if (m_isTwoLeptonAnalysis){
        }
    else{
        // Assign neutrinos from root file
    }

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
    m_HT = 0.0;
    m_ST = 0.0;

    // Get hadronic transverse energy
    if (m_config->useJets()){
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

    // Get MET and lepton transverse energy
    m_ST += m_HT;
    m_ST += m_met.p4.Pt();

    if (m_config->useLeptons()){
        for (const auto& lep : m_leptons)
            m_ST += lep.p4.Pt(); 
    }

    return;
}




/*** GETTER FUNCTIONS ***/
void Event::getDilepton(){
    /* Organize information into struct */
    m_dilepton = {};             // struct of information needed to build neutrinos

    // Leptons
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
    dilep_met.SetX( m_met.p4.X() );
    dilep_met.SetY( m_met.p4.Y() );
    m_dilepton.met = dilep_met;

    // Jets
    std::vector<Jet> bjets;
    for (const auto& j : m_btag_jets_default){
        bjets.push_back(m_jets.at(j));
    }
    m_dilepton.jets  = m_jets;
    m_dilepton.bjets = bjets;

    cma::DEBUG("EVENT : Dilepton");
    cma::DEBUG("EVENT : met       = "+std::to_string(m_met.p4.Pt()));
    cma::DEBUG("EVENT : lepton pT = "+std::to_string(m_leptons.at(0).p4.Pt()));
    cma::DEBUG("EVENT : jet pT    = "+std::to_string(m_jets.at(0).p4.Pt()));

    return;
}


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



/*** RETURN EVENT INFORMATION ***/
void Event::truth(){
    /* Do something with truth information (possibly change type and return information?) */
    return;
}



/*** DELETE VARIABLES ***/
void Event::finalize(){
    // delete variables
    cma::DEBUG("EVENT : Finalize() ");
    //delete m_dileptonTtbar;
    delete m_eventNumber;
    delete m_runNumber;
    delete m_lumiblock;

    if (m_config->useJets()){
      delete m_jet_pt;
      delete m_jet_eta;
      delete m_jet_phi;
      delete m_jet_m;
      delete m_jet_bdisc;
    }

    if (m_config->useLargeRJets()){
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
      delete m_ljet_subjet1_charge;
      delete m_ljet_subjet1_bdisc;
    }

    if (m_config->useLeptons()){
      delete m_el_pt;
      delete m_el_eta;
      delete m_el_phi;
      delete m_el_e;
      delete m_el_charge;
      delete m_el_iso;
      delete m_el_id;

      delete m_mu_pt;
      delete m_mu_eta;
      delete m_mu_phi;
      delete m_mu_e;
      delete m_mu_charge;
      delete m_mu_iso;
      delete m_mu_id;
    }

    delete m_met_met;          // met_Pt
    delete m_met_phi;          // met_Phi

/*
    delete m_HLT_Ele45_WPLoose_Gsf;
    delete m_HLT_Mu50;
    delete m_HLT_TkMu50;
*/
    if (m_isMC){
      //delete m_weight_mc;
      //delete m_weight_pileup;
      //delete m_weight_pileup_UP;
      //delete m_weight_pileup_DOWN;

      if (m_config->useTruth()){
/*        delete m_mc_ht;

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
*/      } // end useTruth
    } // end isMC

    return;
}

// THE END
