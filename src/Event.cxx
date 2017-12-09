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
#include "cms-ttbarAC/CyMiniAna/interface/Event.h"

// constructor
Event::Event( TTreeReader &myReader, configuration &cmaConfig ) :
  m_config(&cmaConfig),
  m_ttree(myReader),
  m_treeName("SetMe"),
  m_fileName("SetMe"),
  m_lwnn(nullptr),
  m_DNN(0.0),
  m_dileptonTtbar(nullptr){
    m_isMC     = m_config->isMC();
    m_grid     = m_config->isGridFile();             // file directly from original analysis team
    m_treeName = m_ttree.GetTree()->GetName();       // for systematics
    m_fileName = m_config->filename();               // for accessing file metadata

    // Neural network and kinematic reconstruction information
    m_isZeroLeptonAnalysis = m_config->isZeroLeptonAnalysis();
    m_isOneLeptonAnalysis  = m_config->isOneLeptonAnalysis();
    m_isTwoLeptonAnalysis  = false; // m_config->isTwoLeptonAnalysis(); (handled elsewhere)

    m_getDNN        = m_config->getDNN();            // build DNN
    m_kinematicReco = m_config->kinematicReco();     // build the ttbar system

    // b-tagging working points
    m_cMVAv2L = m_config->cMVAv2L();
    m_cMVAv2M = m_config->cMVAv2M();
    m_cMVAv2T = m_config->cMVAv2T();


    //** Access branches from Tree **//
    m_eventNumber  = new TTreeReaderValue<int>(m_ttree,"evt_EventNumber");
    m_runNumber    = new TTreeReaderValue<int>(m_ttree,"evt_RunNumber");
    m_rho          = new TTreeReaderValue<double>(m_ttree,"evt_rho");
    m_lumiblock    = new TTreeReaderValue<int>(m_ttree,"evt_LumiBlock");
    m_NGoodVtx     = new TTreeReaderArray<int>(m_ttree,"evt_NGoodVtx");
    m_LHAPDF_ID    = new TTreeReaderArray<int>(m_ttree,"evt_LHA_PDF_ID");
    m_NIsoTrk      = new TTreeReaderArray<int>(m_ttree,"evt_NIsoTrk");
    m_pu_NtrueInt  = new TTreeReaderArray<int>(m_ttree,"pu_NtrueInt");

    m_HLT_Ele45_WPLoose_Gsf          = new TTreeReaderArray<int>(m_ttree,"HLT_Ele45_WPLoose_Gsf");
    m_HLT_Ele45_WPLoose_Gsf_prescale = new TTreeReaderArray<int>(m_ttree,"HLT_Ele45_WPLoose_Gsf_prescale");
    m_HLT_Mu50            = new TTreeReaderArray<int>(m_ttree,"HLT_Mu50");
    m_HLT_Mu50_prescale   = new TTreeReaderArray<int>(m_ttree,"HLT_Mu50_prescale");
    m_HLT_TkMu50          = new TTreeReaderArray<int>(m_ttree,"HLT_TkMu50");
    m_HLT_TkMu50_prescale = new TTreeReaderArray<int>(m_ttree,"HLT_TkMu50_prescale");

    m_scale_size     = new TTreeReaderValue<unsigned int>(m_ttree,"scale_size");
    m_scale_Weights  = new TTreeReaderArray<float>(m_ttree,"scale_Weights");
    m_pdf_size       = new TTreeReaderValue<unsigned int>(m_ttree,"pdf_size");
    m_weights_pdf    = new TTreeReaderArray<float>(m_ttree,"pdf_Weights");
    m_alphas_size    = new TTreeReaderValue<unsigned int>(m_ttree,"alphas_size");
    m_weights_alphas = new TTreeReaderArray<float>(m_ttree,"alphas_Weights");

    if (m_config->useFlags()){
      m_Flag_BadChargedCandidateFilter = new TTreeReaderArray<int>(m_ttree,"Flag_BadChargedCandidateFilter");
      m_Flag_HBHENoiseFilter           = new TTreeReaderArray<int>(m_ttree,"Flag_HBHENoiseFilter");
      m_Flag_HBHENoiseIsoFilter        = new TTreeReaderArray<int>(m_ttree,"Flag_HBHENoiseIsoFilter");
      m_Flag_CSCTightHaloFilter        = new TTreeReaderArray<int>(m_ttree,"Flag_CSCTightHaloFilter");
      m_Flag_HcalStripHaloFilter  = new TTreeReaderArray<int>(m_ttree,"Flag_HcalStripHaloFilter");
      m_Flag_hcalLaserEventFilter = new TTreeReaderArray<int>(m_ttree,"Flag_hcalLaserEventFilter");
      m_Flag_goodVertices         = new TTreeReaderArray<int>(m_ttree,"Flag_goodVertices");
      m_Flag_eeBadScFilter        = new TTreeReaderArray<int>(m_ttree,"Flag_eeBadScFilter");
      m_Flag_ecalLaserCorrFilter  = new TTreeReaderArray<int>(m_ttree,"Flag_ecalLaserCorrFilter");
      m_Flag_CSCTightHaloTrkMuUnvetoFilter  = new TTreeReaderArray<int>(m_ttree,"Flag_CSCTightHaloTrkMuUnvetoFilter");
      m_Flag_CSCTightHalo2015Filter         = new TTreeReaderArray<int>(m_ttree,"Flag_CSCTightHalo2015Filter");
      m_Flag_globalTightHalo2016Filter      = new TTreeReaderArray<int>(m_ttree,"Flag_globalTightHalo2016Filter");
      m_Flag_globalSuperTightHalo2016Filter = new TTreeReaderArray<int>(m_ttree,"Flag_globalSuperTightHalo2016Filter");
      m_Flag_EcalDeadCellTriggerPrimitiveFilter = new TTreeReaderArray<int>(m_ttree,"Flag_EcalDeadCellTriggerPrimitiveFilter");
      m_Flag_EcalDeadCellBoundaryEnergyFilter   = new TTreeReaderArray<int>(m_ttree,"Flag_EcalDeadCellBoundaryEnergyFilter");
      m_Flag_chargedHadronTrackResolutionFilter = new TTreeReaderArray<int>(m_ttree,"Flag_chargedHadronTrackResolutionFilter");

      m_Flag_trkPOGFilters                  = new TTreeReaderArray<int>(m_ttree,"Flag_trkPOGFilters");
      m_Flag_trkPOG_manystripclus53X        = new TTreeReaderArray<int>(m_ttree,"Flag_trkPOG_manystripclus53X");
      m_Flag_trkPOG_toomanystripclus53X     = new TTreeReaderArray<int>(m_ttree,"Flag_trkPOG_toomanystripclus53X");
      m_Flag_trkPOG_logErrorTooManyClusters = new TTreeReaderArray<int>(m_ttree,"Flag_trkPOG_logErrorTooManyClusters");

      m_Flag_METFilters = new TTreeReaderArray<int>(m_ttree,"Flag_METFilters");

      m_Flag_BadPFMuonFilter = new TTreeReaderArray<int>(m_ttree,"Flag_BadPFMuonFilter");
      m_Flag_badMuons        = new TTreeReaderArray<int>(m_ttree,"Flag_badMuons");
      m_Flag_muonBadTrackFilter = new TTreeReaderArray<int>(m_ttree,"Flag_muonBadTrackFilter");
      m_Flag_duplicateMuons  = new TTreeReaderArray<int>(m_ttree,"Flag_duplicateMuons");
      m_Flag_noBadMuons      = new TTreeReaderArray<int>(m_ttree,"Flag_noBadMuons");
    }


    /** JETS **/
    if (m_config->useJets()){
      // small-R jet information
      m_jet_size = new TTreeReaderValue<unsigned int>(m_ttree,"jetAK4CHS_size");
      m_jet_pt   = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_Pt");
      m_jet_eta  = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_Eta");
      m_jet_phi  = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_Phi");
      m_jet_e    = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_E");
      m_jet_charge = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_Charge");
      m_jet_CSVv2  = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_CSVv2");
      m_jet_CMVAv2 = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_CMVAv2");
      m_jet_CvsL   = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_CvsL");
      m_jet_CvsB   = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_CvsB");
      m_jet_partonFlavour = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_PartonFlavour");
      m_jet_hadronFlavour = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_HadronFlavour");
      m_jet_neutralMultiplicity     = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_neutralMultiplicity");
      m_jet_neutralHadronEnergyFrac = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_neutralHadronEnergyFrac");
      m_jet_neutralEmEnergyFrac     = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_neutralEmEnergyFrac");
      m_jet_chargedHadronEnergyFrac = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_chargedHadronEnergyFrac");
      m_jet_chargedEmEnergyFrac     = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_chargedEmEnergyFrac");
      m_jet_chargedMultiplicity     = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_chargedMultiplicity");
      m_jet_jecFactor0 = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_jecFactor0");
      m_jet_jetArea    = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_jetArea");
      m_jet_jecUncertainty = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_jecUncertainty");
      m_jet_ptResolution   = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_PtResolution");
      m_jet_JERSF     = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_JERSF");
      m_jet_JERSFUp   = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_JERSFUp");
      m_jet_JERSFDown = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_JERSFDown");
      m_jet_smearedPt = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_SmearedPt");
    }
    if (m_config->useLargeRJets()){
      // large-R Jet information
      m_ljet_size = new TTreeReaderValue<unsigned int>(m_ttree,"jetAK8CHS_size");
      m_ljet_pt   = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_Pt");
      m_ljet_eta  = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_Eta");
      m_ljet_phi  = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_Phi");
      m_ljet_e    = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_E");
      m_ljet_softDropMass_CHS = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_softDropMassCHS");
      m_ljet_tau1_CHS = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_tau1CHS");
      m_ljet_tau2_CHS = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_tau1CHS");
      m_ljet_tau3_CHS = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_tau1CHS");
      m_ljet_charge = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_Charge");
      m_ljet_CSVv2  = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_CSVv2");
      m_ljet_CMVAv2 = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_CMVAv2");
      m_ljet_CvsL   = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_CvsL");
      m_ljet_CvsB   = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_CvsB");
      m_ljet_PartonFlavour = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_PartonFlavour");
      m_ljet_HadronFlavour = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_HadronFlavour");
      m_ljet_neutralMultiplicity     = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_neutralMultiplicity");
      m_ljet_neutralHadronEnergyFrac = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_neutralHadronEnergyFrac");
      m_ljet_neutralEmEnergyFrac     = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_neutralEmEnergyFrac");
      m_ljet_chargedHadronEnergyFrac = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_chargedHadronEnergyFrac");
      m_ljet_chargedEmEnergyFrac = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_chargedEmEnergyFrac");
      m_ljet_chargedMultiplicity = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_chargedMultiplicity");
      m_ljet_jecFactor0 = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_jecFactor0");
      m_ljet_jetArea    = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_jetArea");
      m_ljet_jecUncertainty = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_jecUncertainty");
      m_ljet_PtResolution   = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_PtResolution");
      m_ljet_JERSF     = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_JERSF");
      m_ljet_JERSFUp   = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_JERSFUp");
      m_ljet_JERSFDown = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_JERSFDown");
      m_ljet_SmearedPt = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_SmearedPt");
      m_ljet_vSubjetIndex0 = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_vSubjetIndex0");
      m_ljet_vSubjetIndex1 = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_vSubjetIndex1");

      m_ljet_subjet_size = new TTreeReaderValue<unsigned int>(m_ttree,"subjetAK8CHS_size");
      m_ljet_subjet_pt   = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_Pt");
      m_ljet_subjet_eta  = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_Eta");
      m_ljet_subjet_phi  = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_Phi");
      m_ljet_subjet_e    = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_E");
      m_ljet_subjet_charge = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_Charge");
      m_ljet_subjet_CSVv2  = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_CSVv2");
      m_ljet_subjet_CMVAv2 = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_CMVAv2");
      m_ljet_subjet_CvsL = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_CvsL");
      m_ljet_subjet_CvsB = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_CvsB");
      m_ljet_subjet_partonFlavour = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_PartonFlavour");
      m_ljet_subjet_hadronFlavour = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_HadronFlavour");
      m_ljet_subjet_neutralMultiplicity     = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_neutralMultiplicity");
      m_ljet_subjet_neutralHadronEnergyFrac = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_neutralHadronEnergyFrac");
      m_ljet_subjet_neutralEmEnergyFrac     = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_neutralEmEnergyFrac");
      m_ljet_subjet_chargedHadronEnergyFrac = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_chargedHadronEnergyFrac");
      m_ljet_subjet_chargedEmEnergyFrac = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_chargedEmEnergyFrac");
      m_ljet_subjet_chargedMultiplicity = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_chargedMultiplicity");
      m_ljet_subjet_jecFactor0 = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_jecFactor0");
      m_ljet_subjet_jetArea    = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_jetArea");
    }


    /** LEPTONS **/
    if (m_config->useLeptons()){
      m_el_size = new TTreeReaderValue<unsigned int>(m_ttree,"el_size");
      m_el_pt   = new TTreeReaderArray<float>(m_ttree,"el_Pt");
      m_el_eta  = new TTreeReaderArray<float>(m_ttree,"el_Eta");
      m_el_phi  = new TTreeReaderArray<float>(m_ttree,"el_Phi");
      m_el_e    = new TTreeReaderArray<float>(m_ttree,"el_E");
      m_el_charge  = new TTreeReaderArray<float>(m_ttree,"el_Charge");
      m_el_key     = new TTreeReaderArray<float>(m_ttree,"el_Key");
      m_el_iso03   = new TTreeReaderArray<float>(m_ttree,"el_Iso03");
      m_el_iso03db = new TTreeReaderArray<float>(m_ttree,"el_Iso03db");
      m_el_miniIso = new TTreeReaderArray<float>(m_ttree,"el_MiniIso");
      m_el_SCEta   = new TTreeReaderArray<float>(m_ttree,"el_SCEta");
      m_el_SCPhi   = new TTreeReaderArray<float>(m_ttree,"el_SCPhi");

      m_el_vidVeto        = new TTreeReaderArray<float>(m_ttree,"el_vidVeto");
      m_el_vidLoose       = new TTreeReaderArray<float>(m_ttree,"el_vidLoose");
      m_el_vidMedium      = new TTreeReaderArray<float>(m_ttree,"el_vidMedium");
      m_el_vidTight       = new TTreeReaderArray<float>(m_ttree,"el_vidTight");
      m_el_vidHEEP        = new TTreeReaderArray<float>(m_ttree,"el_vidHEEP");
      m_el_vidVetonoiso   = new TTreeReaderArray<float>(m_ttree,"el_vidVetonoiso");
      m_el_vidLoosenoiso  = new TTreeReaderArray<float>(m_ttree,"el_vidLoosenoiso");
      m_el_vidMediumnoiso = new TTreeReaderArray<float>(m_ttree,"el_vidMediumnoiso");
      m_el_vidTightnoiso  = new TTreeReaderArray<float>(m_ttree,"el_vidTightnoiso");
      m_el_vidHEEPnoiso   = new TTreeReaderArray<float>(m_ttree,"el_vidHEEPnoiso");
      m_el_vidMvaGPvalue  = new TTreeReaderArray<float>(m_ttree,"el_vidMvaGPvalue");
      m_el_vidMvaGPcateg  = new TTreeReaderArray<float>(m_ttree,"el_vidMvaGPcateg");
      m_el_vidMvaHZZvalue = new TTreeReaderArray<float>(m_ttree,"el_vidMvaHZZvalue");
      m_el_vidMvaHZZcateg = new TTreeReaderArray<float>(m_ttree,"el_vidMvaHZZcateg");

      m_el_veto_NoIsoID   = new TTreeReaderArray<int>(m_ttree,"el_IDVeto_NoIso");
      m_el_loose_NoIsoID  = new TTreeReaderArray<int>(m_ttree,"el_IDLoose_NoIso");
      m_el_medium_NoIsoID = new TTreeReaderArray<int>(m_ttree,"el_IDMedium_NoIso");
      m_el_tight_NoIsoID  = new TTreeReaderArray<int>(m_ttree,"el_IDTight_NoIso");
      m_el_isoVeto   = new TTreeReaderArray<int>(m_ttree,"el_IsoVeto");
      m_el_isoLoose  = new TTreeReaderArray<int>(m_ttree,"el_IsoLoose");
      m_el_isoMedium = new TTreeReaderArray<int>(m_ttree,"el_IsoMedium");
      m_el_isoTight  = new TTreeReaderArray<int>(m_ttree,"el_IsoTight");
      m_el_vetoID    = new TTreeReaderArray<int>(m_ttree,"el_IDVeto");
      m_el_looseID   = new TTreeReaderArray<int>(m_ttree,"el_IDLoose");
      m_el_mediumID  = new TTreeReaderArray<int>(m_ttree,"el_IDMedium");
      m_el_tightID   = new TTreeReaderArray<int>(m_ttree,"el_IDTight");

      m_mu_size = new TTreeReaderValue<int>(m_ttree,"mu_size");
      m_mu_pt   = new TTreeReaderArray<float>(m_ttree,"mu_Pt");
      m_mu_eta  = new TTreeReaderArray<float>(m_ttree,"mu_Eta");
      m_mu_phi  = new TTreeReaderArray<float>(m_ttree,"mu_Phi");
      m_mu_e    = new TTreeReaderArray<float>(m_ttree,"mu_E");
      m_mu_charge  = new TTreeReaderArray<float>(m_ttree,"mu_Charge");
      m_mu_key     = new TTreeReaderArray<float>(m_ttree,"mu_Key");
      m_mu_iso04   = new TTreeReaderArray<float>(m_ttree,"mu_Iso04");
      m_mu_miniIso = new TTreeReaderArray<float>(m_ttree,"mu_MiniIso");
      m_mu_soft       = new TTreeReaderArray<float>(m_ttree,"mu_IsSoftMuon");
      m_mu_loose      = new TTreeReaderArray<float>(m_ttree,"mu_IsLooseMuon");
      m_mu_medium     = new TTreeReaderArray<float>(m_ttree,"mu_IsMediumMuon");
      m_mu_medium2016 = new TTreeReaderArray<float>(m_ttree,"mu_IsMediumMuon2016");
      m_mu_tight   = new TTreeReaderArray<float>(m_ttree,"mu_IsTightMuon");
      m_mu_hightPt = new TTreeReaderArray<float>(m_ttree,"mu_IsHighPtMuon");
    }

    if (!m_kinematicReco && m_config->useNeutrinos()){
      if (m_config->useNeutrinos()){
        // Neutrinos aren't stored in the baseline ntuples, requires 'kinematicReco' to create
        m_nu_pt  = new TTreeReaderArray<float>(m_ttree, "nu_pt");
        m_nu_eta = new TTreeReaderArray<float>(m_ttree, "nu_eta");
        m_nu_phi = new TTreeReaderArray<float>(m_ttree, "nu_phi");
      }

      if (m_config->useTtbar()){
        // ttbar system isn't stored in the baseline ntuples, requires 'kinematicReco' to create
        m_top_pt  = new TTreeReaderValue<float>(m_ttree, "top_pt");
        m_top_eta = new TTreeReaderValue<float>(m_ttree, "top_eta");
        m_top_phi = new TTreeReaderValue<float>(m_ttree, "top_phi");
        m_top_e   = new TTreeReaderValue<float>(m_ttree, "top_e");
        m_lepton_top_index = new TTreeReaderValue<int>(m_ttree, "lepton_top_index");
        m_jet_top_index    = new TTreeReaderValue<int>(m_ttree, "jet_top_index");
        m_nu_top_index     = new TTreeReaderValue<int>(m_ttree, "nu_top_index");

        m_antitop_pt  = new TTreeReaderValue<float>(m_ttree, "antitop_pt");
        m_antitop_eta = new TTreeReaderValue<float>(m_ttree, "antitop_eta");
        m_antitop_phi = new TTreeReaderValue<float>(m_ttree, "antitop_phi");
        m_antitop_e   = new TTreeReaderValue<float>(m_ttree, "antitop_e");
        m_lepton_antitop_index = new TTreeReaderValue<int>(m_ttree, "lepton_antitop_index");
        m_jet_antitop_index    = new TTreeReaderValue<int>(m_ttree, "jet_antitop_index");
        m_nu_antitop_index     = new TTreeReaderValue<int>(m_ttree, "nu_antitop_index");

        m_dileptonTtbarWeight = new TTreeReaderValue<float>(m_ttree, "dileptonTtbar");
      }
    }


    m_met_size = new TTreeReaderValue<unsigned int>(m_ttree,"m_met_size");
    m_met_met  = new TTreeReaderArray<float>(m_ttree,"met_Pt");
    m_met_phi  = new TTreeReaderArray<float>(m_ttree,"met_Phi");
    m_met_met_uncor = new TTreeReaderArray<float>(m_ttree,"met_uncorPt");
    m_met_phi_uncor = new TTreeReaderArray<float>(m_ttree,"met_uncorPhi");

    m_met_muCleanOnly_size = new TTreeReaderValue<unsigned int>(m_ttree,"met_muCleanOnly_size");
    m_met_muCleanOnly_met  = new TTreeReaderArray<float>(m_ttree,"met_MuCleanOnly_Pt");
    m_met_muCleanOnly_phi  = new TTreeReaderArray<float>(m_ttree,"met_MuCleanOnly_Phi");
    m_met_muCleanOnly_met_uncor = new TTreeReaderArray<float>(m_ttree,"met_MuCleanOnly_uncorPt");
    m_met_muCleanOnly_phi_uncor = new TTreeReaderArray<float>(m_ttree,"met_MuCleanOnly_uncorPhi");

    m_met_syst_size = new TTreeReaderValue<unsigned int>(m_ttree,"met_syst_size");
    m_met_syst_met  = new TTreeReaderArray<float>(m_ttree,"metsyst_Pt");
    m_met_syst_phi  = new TTreeReaderArray<float>(m_ttree,"metsyst_Phi");
    m_met_syst_muCleanOnly_met = new TTreeReaderArray<float>(m_ttree,"metsyst_MuCleanOnly_Pt");
    m_met_syst_muCleanOnly_phi = new TTreeReaderArray<float>(m_ttree,"metsyst_MuCleanOnly_Phi");


    // set some event weights and access necessary branches
    m_xsection       = 1.0;
    m_kfactor        = 1.0;
    m_sumOfWeights   = 1.0;
    m_LUMI           = m_config->LUMI();

    // MC information
    if (m_isMC){
      m_weight_mc   = new TTreeReaderValue<float>(m_ttree,"evt_Gen_Weight");
      m_xsection     = m_config->XSectionMap( m_fileName );
      m_kfactor      = m_config->KFactorMap(  m_fileName );
      m_sumOfWeights = m_config->sumWeightsMap( m_fileName );

      m_mc_ht = new TTreeReaderValue<float>(m_ttree,"evt_Gen_Ht");

      m_MC_part1_factor = new TTreeReaderArray<float>(m_ttree,"MC_part1_factor");
      m_MC_part1_ID = new TTreeReaderArray<float>(m_ttree,"MC_part1_ID");
      m_MC_part2_factor = new TTreeReaderArray<float>(m_ttree,"MC_part2_factor");
      m_MC_part2_ID = new TTreeReaderArray<float>(m_ttree,"MC_part2_ID");

      m_MC_t_pt  = new TTreeReaderArray<float>(m_ttree,"MC_t_pt");
      m_MC_t_eta = new TTreeReaderArray<float>(m_ttree,"MC_t_eta");
      m_MC_t_phi = new TTreeReaderArray<float>(m_ttree,"MC_t_phi");
      m_MC_t_e   = new TTreeReaderArray<float>(m_ttree,"MC_t_E");
      m_MC_tbar_pt  = new TTreeReaderArray<float>(m_ttree,"MC_tbar_pt");
      m_MC_tbar_eta = new TTreeReaderArray<float>(m_ttree,"MC_tbar_eta");
      m_MC_tbar_phi = new TTreeReaderArray<float>(m_ttree,"MC_tbar_phi");
      m_MC_tbar_e = new TTreeReaderArray<float>(m_ttree,"MC_tbar_E");
      m_MC_lep_pt = new TTreeReaderArray<float>(m_ttree,"MC_lep_pt");
      m_MC_lep_eta = new TTreeReaderArray<float>(m_ttree,"MC_lep_eta");
      m_MC_lep_phi = new TTreeReaderArray<float>(m_ttree,"MC_lep_phi");
      m_MC_lep_e  = new TTreeReaderArray<float>(m_ttree,"MC_lep_E");
      m_MC_lep_ID = new TTreeReaderArray<float>(m_ttree,"MC_lep_ID");
      m_MC_nu_pt  = new TTreeReaderArray<float>(m_ttree,"MC_nu_pt");
      m_MC_nu_eta = new TTreeReaderArray<float>(m_ttree,"MC_nu_eta");
      m_MC_nu_phi = new TTreeReaderArray<float>(m_ttree,"MC_nu_phi");
      m_MC_nu_e   = new TTreeReaderArray<float>(m_ttree,"MC_nu_E");
      m_MC_lepb_pt  = new TTreeReaderArray<float>(m_ttree,"MC_lepb_pt");
      m_MC_lepb_eta = new TTreeReaderArray<float>(m_ttree,"MC_lepb_eta");
      m_MC_lepb_phi = new TTreeReaderArray<float>(m_ttree,"MC_lepb_phi");
      m_MC_lepb_e   = new TTreeReaderArray<float>(m_ttree,"MC_lepb_E");
      m_MC_hadW_pt  = new TTreeReaderArray<float>(m_ttree,"MC_hadW_pt");
      m_MC_hadW_eta = new TTreeReaderArray<float>(m_ttree,"MC_hadW_eta");
      m_MC_hadW_phi = new TTreeReaderArray<float>(m_ttree,"MC_hadW_phi");
      m_MC_hadW_e   = new TTreeReaderArray<float>(m_ttree,"MC_hadW_E");
      m_MC_hadb_pt  = new TTreeReaderArray<float>(m_ttree,"MC_hadb_pt");
      m_MC_hadb_eta = new TTreeReaderArray<float>(m_ttree,"MC_hadb_eta");
      m_MC_hadb_phi = new TTreeReaderArray<float>(m_ttree,"MC_hadb_phi");
      m_MC_hadb_e = new TTreeReaderArray<float>(m_ttree,"MC_hadb_E");
      m_MC_cstar  = new TTreeReaderArray<float>(m_ttree,"MC_cstar");
      m_MC_x_F    = new TTreeReaderArray<float>(m_ttree,"MC_x_F");
      m_MC_Mtt    = new TTreeReaderArray<float>(m_ttree,"MC_Mtt");

      m_truth_jet_pt  = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_GenJetPt");
      m_truth_jet_eta = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_GenJetEta");
      m_truth_jet_phi = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_GenJetPhi");
      m_truth_jet_e   = new TTreeReaderArray<float>(m_ttree,"jetAK4CHS_GenJetCharge");

      m_truth_ljet_pt  = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_GenJetPt");
      m_truth_ljet_eta = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_GenJetEta");
      m_truth_ljet_phi = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_GenJetPhi");
      m_truth_ljet_e   = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_GenJetE");
      m_truth_ljet_charge     = new TTreeReaderArray<float>(m_ttree,"jetAK8CHS_GenJetCharge");
      m_truth_ljet_subjet_pt  = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_GenJetPt");
      m_truth_ljet_subjet_eta = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_GenJetEta");
      m_truth_ljet_subjet_phi = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_GenJetPhi");
      m_truth_ljet_subjet_e   = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_GenJetE");
      m_truth_ljet_subjet_charge = new TTreeReaderArray<float>(m_ttree,"subjetAK8CHS_GenJetCharge");

      if (m_config->useTruth()){
        //TTree* truth_tree = (TTree*)m_ttree.GetTree()->GetCurrentFile()->Get("truth");
        //m_truth_tree.SetTree(truth_tree);
        //m_truthEventNumber = new TTreeReaderValue<unsigned long long>(m_truth_tree,  "eventNumber");
        //m_truthRunNumber   = new TTreeReaderValue<unsigned int>(m_truth_tree,"runNumber");
        //m_truth_weight_mc  = new TTreeReaderValue<float>(m_truth_tree,"weight_mc");

        m_truth_ljet_pt  = new TTreeReaderArray<float>(m_ttree,"truth_ljet_pt");
        m_truth_ljet_eta = new TTreeReaderArray<float>(m_ttree,"truth_ljet_eta");
        m_truth_ljet_phi = new TTreeReaderArray<float>(m_ttree,"truth_ljet_phi");
        m_truth_ljet_e   = new TTreeReaderArray<float>(m_ttree,"truth_ljet_e");
        m_truth_ljet_Qw  = new TTreeReaderArray<float>(m_ttree,"truth_ljet_Qw");
        m_truth_ljet_tau32_wta = new TTreeReaderArray<float>(m_ttree,"truth_ljet_tau32_wta");
        m_truth_ljet_split23   = new TTreeReaderArray<float>(m_ttree,"truth_ljet_split23");
      } // end useTruth
    } // end isMC

    // DNN material
    bool useDNN(false);
    if (!m_getDNN && useDNN)  // always false for now
        m_dnn_score = new TTreeReaderValue<float>(m_ttree,"dnn_score");

    std::ifstream input_cfg = cma::open_file( m_config->dnnFile() );
    lwt::JSONConfig cfg     = lwt::parse_json( input_cfg );
    m_lwnn   = new lwt::LightweightNeuralNetwork(cfg.inputs, cfg.layers, cfg.outputs);
    m_dnnKey = m_config->dnnKey();

    // Kinematic reconstruction algorithms
    // m_semileptonTtbar = new semileptonTtbarReco(m_config);  // semi-leptonic ttbar kinematic reco
    // m_allhadTtbar = new allhadTtbarReco(m_config);          // all-hadronic ttbar kinematic reco
    m_dileptonTtbar = new dileptonTtbarReco(cmaConfig, configuration::run2_13tev_2016_25ns, 2, true);
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
        m_weightSystematicsVectorFloats[syst.first] = new TTreeReaderArray<float>(m_ttree,syst.first.c_str());

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


bool Event::isValidRecoEntry(){
    return (m_entry > (long long)-1);
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

    // Get the event weights (for cutflow & histograms)
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

    // Get some kinematic variables (MET, HT, ST)
    initialize_kinematics();
    cma::DEBUG("EVENT : Setup kinematic variables ");



    // kinematic reconstruction -- usually needed for neutrinos!
    if (m_kinematicReco){
        // if 0/1-lepton:
        if (m_isTwoLeptonAnalysis){
            getDilepton();
            cma::DEBUG("EVENT : Setup the dilepton struct. ");
        }
    }
    // build the ttbar system (depends on analysis!)
    // use the reconstruction, or load from root file; depends on m_kinematicReco value
    buildTtbar();
    cma::DEBUG("EVENT : Ttbar system constructed & defined");

    // Neutrinos
    if (m_config->useNeutrinos()){
        // relies on kinematic reconstruction, unless the information is saved in root file
        initialize_neutrinos();
        cma::DEBUG("EVENT : Setup neutrinos ");
    }


    // DNN
    if (m_getDNN){
        cma::DEBUG("EVENT : Calculate DNN ");
        getDNN();
    }
    else{
        // load from ntuple
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
    m_jets.resize( **m_jet_size );   // (*m_jet_pt)->size());
    for (const auto& btagWP : m_config->btagWkpts() ){
        m_btag_jets[btagWP].clear();
    }

    for (unsigned int i=0,size=**m_jet_size; i<size; i++){
        Jet jet;
        jet.p4.SetPtEtaPhiE(m_jet_pt->At(i),m_jet_eta->At(i),m_jet_phi->At(i),m_jet_e->At(i));

        jet.cMVAv2 = m_jet_CMVAv2->At(i);
        jet.index = i;

        getBtaggedJets(jet);

        m_jets[i] = jet;
    }

    m_btag_jets_default = m_btag_jets.at(m_config->jet_btagWkpt());

    return;
}


void Event::initialize_ljets(){
    /* Setup struct of large-R jets and relevant information */
    m_ljets.resize(**m_ljet_size);

    for (unsigned int i=0,size=**m_ljet_size; i<size; i++){
        // pre-selection (should be done in AnalysisTop?)
        // leading pT>500 GeV; others pT>350 GeV

        Ljet ljet;
        ljet.p4.SetPtEtaPhiE( m_ljet_pt->At(i),m_ljet_eta->At(i),m_ljet_phi->At(i),m_ljet_e->At(i));
        ljet.charge    = m_ljet_charge->At(i);

        ljet.tau1_CHS  = m_ljet_tau1_CHS->At(i);
        ljet.tau2_CHS  = m_ljet_tau2_CHS->At(i);
        ljet.tau3_CHS  = m_ljet_tau3_CHS->At(i);
        ljet.tau21_CHS = ljet.tau2_CHS / ljet.tau1_CHS;
        ljet.tau32_CHS = ljet.tau3_CHS / ljet.tau2_CHS;

        ljet.softDropMass_CHS = m_ljet_softDropMass_CHS->At(i);
        ljet.vSubjetIndex0    = m_ljet_vSubjetIndex0->At(i);
        ljet.vSubjetIndex1    = m_ljet_vSubjetIndex1->At(i);

        ljet.isGood    = (ljet.p4.Pt()>200000. && fabs(ljet.p4.Eta())<2.0) ? 1 : 0;

        m_ljets[i] = ljet;
    }

    return;
}


void Event::initialize_leptons(){
    /* Setup struct of lepton and relevant information */
    m_ee   = false;
    m_mumu = false;
    m_emu  = false;

    m_leptons.clear();

    Lepton lep;
    m_leptons.push_back(lep);

    return;
}


void Event::initialize_neutrinos(){
    /* Build the neutrinos */
    m_neutrinos.clear();

    Neutrino nu1;
    Neutrino nu2;

    // Assign neutrinos from 'm_ttbar'
    if (m_isOneLeptonAnalysis){
        nu1 = m_ttbar["top"].neutrino;
        m_neutrinos.push_back( nu1 );
    }
    else if (m_isTwoLeptonAnalysis){
        nu1 = m_ttbar["top"].neutrino;
        nu2 = m_ttbar["antitop"].neutrino;

        // pT-ordering
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
        cma::ERROR("EVENT : 'initialize_neutrinos()' called for analysis without neutrinos!");
    }

    return;
}


void Event::initialize_weights(){
    /* Event weights */
    m_nominal_weight = 1.0;

    m_weight_btag.clear();
    if (m_isMC){
        m_nominal_weight  = (**m_weight_pileup) * (**m_weight_mc);
        m_nominal_weight *= (m_xsection) * (m_kfactor) * m_LUMI / (m_sumOfWeights);
/*      // weights not in CMS (so far):
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
    m_metmet = m_met_met->At(0);
    m_metphi = m_met_phi->At(0);

    // Get MET and lepton transverse energy
    m_ST += m_HT;
    m_ST += m_metmet;
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
    dilep_met.SetX( m_metmet*cos(m_metphi) );
    dilep_met.SetY( m_metmet*sin(m_metphi) );
    m_dilepton.met = dilep_met;

    // Jets
    std::vector<Jet> bjets;
    for (const auto& j : m_btag_jets_default){
        bjets.push_back(m_jets.at(j));
    }
    m_dilepton.jets  = m_jets;
    m_dilepton.bjets = bjets;

    cma::DEBUG("EVENT : Dilepton");
    cma::DEBUG("EVENT : met       = "+std::to_string(m_metmet));
    cma::DEBUG("EVENT : lepton pT = "+std::to_string(m_leptons.at(0).p4.Pt()));
    cma::DEBUG("EVENT : jet pT    = "+std::to_string(m_jets.at(0).p4.Pt()));

    return;
}



//BUILD TTBAR
//    USE ROOT FILE INFORMATION IF "m_kinematicReco"==false
//    ELSE USE ALGORITHM (DEPENDS ON 0/1/2-lepton ANALYSIS

void Event::buildTtbar(){
    /* Build ttbar system */
    if (m_kinematicReco){
        getDilepton();
        m_ttbar = m_dileptonTtbar->execute(m_dilepton);
    }
    else{
        // no kinematic reconstruction performed.
        // re-define the ttbar system (map named 'm_ttbar')
        int nu_t = **m_nu_top_index;
        Neutrino nu_top; // neutrino
        nu_top.p4.SetPtEtaPhiM( m_nu_pt->At(nu_t), m_nu_eta->At(nu_t), m_nu_phi->At(nu_t), 0.0 );

        int nu_atop = **m_nu_antitop_index;
        Neutrino nu_antitop; // anti-neutrino
        nu_antitop.p4.SetPtEtaPhiM( m_nu_pt->At(nu_atop), m_nu_eta->At(nu_atop), m_nu_phi->At(nu_atop), 0.0 );

        // fill m_ttbar with values from TTree (top, lepton, jet, neutrino)
        LepTop top;
        top.p4.SetPtEtaPhiE(*(*m_top_pt), *(*m_top_eta), *(*m_top_phi), *(*m_top_e));
        top.lepton = m_leptons.at( *(*m_lepton_top_index) );
        top.jet    = m_jets.at( *(*m_jet_top_index) );
        top.neutrino = nu_top;
        top.weight = 0.0;
        top.weight_ES = 0.0;
        top.weight_tt = *(*m_dileptonTtbarWeight);

        LepTop antitop;
        antitop.p4.SetPtEtaPhiE(*(*m_antitop_pt), *(*m_antitop_eta), *(*m_antitop_phi), *(*m_antitop_e));
        antitop.lepton = m_leptons.at( *(*m_lepton_antitop_index) );
        antitop.jet    = m_jets.at( *(*m_jet_antitop_index) );
        antitop.neutrino = nu_antitop;
        antitop.weight = 0.0;
        antitop.weight_ES = 0.0;
        antitop.weight_tt = *(*m_dileptonTtbarWeight);

        m_ttbar = {};
        m_ttbar["top"]     = top;
        m_ttbar["antitop"] = antitop;
    }

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


void Event::getBtaggedJets( Jet& jet ){
    /* Determine the b-tagging */
    jet.isbtagged["L"] = false;
    jet.isbtagged["M"] = false;
    jet.isbtagged["T"] = false;

    if (jet.cMVAv2 > m_cMVAv2L){
        jet.isbtagged["L"] = true;
        m_btag_jets["L"].push_back(jet.index);  // 0 = index of this jet
        if (jet.cMVAv2 > m_cMVAv2M){
            jet.isbtagged["M"] = true;
            m_btag_jets["M"].push_back(jet.index);
            if (jet.cMVAv2 > m_cMVAv2T){
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
    else if (syst.find("jvt")!=std::string::npos){
        // pileup event weight
        syst_event_weight  = (**m_weight_pileup) * (**m_weight_mc);
        syst_event_weight *= m_weight_btag_default;
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);

        syst_event_weight *= **m_weightSystematicsFloats.at(syst);
    }
    else if (syst.find("pileup")!=std::string::npos){
        // pileup event weight
        syst_event_weight  = (**m_weight_mc);
        syst_event_weight *= m_weight_btag_default;
        syst_event_weight *= (m_xsection) * (m_kfactor) * (m_LUMI);
        syst_event_weight /= (m_sumOfWeights);

        syst_event_weight *= **m_weightSystematicsFloats.at(syst);
    }
    else if (syst.find("leptonSF")!=std::string::npos){
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

    return syst_event_weight;
}




/*** RETURN PHYSICS INFORMATION ***/
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

std::map<std::string,LepTop> Event::ttbar(){
    return m_ttbar;
}

double Event::DNN(){
    /* Return the DNN value */
    return m_DNN;
}


/*** RETURN WEIGHTS ***/
float Event::nominal_weight(){
    return m_nominal_weight;
}
float Event::truth_weight_mc(){
    return **m_truth_weight_mc;
}
float Event::weight_mc(){
    return **m_weight_mc;
}
float Event::weight_pileup(){
    return **m_weight_pileup;
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

float Event::mu(){
    return **m_mu;
}

int Event::lumiblock(){
    return **m_lumiblock;
}




/*** DELETE VARIABLES ***/
void Event::finalize(){
    // delete variables
    cma::DEBUG("EVENT : Finalize() ");
    delete m_lwnn;
    delete m_dileptonTtbar;
    delete m_eventNumber;
    delete m_runNumber;
//    delete m_mcChannelNumber;
//    delete m_mu;
//    delete m_lumiblock;
    if (m_config->useLargeRJets()){
      delete m_ljet_size;        // jetAK8CHS_size
      delete m_ljet_pt;          // jetAK8CHS_Pt
      delete m_ljet_eta;         // jetAK8CHS_Eta
      delete m_ljet_phi;         // jetAK8CHS_Phi
      delete m_ljet_e;           // jetAK8CHS_E
      delete m_ljet_tau1_CHS;    // jetAK8CHS_tau1CHS
      delete m_ljet_tau2_CHS;    // jetAK8CHS_tau1CHS
      delete m_ljet_tau3_CHS;    // jetAK8CHS_tau1CHS
      delete m_ljet_charge;      // jetAK8CHS_Charge
      delete m_ljet_softDropMass_CHS;  // jetAK8CHS_softDropMassCHS
      delete m_ljet_CSVv2;                   // jetAK8CHS_CSVv2
      delete m_ljet_CMVAv2;                  // jetAK8CHS_CMVAv2
      delete m_ljet_CvsL;                    // jetAK8CHS_CvsL
      delete m_ljet_CvsB;                    // jetAK8CHS_CvsB
      delete m_ljet_PartonFlavour;           // jetAK8CHS_PartonFlavour
      delete m_ljet_HadronFlavour;           // jetAK8CHS_HadronFlavour
      delete m_ljet_neutralMultiplicity;     // jetAK8CHS_neutralMultiplicity
      delete m_ljet_neutralHadronEnergyFrac; // jetAK8CHS_neutralHadronEnergyFrac
      delete m_ljet_neutralEmEnergyFrac;     // jetAK8CHS_neutralEmEnergyFrac
      delete m_ljet_chargedHadronEnergyFrac; // jetAK8CHS_chargedHadronEnergyFrac
      delete m_ljet_chargedEmEnergyFrac;     // jetAK8CHS_chargedEmEnergyFrac
      delete m_ljet_chargedMultiplicity;     // jetAK8CHS_chargedMultiplicity
      delete m_ljet_jecFactor0;              // jetAK8CHS_jecFactor0
      delete m_ljet_jetArea;                 // jetAK8CHS_jetArea
      delete m_ljet_jecUncertainty;          // jetAK8CHS_jecUncertainty
      delete m_ljet_PtResolution;            // jetAK8CHS_PtResolution
      delete m_ljet_JERSF;                   // jetAK8CHS_JERSF
      delete m_ljet_JERSFUp;                 // jetAK8CHS_JERSFUp
      delete m_ljet_JERSFDown;               // jetAK8CHS_JERSFDown
      delete m_ljet_SmearedPt;               // jetAK8CHS_SmearedPt
      delete m_ljet_vSubjetIndex0;           // jetAK8CHS_vSubjetIndex0
      delete m_ljet_vSubjetIndex1;           // jetAK8CHS_vSubjetIndex1
      delete m_ljet_keys;                    // jetAK8CHS_Keys
      delete m_ljet_subjet_size;                    // subjetAK8CHS_size
      delete m_ljet_subjet_pt;                      // subjetAK8CHS_Pt
      delete m_ljet_subjet_eta;                     // subjetAK8CHS_Eta
      delete m_ljet_subjet_phi;                     // subjetAK8CHS_Phi
      delete m_ljet_subjet_e;                       // subjetAK8CHS_E
      delete m_ljet_subjet_charge;                  // subjetAK8CHS_Charge
      delete m_ljet_subjet_CSVv2;                   // subjetAK8CHS_CSVv2
      delete m_ljet_subjet_CMVAv2;                  // subjetAK8CHS_CMVAv2
      delete m_ljet_subjet_CvsL;                    // subjetAK8CHS_CvsL
      delete m_ljet_subjet_CvsB;                    // subjetAK8CHS_CvsB
      delete m_ljet_subjet_partonFlavour;           // subjetAK8CHS_PartonFlavour
      delete m_ljet_subjet_hadronFlavour;           // subjetAK8CHS_HadronFlavour
      delete m_ljet_subjet_neutralMultiplicity;     // subjetAK8CHS_neutralMultiplicity
      delete m_ljet_subjet_neutralHadronEnergyFrac; // subjetAK8CHS_neutralHadronEnergyFrac
      delete m_ljet_subjet_neutralEmEnergyFrac;     // subjetAK8CHS_neutralEmEnergyFrac
      delete m_ljet_subjet_chargedHadronEnergyFrac; // subjetAK8CHS_chargedHadronEnergyFrac
      delete m_ljet_subjet_chargedEmEnergyFrac;     // subjetAK8CHS_chargedEmEnergyFrac
      delete m_ljet_subjet_chargedMultiplicity;     // subjetAK8CHS_chargedMultiplicity
      delete m_ljet_subjet_jecFactor0;              // subjetAK8CHS_jecFactor0
      delete m_ljet_subjet_jetArea;                 // subjetAK8CHS_jetArea
      delete m_ljet_subjet_keys;                    // subjetAK8CHS_Keys
    }
    if (m_config->useLeptons()){
      delete m_el_size;    // el_size
      delete m_el_pt;             // el_Pt
      delete m_el_eta;            // el_Eta
      delete m_el_phi;            // el_Phi
      delete m_el_e;              // el_E
      delete m_el_charge;         // el_Charge
      delete m_el_key;            // el_Key
      delete m_el_iso03;          // el_Iso03
      delete m_el_iso03db;        // el_Iso03db
      delete m_el_miniIso;        // el_MiniIso
      delete m_el_SCEta;          // el_SCEta
      delete m_el_SCPhi;          // el_SCPhi
      delete m_el_vidVeto;        // el_vidVeto
      delete m_el_vidLoose;       // el_vidLoose
      delete m_el_vidMedium;      // el_vidMedium
      delete m_el_vidTight;       // el_vidTight
      delete m_el_vidHEEP;        // el_vidHEEP
      delete m_el_vidVetonoiso;   // el_vidVetonoiso
      delete m_el_vidLoosenoiso;  // el_vidLoosenoiso
      delete m_el_vidMediumnoiso; // el_vidMediumnoiso
      delete m_el_vidTightnoiso;  // el_vidTightnoiso
      delete m_el_vidHEEPnoiso;   // el_vidHEEPnoiso
      delete m_el_vidMvaGPvalue;  // el_vidMvaGPvalue
      delete m_el_vidMvaGPcateg;  // el_vidMvaGPcateg
      delete m_el_vidMvaHZZvalue; // el_vidMvaHZZvalue
      delete m_el_vidMvaHZZcateg; // el_vidMvaHZZcateg
      delete m_el_veto_NoIsoID;     // el_IDVeto_NoIso
      delete m_el_loose_NoIsoID;    // el_IDLoose_NoIso
      delete m_el_medium_NoIsoID;   // el_IDMedium_NoIso
      delete m_el_tight_NoIsoID;    // el_IDTight_NoIso
      delete m_el_isoVeto;          // el_IsoVeto
      delete m_el_isoLoose;         // el_IsoLoose
      delete m_el_isoMedium;        // el_IsoMedium
      delete m_el_isoTight;         // el_IsoTight
      delete m_el_vetoID;           // el_IDVeto
      delete m_el_looseID;          // el_IDLoose
      delete m_el_mediumID;         // el_IDMedium
      delete m_el_tightID;          // el_IDTight


      delete m_mu_size;         // mu_size
      delete m_mu_pt;         // mu_Pt
      delete m_mu_eta;        // mu_Eta
      delete m_mu_phi;        // mu_Phi
      delete m_mu_e;          // mu_E
      delete m_mu_charge;     // mu_Charge
      delete m_mu_key;        // mu_Key
      delete m_mu_iso04;      // mu_Iso04
      delete m_mu_miniIso;    // mu_MiniIso
      delete m_mu_soft;       // mu_IsSoftMuon
      delete m_mu_loose;      // mu_IsLooseMuon
      delete m_mu_medium;     // mu_IsMediumMuon
      delete m_mu_medium2016; // mu_IsMediumMuon2016
      delete m_mu_tight;      // mu_IsTightMuon
      delete m_mu_hightPt;    // mu_IsHighPtMuon
    }
    if (m_config->useJets()){

    }
    delete m_met_size;  // met_size
    delete m_met_met;          // met_Pt
    delete m_met_phi;          // met_Phi
    delete m_met_met_uncor;    // met_uncorPt
    delete m_met_phi_uncor;    // met_uncorPhi

    delete m_met_muCleanOnly_size;  // met_MuCleanOnly_size
    delete m_met_muCleanOnly_met;          // met_MuCleanOnly_Pt
    delete m_met_muCleanOnly_phi;          // met_MuCleanOnly_Phi
    delete m_met_muCleanOnly_met_uncor;    // met_MuCleanOnly_uncorPt
    delete m_met_muCleanOnly_phi_uncor;    // met_MuCleanOnly_uncorPhi

    delete m_met_syst_size;  // metsyst_size
    delete m_met_syst_met;          // metsyst_Pt
    delete m_met_syst_phi;          // metsyst_Phi
    delete m_met_syst_muCleanOnly_met;    // metsyst_MuCleanOnly_Pt
    delete m_met_syst_muCleanOnly_phi;    // metsyst_MuCleanOnly_Phi

    if (m_config->useFlags()){
      delete m_Flag_BadPFMuonFilter;                    // Flag_BadPFMuonFilter
      delete m_Flag_BadChargedCandidateFilter;          // Flag_BadChargedCandidateFilter
      delete m_Flag_HBHENoiseFilter;                    // Flag_HBHENoiseFilter
      delete m_Flag_HBHENoiseIsoFilter;                 // Flag_HBHENoiseIsoFilter
      delete m_Flag_CSCTightHaloFilter;                 // Flag_CSCTightHaloFilter
      delete m_Flag_CSCTightHaloTrkMuUnvetoFilter;      // Flag_CSCTightHaloTrkMuUnvetoFilter
      delete m_Flag_CSCTightHalo2015Filter;             // Flag_CSCTightHalo2015Filter
      delete m_Flag_globalTightHalo2016Filter;          // Flag_globalTightHalo2016Filter
      delete m_Flag_globalSuperTightHalo2016Filter;     // Flag_globalSuperTightHalo2016Filter
      delete m_Flag_HcalStripHaloFilter;                // Flag_HcalStripHaloFilter
      delete m_Flag_hcalLaserEventFilter;               // Flag_hcalLaserEventFilter
      delete m_Flag_EcalDeadCellTriggerPrimitiveFilter; // Flag_EcalDeadCellTriggerPrimitiveFilter
      delete m_Flag_EcalDeadCellBoundaryEnergyFilter;   // Flag_EcalDeadCellBoundaryEnergyFilter
      delete m_Flag_goodVertices;                       // Flag_goodVertices
      delete m_Flag_eeBadScFilter;                      // Flag_eeBadScFilter
      delete m_Flag_ecalLaserCorrFilter;                // Flag_ecalLaserCorrFilter
      delete m_Flag_trkPOGFilters;                      // Flag_trkPOGFilters
      delete m_Flag_chargedHadronTrackResolutionFilter; // Flag_chargedHadronTrackResolutionFilter
      delete m_Flag_muonBadTrackFilter;                 // Flag_muonBadTrackFilter
      delete m_Flag_trkPOG_manystripclus53X;            // Flag_trkPOG_manystripclus53X
      delete m_Flag_trkPOG_toomanystripclus53X;         // Flag_trkPOG_toomanystripclus53X
      delete m_Flag_trkPOG_logErrorTooManyClusters;     // Flag_trkPOG_logErrorTooManyClusters
      delete m_Flag_METFilters;                         // Flag_METFilters
      delete m_Flag_badMuons;                           // Flag_badMuons
      delete m_Flag_duplicateMuons;                     // Flag_duplicateMuons
      delete m_Flag_noBadMuons;                         // Flag_noBadMuons
    }

    // Event info 
    delete m_eventNumber;          // evt_EventNumber
    delete m_runNumber;            // evt_RunNumber
    delete m_mu;
    delete m_rho;               // evt_rho
    delete m_lumiblock;            // evt_LumiBlock
    delete m_treeXSection;       // evt_XSec
    delete m_treeKFactor;
    delete m_treeSumOfWeights;
    delete m_NGoodVtx;             // evt_NGoodVtx
    delete m_LHAPDF_ID;            // evt_LHA_PDF_ID
    delete m_NIsoTrk;              // evt_NIsoTrk
    delete m_pu_NtrueInt;          // pu_NtrueInt

    delete m_HLT_Ele45_WPLoose_Gsf;          // HLT_Ele45_WPLoose_Gsf
    delete m_HLT_Ele45_WPLoose_Gsf_prescale; // HLT_Ele45_WPLoose_Gsf_prescale
    delete m_HLT_Mu50;            // HLT_Mu50
    delete m_HLT_Mu50_prescale;   // HLT_Mu50_prescale
    delete m_HLT_TkMu50;          // HLT_TkMu50
    delete m_HLT_TkMu50_prescale; // HLT_TkMu50_prescale

    delete m_scale_size;         // scale_size
    delete m_scale_Weights;      // scale_Weights
    delete m_pdf_size;           // pdf_size
    delete m_weights_pdf;        // pdf_Weights
    delete m_alphas_size;        // alphas_size
    delete m_weights_alphas;     // alphas_Weights

    if (m_isMC){
      delete m_weight_mc;
      delete m_weight_pileup;
      delete m_weight_pileup_UP;
      delete m_weight_pileup_DOWN;

      if (m_config->useTruth()){
        delete m_truthEventNumber;
        delete m_truthRunNumber;
        delete m_truth_weight_mc;
        delete m_mc_ht;

        delete m_truth_jet_pt;           // jetAK4CHS_GenJetPt
        delete m_truth_jet_eta;          // jetAK4CHS_GenJetEta
        delete m_truth_jet_phi;          // jetAK4CHS_GenJetPhi
        delete m_truth_jet_e;            // jetAK4CHS_GenJetCharge

        delete m_truth_ljet_pt;            // jetAK8CHS_GenJetPt
        delete m_truth_ljet_eta;           // jetAK8CHS_GenJetEta
        delete m_truth_ljet_phi;           // jetAK8CHS_GenJetPhi
        delete m_truth_ljet_e;             // jetAK8CHS_GenJetE
        delete m_truth_ljet_charge;        // jetAK8CHS_GenJetCharge
        delete m_truth_ljet_subjet_pt;     // subjetAK8CHS_GenJetPt
        delete m_truth_ljet_subjet_eta;    // subjetAK8CHS_GenJetEta
        delete m_truth_ljet_subjet_phi;    // subjetAK8CHS_GenJetPhi
        delete m_truth_ljet_subjet_e;      // subjetAK8CHS_GenJetE
        delete m_truth_ljet_subjet_charge; // subjetAK8CHS_GenJetCharge

        delete m_MC_part1_factor;      // MC_part1_factor
        delete m_MC_part1_ID;          // MC_part1_ID
        delete m_MC_part2_factor;      // MC_part2_factor
        delete m_MC_part2_ID;          // MC_part2_ID
        delete m_MC_t_pt;              // MC_t_pt
        delete m_MC_t_eta;             // MC_t_eta
        delete m_MC_t_phi;             // MC_t_phi
        delete m_MC_t_e;               // MC_t_E
        delete m_MC_tbar_pt;           // MC_tbar_pt
        delete m_MC_tbar_eta;          // MC_tbar_eta
        delete m_MC_tbar_phi;          // MC_tbar_phi
        delete m_MC_tbar_e;            // MC_tbar_E
        delete m_MC_lep_pt;            // MC_lep_pt
        delete m_MC_lep_eta;           // MC_lep_eta
        delete m_MC_lep_phi;           // MC_lep_phi
        delete m_MC_lep_e;             // MC_lep_E
        delete m_MC_lep_ID;            // MC_lep_ID
        delete m_MC_nu_pt;             // MC_nu_pt
        delete m_MC_nu_eta;            // MC_nu_eta
        delete m_MC_nu_phi;            // MC_nu_phi
        delete m_MC_nu_e;              // MC_nu_E
        delete m_MC_lepb_pt;           // MC_lepb_pt
        delete m_MC_lepb_eta;          // MC_lepb_eta
        delete m_MC_lepb_phi;          // MC_lepb_phi
        delete m_MC_lepb_e;            // MC_lepb_E
        delete m_MC_hadW_pt;           // MC_hadW_pt
        delete m_MC_hadW_eta;          // MC_hadW_eta
        delete m_MC_hadW_phi;          // MC_hadW_phi
        delete m_MC_hadW_e;            // MC_hadW_E
        delete m_MC_hadb_pt;           // MC_hadb_pt
        delete m_MC_hadb_eta;          // MC_hadb_eta
        delete m_MC_hadb_phi;          // MC_hadb_phi
        delete m_MC_hadb_e;            // MC_hadb_E
        delete m_MC_cstar;             // MC_cstar
        delete m_MC_x_F;               // MC_x_F
        delete m_MC_Mtt;               // MC_Mtt
      } // end useTruth
    } // end isMC

    return;
}

// THE END
