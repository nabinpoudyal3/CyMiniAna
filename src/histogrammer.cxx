/*
Created:        --
Last Updated:   29 August    2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Make histograms for systematic uncertainties (& nominal) 
to go into plots || TRexFitter

*/
#include "Analysis/CyMiniAna/interface/histogrammer.h"


histogrammer::histogrammer( configuration& cmaConfig, std::string name ) :
  m_config(&cmaConfig),
  m_name(name),
  m_putOverflowInLastBin(true),
  m_putUnderflowInFirstBin(true){
    m_map_histograms1D.clear();
    m_map_histograms2D.clear();
    m_map_histograms3D.clear();

    m_isMC  = m_config->isMC();

    m_useJets      = m_config->useJets();
    m_useLjets     = m_config->useLargeRJets();
    m_useLeptons   = m_config->useLeptons();
    m_useNeutrinos = m_config->useNeutrinos();

    if (m_name.length()>0  && m_name.substr(m_name.length()-1,1).compare("_")!=0)
        m_name = m_name+"_"; // add '_' to end of string, if needed

    m_mapContainmentRev = m_config->mapOfPartonContainmentRev();
    m_containments.clear();
    for (const auto& x : m_config->mapOfPartonContainment())
        m_containments.push_back(x.first);
  }

histogrammer::~histogrammer() {}


/**** INITIALIZE HISTOGRAMS ****/

// -- 1D Histograms
void histogrammer::init_hist( const std::string& name, const unsigned int nBins, const double x_min, const double x_max ){
    /* Initialize histogram -- equal bins */
    m_map_histograms1D["h_"+name] = new TH1D(("h_"+name).c_str(), ("h_"+name).c_str(),nBins,x_min,x_max);
    m_map_histograms1D["h_"+name]->Sumw2();

    return;
}
void histogrammer::init_hist( const std::string& name, const unsigned int nBins, const double *xbins ){
    /* Initialize histogram -- variable bins */
    m_map_histograms1D["h_"+name] = new TH1D(("h_"+name).c_str(), ("h_"+name).c_str(),nBins,xbins);
    m_map_histograms1D["h_"+name]->Sumw2();

    return;
}
// -- 2D Histograms
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double x_min, const double x_max,
                              const unsigned int nBinsY, const double y_min, const double y_max ){
    /* Initialize histogram -- equal bins */
    m_map_histograms2D["h_"+name] = new TH2D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                            nBinsX,x_min,x_max,nBinsY,y_min,y_max);
    m_map_histograms2D["h_"+name]->Sumw2();

    return;
}
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double *xbins,
                              const unsigned int nBinsY, const double *ybins ){
    /* Initialize histogram -- variable bins */
    m_map_histograms2D["h_"+name] = new TH2D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                           nBinsX,xbins,nBinsY,ybins);
    m_map_histograms2D["h_"+name]->Sumw2();

    return;
}
// -- 3D Histograms
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double x_min, const double x_max,
                              const unsigned int nBinsY, const double y_min, const double y_max,
                              const unsigned int nBinsZ, const double z_min, const double z_max ){
    /* Initialize histogram -- equal bins */
    m_map_histograms3D["h_"+name] = new TH3D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                            nBinsX,x_min,x_max,nBinsY,y_min,y_max,nBinsZ,z_min,z_max);
    m_map_histograms3D["h_"+name]->Sumw2();

    return;
}
void histogrammer::init_hist( const std::string& name, const unsigned int nBinsX, const double *xbins,
                              const unsigned int nBinsY, const double *ybins,
                              const unsigned int nBinsZ, const double *zbins ){
    /* Initialize histogram -- variable bins */
    m_map_histograms3D["h_"+name] = new TH3D(("h_"+name).c_str(), ("h_"+name).c_str(),
                                           nBinsX,xbins,nBinsY,ybins,nBinsZ,zbins);
    m_map_histograms3D["h_"+name]->Sumw2();

    return;
}


void histogrammer::initialize( TFile& outputFile, bool doSystWeights ){
    /* Setup some values and book histograms */
    m_doSystWeights = doSystWeights;
    outputFile.cd();


    // loop over selections (typically only one treename)
    for (const auto& sel : m_config->selections() ){
        bookHists( m_name+sel );

        // weight systematics
        if (m_isMC && m_doSystWeights){
            for (const auto& syst : m_config->listOfWeightSystematics()){
                bookHists( m_name+sel+"_"+syst );
            } // end weight systematics

            // vector weight systematics
            for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
                for (unsigned int el=0;el<syst.second;++el){
                    std::string weightIndex = std::to_string(el);
                    bookHists( m_name+sel+"_"+weightIndex+"_"+syst.first );
                } // end components of vector
            } // end vector weight systematics
        } // end if MC and save weight systematics
    } // end loop over selections

    return;
}

void histogrammer::bookHists( std::string name ){
    /* 
      Book histograms -- modify/inherit this function for analysis-specific hists 

      @param name   This is the string used to identify histograms for different systematics/event weights
    */
    m_names.resize(0); // append names to this to keep track of later

    if (m_useJets){
        init_hist("n_jets_"+name,   31, -0.5,  30.5);
        init_hist("n_btags_"+name,  11, -0.5,  10.5);

        init_hist("jet_pt_"+name,     2000,  0.0, 2000.0);
        init_hist("jet_eta_"+name,      50, -2.5,    2.5);
        init_hist("jet_phi_"+name,      64, -3.2,    3.2);
        init_hist("jet_bdisc_"+name,   100,  0.0,    1.0);
    }

    if (m_useLjets){
        init_hist("n_ljets_"+name,       31, -0.5,   30.5);
        init_hist("ljet_pt_"+name,     2000,  0.0, 2000.0);
        init_hist("ljet_eta_"+name,      50, -2.5,    2.5);
        init_hist("ljet_phi_"+name,      64, -3.2,    3.2);
        init_hist("ljet_SDmass_"+name,  500,  0.0,  500.0);
        init_hist("ljet_charge_"+name, 1000, -5.0,    5.0);
        init_hist("ljet_BEST_t_"+name,  100,  0.0,    1.0);
        init_hist("ljet_BEST_w_"+name,  100,  0.0,    1.0);
        init_hist("ljet_BEST_z_"+name,  100,  0.0,    1.0);
        init_hist("ljet_BEST_h_"+name,  100,  0.0,    1.0);
        init_hist("ljet_BEST_j_"+name,  100,  0.0,    1.0);
        init_hist("ljet_BEST_t_j_"+name,100,  0.0,    1.0);
        init_hist("ljet_tau1_"+name,    200,  0.0,    2.0);
        init_hist("ljet_tau2_"+name,    200,  0.0,    2.0);
        init_hist("ljet_tau3_"+name,    200,  0.0,    2.0);
        init_hist("ljet_tau21_"+name,   100,  0.0,    1.0);
        init_hist("ljet_tau32_"+name,   100,  0.0,    1.0);
        init_hist("ljet_subjet0_bdisc_"+name, 100, 0.0, 1.0);
        init_hist("ljet_subjet1_bdisc_"+name, 100, 0.0, 1.0);
        init_hist("ljet_subjet0_charge_"+name,1000, -5.0, 5.0);
        init_hist("ljet_subjet1_charge_"+name,1000, -5.0, 5.0);
        init_hist("ljet_subjet0_tau21_"+name,100, 0.0, 1.0);
        init_hist("ljet_subjet1_tau21_"+name,100, 0.0, 1.0);
        init_hist("ljet_subjet0_tau32_"+name,100, 0.0, 1.0);
        init_hist("ljet_subjet1_tau32_"+name,100, 0.0, 1.0);

        init_hist("ljet_subjet-b_subjet-w_tau21_"+name,100, 0.0,1.0, 100,0.0,1.0);  // subjet b (=higher bdisc) tau21 vs subjet w (=lower bdisc) tau21
        init_hist("ljet_subjet0_subjet1_tau21_"+name,100, 0.0,1.0, 100,0.0,1.0);  // subjet0 tau21 vs subjet1 tau21
        init_hist("ljet_subjet0_subjet1_tau32_"+name,100, 0.0,1.0, 100,0.0,1.0);  // subjet0 tau32 vs subjet1 tau32
        init_hist("ljet_subjet0_charge_bdisc_"+name,1000,-5.0,5.0, 100,0.0,1.0);  // charge vs bdisc (charge=x-axis)
        init_hist("ljet_subjet1_charge_bdisc_"+name,1000,-5.0,5.0, 100,0.0,1.0);  // charge vs bdisc (charge=x-axis)

        init_hist("ljet_pt_eta_"+name,    200,  0.0, 2000.0,  50, -2.5, 2.5);  // pt vs eta (pt=x-axis)
        init_hist("ljet_pt_SDmass_"+name, 200,  0.0, 2000.0,  50,  0, 500);    // pt vs SDmass (pt=x-axis)

        init_hist("ljet_SDmass_tau32_"+name, 500,  0.0,  500.0, 100,  0,   1);    // SDmass  vs tau32  (SDmass=x-axis)
        init_hist("ljet_BEST_t_SDmass_"+name,100,  0.0,    1.0, 500,  0, 500);    // BEST(t) vs SDmass (BEST=x-axis)
        init_hist("ljet_BEST_t_tau32_"+name, 100,  0.0,    1.0, 100,  0,   1);    // BEST(t) vs tau32  (BEST=x-axis)

        // plots of subjet charge/bdisc for different lepton charges
        if (m_useLeptons){
            init_hist("ljet_charge_Qpos_"+name,1000, -5.0, 5.0);
            init_hist("ljet_charge_Qneg_"+name,1000, -5.0, 5.0);
            init_hist("ljet_subjet_0_charge_Qpos_"+name,1000, -5.0, 5.0);
            init_hist("ljet_subjet_0_bdisc_Qpos_"+name,  100,  0.0, 1.0);
            init_hist("ljet_subjet_1_charge_Qpos_"+name,1000, -5.0, 5.0);
            init_hist("ljet_subjet_1_bdisc_Qpos_"+name,  100,  0.0, 1.0);
            init_hist("ljet_subjet_0_charge_Qneg_"+name,1000, -5.0, 5.0);
            init_hist("ljet_subjet_0_bdisc_Qneg_"+name,  100,  0.0, 1.0);
            init_hist("ljet_subjet_1_charge_Qneg_"+name,1000, -5.0, 5.0);
            init_hist("ljet_subjet_1_bdisc_Qneg_"+name,  100,  0.0, 1.0);
        }

        if (m_config->isTtbar()){
            for (const auto& c : m_containments){
                std::string cname = c+"_"+name;
                init_hist("ljet_pt_"+cname,     2000,  0.0, 2000.0);
                init_hist("ljet_SDmass_"+cname,  500,  0.0,  500.0);
                init_hist("ljet_BEST_t_"+cname,  100,  0.0,    1.0);
                init_hist("ljet_BEST_w_"+cname,  100,  0.0,    1.0);
                init_hist("ljet_BEST_z_"+cname,  100,  0.0,    1.0);
                init_hist("ljet_BEST_h_"+cname,  100,  0.0,    1.0);
                init_hist("ljet_BEST_j_"+cname,  100,  0.0,    1.0);
                init_hist("ljet_BEST_t_j_"+cname,100,  0.0,    1.0);
                init_hist("ljet_tau21_"+cname,   100,  0.0,    1.0);
                init_hist("ljet_tau32_"+cname,   100,  0.0,    1.0);

                init_hist("ljet_pt_SDmass_"+cname,    200,  0.0, 2000.0,  50,  0, 500);    // pt vs SDmass (pt=x-axis)
                init_hist("ljet_SDmass_tau32_"+cname, 500,  0.0,  500.0, 100,  0,   1);    // SDmass  vs tau32  (SDmass=x-axis)
                init_hist("ljet_BEST_t_SDmass_"+cname,100,  0.0,    1.0, 500,  0, 500);    // BEST(t) vs SDmass (BEST=x-axis)
                init_hist("ljet_BEST_t_tau32_"+cname, 100,  0.0,    1.0, 100,  0,   1);    // BEST(t) vs tau32  (BEST=x-axis)
            } // end loop over containments
        } // end if isttbar
    }

    if (m_useLeptons){
        init_hist("el_pt_"+name,    500, 0.0,2000);
        init_hist("el_eta_"+name,    50,-2.5, 2.5);
        init_hist("el_phi_"+name,    64,-3.2, 3.2);
        init_hist("el_charge_"+name,240,-1.2, 1.2);
        init_hist("el_ptrel_"+name, 500, 0.0, 500);
        init_hist("el_drmin_"+name,  50, 0.0,   5);
        init_hist("n_el_"+name,      10,   0,  10);

        init_hist("mu_pt_"+name,    500, 0.0,2000);
        init_hist("mu_eta_"+name,    50,-2.5, 2.5);
        init_hist("mu_phi_"+name,    64,-3.2, 3.2);
        init_hist("mu_charge_"+name,240,-1.2, 1.2);
        init_hist("mu_ptrel_"+name, 500, 0.0, 500);
        init_hist("mu_drmin_"+name,  50, 0.0,   5);
        init_hist("n_mu_"+name,      10,   0,  10);
    }

    if (m_useNeutrinos){
        init_hist("nu_pt_"+name,  500, 0.0, 2000);
        init_hist("nu_eta_"+name,  50, -2.5, 2.5);
        init_hist("nu_phi_"+name,  64, -3.2, 3.2);
    }

    // kinematics
    init_hist("met_met_"+name, 500,  0.0,  500);
    init_hist("met_phi_"+name,  64, -3.2,  3.2);
    init_hist("ht_"+name,     5000,  0.0, 5000);
    init_hist("mtw_"+name,     500,  0.0,  500);

    // ttbar reconstruction + asymmetry values
    if (m_config->isOneLeptonAnalysis()){
        init_hist("deltaR_lep_ak4_"+name, 50, 0.0, 5.0);
        init_hist("pTrel_lep_ak4_"+name, 100, 0.0, 500);
        init_hist("deltaR_lep_ak8_"+name, 50, 0.0, 5.0);
        init_hist("deltaR_ak4_ak8_"+name, 50, 0.0, 5.0);
        init_hist("deltaR_pTrel_lep_ak4_"+name, 50,0.0,5.0, 100,0.0,500);

        init_hist("hadtop_pt_"+name,     2000,  0.0, 2000.0);
        init_hist("hadtop_eta_"+name,      50, -2.5,    2.5);
        init_hist("hadtop_SDmass_"+name,  500,  0.0,  500.0);
        init_hist("hadtop_tau21_"+name,   100,  0.0,    1.0);
        init_hist("hadtop_tau32_"+name,   100,  0.0,    1.0);

        init_hist("deltay_"+name,  1000,-5.0,  5.0);
        init_hist("mttbar_"+name,  5000, 0.0, 5000);
        init_hist("pTttbar_"+name,  300, 0.0,  600);
        init_hist("yttbar_"+name,   100,  0.,   10);
        init_hist("betattbar_"+name,100,  0.,    1);

        init_hist("mttbar_deltay_"+name,  5000, 0.0, 5000, 1000,-5.0,  5.0);  // m_ttbar = x-axis; deltay = y-axis
        init_hist("pTttbar_deltay_"+name,  300, 0.0,  600, 1000,-5.0,  5.0);
        init_hist("yttbar_deltay_"+name,   100,  0.,   10, 1000,-5.0,  5.0);
        init_hist("betattbar_deltay_"+name, 100, 0.0, 1.0, 1000,-5.0,  5.0);

        if (m_config->isTtbar()){
            init_hist("resmat_"+name,          100,-5.0,  5.0, 100,-5.0,  5.0);  // reco=x-axis; truth=y-axis
            init_hist("deltay_dyres_"+name,   1000,-5.0,  5.0, 200,-10,10);
            init_hist("mttbar_dyres_"+name,   5000, 0.0, 5000, 200,-10,10);
            init_hist("pTttbar_dyres_"+name,   300, 0.0,  600, 200,-10,10);
            init_hist("betattbar_dyres_"+name, 100,  0.,    1, 200,-10,10);
            init_hist("yttbar_dyres_"+name,    100,  0.,   10, 200,-10,10);
        }
    }

    return;
}




/**** FILL HISTOGRAMS ****/

void histogrammer::fill( const std::string& name, const double& value, const double& weight ){
    /* TH1D */
    if (!m_map_histograms1D.count("h_"+name)){
        cma::ERROR("HISTOGRAMMER : Filling 1D histogram with key '"+name+"': KEY DOES NOT EXIST");
        cma::ERROR("HISTOGRAMMER : Please check 'init_hist' and 'fill' functions");
    }

    TH1D* this_hist = m_map_histograms1D.at("h_"+name);
    this_hist->Fill(value,weight);

    return;
}
void histogrammer::fill( const std::string& name, 
                         const double& xvalue, const double& yvalue, const double& weight ){
    /* TH2D */
    if (!m_map_histograms2D.count("h_"+name)){
        cma::ERROR("HISTOGRAMMER : Filling 2D histogram with key '"+name+"': KEY DOES NOT EXIST");
        cma::ERROR("HISTOGRAMMER : Please check 'init_hist' and 'fill' functions");
    }

    TH2D* this_hist = m_map_histograms2D.at("h_"+name);
    this_hist->Fill(xvalue,yvalue,weight);

    return;
}
void histogrammer::fill( const std::string& name, 
                         const double& xvalue, const double& yvalue, const double& zvalue, const double& weight ){
    /* TH3D */
    if (!m_map_histograms3D.count("h_"+name)){
        cma::ERROR("HISTOGRAMMER : Filling 3D histogram with key '"+name+"': KEY DOES NOT EXIST");
        cma::ERROR("HISTOGRAMMER : Please check 'init_hist' and 'fill' functions");
    }

    TH3D* this_hist = m_map_histograms3D.at("h_"+name);
    this_hist->Fill(xvalue,yvalue,zvalue,weight);

    return;
}


void histogrammer::fill( Event& event, const std::vector<unsigned int>& evtsel_decisions ){
    /* Fill histograms -- fill histograms based on selection, tree, or systematic weights ("nominal" but different weight)
       This is the function to modify / inherit for analysis-specific purposes
    */
    double event_weight = event.nominal_weight();

    std::vector<std::string> selections = m_config->selections();
    for (unsigned int ss=0, size=selections.size(); ss<size; ss++){
        std::string sel( selections.at(ss) );
        if (!evtsel_decisions.at(ss)) continue;
        fill( m_name+sel, event, event_weight );

        // if there are systematics stored as weights (e.g., b-tagging, pileup, etc.)
        // the following calls the fill() function with different event weights
        // to make histograms
        bool isNominal = m_config->isNominalTree( event.treeName() );
        if (m_isMC && isNominal && m_doSystWeights){
            // weight systematics
            event_weight = 1.0;
            for (const auto& syst : m_config->listOfWeightSystematics()){
                event_weight = event.getSystEventWeight( syst );
                fill( m_name+sel+"_"+syst, event, event_weight );
            } // end weight systematics

            // vector weight systematics
            event_weight = 1.0;
            for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
                for (unsigned int el=0;el<syst.second;++el){
                    event_weight = event.getSystEventWeight( syst.first, el );
                    std::string weightIndex = std::to_string(el);

                    fill( m_name+sel+"_"+weightIndex+"_"+syst.first, event, event_weight );
                } // end components of vector
            } // end vector weight systematics
        } // end if nominal and doSystWeights
    } // end loop over selections

    return;
}


void histogrammer::fill( const std::string& name, Event& event, double event_weight){
    /* Fill histograms -- just use information from the event and fill histogram
       This is the function to modify / inherit for analysis-specific purposes

       //std::vector<Muon> muons = event.muons();  -- merged into "Lepton"
       //std::vector<Electron> electrons = event.electrons();
    */
    cma::DEBUG("HISTOGRAMMER : Fill histograms.");
    cma::DEBUG("HISTOGRAMMER : event weight = "+std::to_string(event_weight) );

    // physics information
    std::vector<Jet> jets = event.jets();
    std::vector<Ljet> ljets = event.ljets();
    std::vector<Lepton> leptons = event.leptons();
    std::vector<Neutrino> neutrinos = event.neutrinos();
    MET met = event.met();

    // fill histograms!
    // -- only filling the "good" objects

    if (m_useJets){
        cma::DEBUG("HISTOGRAMMER : Fill small-R jets");
        fill("n_btags_"+name, event.btag_jets().size(), event_weight );

        fill("n_jets_"+name, jets.size(), event_weight );

        for (const auto& jet : jets){
            if (!jet.isGood) continue;
            fill("jet_pt_"+name,  jet.p4.Pt(),   event_weight);
            fill("jet_eta_"+name, jet.p4.Eta(),  event_weight);
            fill("jet_phi_"+name, jet.p4.Phi(),  event_weight);
            fill("jet_bdisc_"+name, jet.bdisc,  event_weight);
        }
    }


    if (m_useLjets){
        cma::DEBUG("HISTOGRAMMER : Fill large-R jets");
        fill("n_ljets_"+name, ljets.size(), event_weight );

        for (const auto& ljet : ljets){
            if (!ljet.isGood) continue;
            fill("ljet_pt_"+name,    ljet.p4.Pt(),  event_weight);
            fill("ljet_eta_"+name,   ljet.p4.Eta(), event_weight);
            fill("ljet_phi_"+name,   ljet.p4.Phi(), event_weight);
            fill("ljet_SDmass_"+name,ljet.softDropMass, event_weight);
            fill("ljet_charge_"+name,ljet.charge,event_weight);

            fill("ljet_BEST_t_"+name,  ljet.BEST_t,  event_weight);
            fill("ljet_BEST_w_"+name,  ljet.BEST_w,  event_weight);
            fill("ljet_BEST_z_"+name,  ljet.BEST_z,  event_weight);
            fill("ljet_BEST_h_"+name,  ljet.BEST_h,  event_weight);
            fill("ljet_BEST_j_"+name,  ljet.BEST_j,  event_weight);
            fill("ljet_BEST_t_j_"+name,  ljet.BEST_t / (ljet.BEST_t+ljet.BEST_j),  event_weight);

            fill("ljet_tau1_"+name,  ljet.tau1,  event_weight);
            fill("ljet_tau2_"+name,  ljet.tau2,  event_weight);
            fill("ljet_tau3_"+name,  ljet.tau3,  event_weight);
            fill("ljet_tau21_"+name, ljet.tau21, event_weight);
            fill("ljet_tau32_"+name, ljet.tau32, event_weight);
            fill("ljet_subjet0_bdisc_"+name, ljet.subjet0_bdisc, event_weight);
            fill("ljet_subjet1_bdisc_"+name, ljet.subjet1_bdisc, event_weight);
            fill("ljet_subjet0_charge_"+name,ljet.subjet0_charge,event_weight);
            fill("ljet_subjet1_charge_"+name,ljet.subjet1_charge,event_weight);

            float subjet0_tau21 = ljet.subjet0_tau2/ljet.subjet0_tau1;
            float subjet0_tau32 = ljet.subjet0_tau3/ljet.subjet0_tau2;
            float subjet1_tau21 = ljet.subjet1_tau2/ljet.subjet1_tau1;
            float subjet1_tau32 = ljet.subjet1_tau3/ljet.subjet1_tau2;
            fill("ljet_subjet0_tau21_"+name, subjet0_tau21, event_weight);
            fill("ljet_subjet1_tau21_"+name, subjet1_tau21, event_weight);
            fill("ljet_subjet0_tau32_"+name, subjet0_tau32, event_weight);
            fill("ljet_subjet1_tau32_"+name, subjet1_tau32, event_weight);

            // fill based on which subjet has a higher b-disc value
            if (ljet.subjet0_bdisc>ljet.subjet1_bdisc)
                fill("ljet_subjet-b_subjet-w_tau21_"+name, subjet0_tau21, subjet1_tau21, event_weight);
            else
                fill("ljet_subjet-b_subjet-w_tau21_"+name, subjet1_tau21, subjet0_tau21, event_weight);

            fill("ljet_subjet0_subjet1_tau21_"+name,subjet0_tau21,subjet1_tau21,event_weight);  // subjet0 tau21 vs subjet1 tau21
            fill("ljet_subjet0_subjet1_tau32_"+name,subjet0_tau32,subjet1_tau32,event_weight);  // subjet0 tau32 vs subjet1 tau32

            fill("ljet_subjet0_charge_bdisc_"+name, ljet.subjet0_charge, ljet.subjet0_bdisc, event_weight);
            fill("ljet_subjet1_charge_bdisc_"+name, ljet.subjet1_charge, ljet.subjet1_bdisc, event_weight);

            fill("ljet_pt_eta_"+name,   ljet.p4.Pt(), ljet.p4.Eta(), event_weight);
            fill("ljet_pt_SDmass_"+name,ljet.p4.Pt(), ljet.softDropMass, event_weight);

            fill("ljet_SDmass_tau32_"+name, ljet.softDropMass, ljet.tau32,  event_weight);    // SDmass  vs tau32  (SDmass=x-axis)
            fill("ljet_BEST_t_SDmass_"+name,ljet.BEST_t, ljet.softDropMass, event_weight);    // BEST(t) vs SDmass (BEST=x-axis)
            fill("ljet_BEST_t_tau32_"+name, ljet.BEST_t, ljet.tau32,        event_weight);    // BEST(t) vs tau32  (BEST=x-axis)

            int charge(-999);
            if (m_useLeptons && leptons.size()>0) {
                // only interested in these plots for lepton+jets channel (1 lepton reconstructed)
                charge = leptons.at(0).charge;

                if (charge>0) {
                    fill("ljet_charge_Qpos_"+name, ljet.charge, event_weight);
                    fill("ljet_subjet_0_charge_Qpos_"+name, ljet.subjet0_charge,event_weight);
                    fill("ljet_subjet_0_bdisc_Qpos_"+name,  ljet.subjet0_bdisc, event_weight);
                    fill("ljet_subjet_1_charge_Qpos_"+name, ljet.subjet1_charge,event_weight);
                    fill("ljet_subjet_1_bdisc_Qpos_"+name,  ljet.subjet1_bdisc, event_weight);
                }
                else {
                    fill("ljet_charge_Qneg_"+name, ljet.charge, event_weight);
                    fill("ljet_subjet_0_charge_Qneg_"+name, ljet.subjet0_charge,event_weight);
                    fill("ljet_subjet_0_bdisc_Qneg_"+name,  ljet.subjet0_bdisc, event_weight);
                    fill("ljet_subjet_1_charge_Qneg_"+name, ljet.subjet1_charge,event_weight);
                    fill("ljet_subjet_1_bdisc_Qneg_"+name,  ljet.subjet1_bdisc, event_weight);
                }
            } // end if use leptons

            if (m_config->isTtbar()){
                std::string cname = m_mapContainmentRev[ std::abs(ljet.containment) ]+"_"+name;
                fill("ljet_pt_"+cname,     ljet.p4.Pt(),      event_weight);
                fill("ljet_SDmass_"+cname, ljet.softDropMass, event_weight);
                fill("ljet_BEST_t_"+cname,  ljet.BEST_t,       event_weight);
                fill("ljet_BEST_w_"+cname,  ljet.BEST_w,       event_weight);
                fill("ljet_BEST_z_"+cname,  ljet.BEST_z,       event_weight);
                fill("ljet_BEST_h_"+cname,  ljet.BEST_h,       event_weight);
                fill("ljet_BEST_j_"+cname,  ljet.BEST_j,       event_weight);
                fill("ljet_BEST_t_j_"+cname,ljet.BEST_t / (ljet.BEST_t+ljet.BEST_j), event_weight);
                fill("ljet_tau21_"+cname, ljet.tau21, event_weight);
                fill("ljet_tau32_"+cname, ljet.tau32, event_weight);

                fill("ljet_pt_SDmass_"+cname, ljet.p4.Pt(),ljet.softDropMass,event_weight);

                fill("ljet_SDmass_tau32_"+cname, ljet.softDropMass, ljet.tau32,  event_weight);    // SDmass  vs tau32  (SDmass=x-axis)
                fill("ljet_BEST_t_SDmass_"+cname,ljet.BEST_t, ljet.softDropMass, event_weight);    // BEST(t) vs SDmass (BEST=x-axis)
                fill("ljet_BEST_t_tau32_"+cname, ljet.BEST_t, ljet.tau32,        event_weight);    // BEST(t) vs tau32  (BEST=x-axis)
            } // end if isttbar
        } // end loop over ljets
    } // end if use ljets

    if (m_useLeptons){
        cma::DEBUG("HISTOGRAMMER : Fill leptons");
        unsigned int n_electrons(0);
        unsigned int n_muons(0);
        for (const auto& lep : leptons){
            if (!lep.isGood) continue;

            if (lep.isElectron){
                fill("el_pt_"+name,  lep.p4.Pt(),  event_weight);
                fill("el_eta_"+name, lep.p4.Eta(), event_weight);
                fill("el_phi_"+name, lep.p4.Phi(), event_weight);
                fill("el_charge_"+name, lep.charge, event_weight);
                fill("el_ptrel_"+name,  lep.ptrel,  event_weight);
                fill("el_drmin_"+name,  lep.drmin,  event_weight);
                n_electrons++;
            }
            else if (lep.isMuon){
                fill("mu_pt_"+name,  lep.p4.Pt(),  event_weight);
                fill("mu_eta_"+name, lep.p4.Eta(), event_weight);
                fill("mu_phi_"+name, lep.p4.Phi(), event_weight);
                fill("mu_charge_"+name, lep.charge, event_weight);
                fill("mu_ptrel_"+name,  lep.ptrel,  event_weight);
                fill("mu_drmin_"+name,  lep.drmin,  event_weight);
                n_muons++;
            }
        } // end loop over leptons
        fill("n_el_"+name,  n_electrons, event_weight);
        fill("n_mu_"+name,  n_muons,     event_weight);
    } // end if use leptons

    if (m_useNeutrinos){
        cma::DEBUG("HISTOGRAMMER : Fill neutrinos");
        for (const auto& nu : neutrinos){
            fill("nu_pt_"+name,  nu.p4.Pt(),  event_weight);
            fill("nu_eta_"+name, nu.p4.Eta(), event_weight);
            fill("nu_phi_"+name, nu.p4.Phi(), event_weight);
        }
    }

    // kinematics
    cma::DEBUG("HISTOGRAMMER : Fill kinematics");
    fill("met_met_"+name, met.p4.Pt(),  event_weight);
    fill("met_phi_"+name, met.p4.Phi(), event_weight);
    fill("ht_"+name,      event.HT(),   event_weight);
    fill("mtw_"+name,     met.mtw,      event_weight);


    if (m_config->isOneLeptonAnalysis()){
        // have the asymmetry readily available in boosted single lepton
        cma::DEBUG("HISTOGRAMMER : Fill ttbar AC values");
        Ttbar1L ttbarSL = event.ttbar1L();

        Jet tt_jet     = ttbarSL.jet;
        Ljet tt_ljet   = ttbarSL.ljet;
        Lepton tt_lep  = ttbarSL.lepton;
        Neutrino tt_nu = ttbarSL.neutrino;

        TLorentzVector top_lep;
        TLorentzVector top_had;
        TLorentzVector ttbar;
        top_lep = (tt_nu.p4 + tt_lep.p4 + tt_jet.p4);
        top_had = tt_ljet.p4;
        ttbar   = top_had+top_lep;

        float dr_lep_ak4    = tt_jet.p4.DeltaR( tt_lep.p4 );
        float ptrel_lep_ak4 = cma::ptrel( tt_lep.p4, tt_jet.p4);

        fill("deltaR_lep_ak4_"+name, dr_lep_ak4,    event_weight);
        fill("pTrel_lep_ak4_"+name,  ptrel_lep_ak4, event_weight);
        fill("deltaR_lep_ak8_"+name, top_had.DeltaR( tt_lep.p4 ), event_weight);
        fill("deltaR_ak4_ak8_"+name, top_had.DeltaR( tt_jet.p4 ), event_weight);
        fill("deltaR_pTrel_lep_ak4_"+name, dr_lep_ak4, ptrel_lep_ak4, event_weight);

        fill("hadtop_pt_"+name,     top_had.Pt(), event_weight);
        fill("hadtop_eta_"+name,    top_had.Eta(), event_weight);
        fill("hadtop_SDmass_"+name, tt_ljet.softDropMass, event_weight);
        fill("hadtop_tau21_"+name,  tt_ljet.tau21, event_weight);
        fill("hadtop_tau32_"+name,  tt_ljet.tau32, event_weight);

        // asymmetry
        float dy     = tt_lep.charge * ( std::abs(top_lep.Rapidity()) - std::abs(top_had.Rapidity()) );
        float mtt    = ttbar.M();
        float pttt   = ttbar.Pt();
        float ytt    = std::abs(ttbar.Rapidity());
        float betatt = std::abs(top_had.Pz() + top_lep.Pz()) / (top_had.E() + top_lep.E());

        fill("deltay_"+name,  dy,   event_weight);
        fill("mttbar_"+name,  mtt,  event_weight);
        fill("pTttbar_"+name, pttt, event_weight);
        fill("yttbar_"+name,  ytt,  event_weight);
        fill("betattbar_"+name, betatt, event_weight);

        fill("mttbar_deltay_"+name,  mtt,  dy, event_weight);
        fill("pTttbar_deltay_"+name, pttt, dy, event_weight);
        fill("yttbar_deltay_"+name,  ytt,  dy, event_weight);
        fill("betattbar_deltay_"+name, betatt, dy, event_weight);

        if (m_config->isTtbar()){
            float true_dy(0);
            std::vector<TruthTop> truth_tops = event.truth();
            std::vector<Parton> truth_partons = event.truth_partons();
            if (truth_tops.size()==2){
                TruthTop top0 = truth_tops.at(0);
                TruthTop top1 = truth_tops.at(1);
                Parton ptop0  = truth_partons.at( top0.Top );
                Parton ptop1  = truth_partons.at( top1.Top );

                true_dy = (top0.isTop && top1.isAntiTop) ? std::abs(ptop0.p4.Rapidity()) - std::abs(ptop1.p4.Rapidity()) : 
                                                           std::abs(ptop1.p4.Rapidity()) - std::abs(ptop0.p4.Rapidity());
                fill("resmat_"+name, dy, true_dy, event_weight);  // x=reco, y=truth

                float dyres = true_dy - dy;
                fill("deltay_dyres_"+name,   dy,    dyres, event_weight);
                fill("mttbar_dyres_"+name,   mtt,   dyres, event_weight);
                fill("pTttbar_dyres_"+name,  pttt,  dyres, event_weight);
                fill("betattbar_dyres_"+name,betatt,dyres, event_weight);
                fill("yttbar_dyres_"+name,   ytt,   dyres, event_weight);
            }
        } // end response matrix creation
    }

    cma::DEBUG("HISTOGRAMMER : End histograms");

    return;
}





/**** OVER/UNDERFLOW ****/

void histogrammer::overUnderFlow(){
    /* Call overflow and underflow functions at once */
    overFlow();
    underFlow();
    return;
}


void histogrammer::overFlow() {
    /* Add overflow to last bin */
    if (!m_putOverflowInLastBin){
        cma::INFO("HISTOGRAMMER : Not putting overflow in last bin(s)");
        return;
    }
    else{
        // loop over 1D histograms
        for (const auto& hist : m_map_histograms1D){
            unsigned int numBins = hist.second->GetNbinsX();
            double overflowContent = hist.second->GetBinContent(numBins + 1);

            hist.second->AddBinContent(numBins,overflowContent); // add overflow to last bin
            hist.second->SetBinContent(numBins+1, 0);            // set overflow to 0
        }
        // loop over 2D histograms
        for (const auto& hist : m_map_histograms2D){
            unsigned int numBinsX = hist.second->GetXaxis()->GetNbins();
            unsigned int numBinsY = hist.second->GetYaxis()->GetNbins();

            // x-axis :: overflow in y
            for (unsigned int xx=1;xx<numBinsX+1;++xx){
                double overflowContent = hist.second->GetBinContent(xx,numBinsY+1);

                int lastBin     = hist.second->GetBin(xx,numBinsY);
                int overflowBin = hist.second->GetBin(xx,numBinsY+1);
                hist.second->AddBinContent(lastBin,overflowContent); // add overflow to last bin
                hist.second->SetBinContent(overflowBin,0);           // set overflow to 0
            }
            // y-axis :: overflow in x
            for (unsigned int yy=1;yy<numBinsY;++yy){
                double overflowContent = hist.second->GetBinContent(numBinsX+1,yy);

                int lastBin     = hist.second->GetBin(numBinsX,yy);
                int overflowBin = hist.second->GetBin(numBinsX+1,yy);
                hist.second->AddBinContent(lastBin,overflowContent); // add overflow to last bin
                hist.second->SetBinContent(overflowBin,0);           // set overflow to 0
            }
        } // end 2D histogram overflow
    } // end else put overflow in first bin

    return;
}

void histogrammer::underFlow() {
    /* Add underflow to first bin */
    if (!m_putUnderflowInFirstBin){
        cma::INFO("HISTOGRAMMER : Not putting underflow in first bin(s)");
        return;
    }
    else{
        // loop over 1D histograms
        for (const auto& hist : m_map_histograms1D){
            double underflowContent = hist.second->GetBinContent(0);

            hist.second->AddBinContent(1, underflowContent);  // add underflow to first bin
            hist.second->SetBinContent(0, 0);                 // set underflow to 0
        }
        // loop over 2D histograms
        for (const auto& hist : m_map_histograms2D){
            unsigned int numBinsX = hist.second->GetXaxis()->GetNbins();
            unsigned int numBinsY = hist.second->GetYaxis()->GetNbins();

            // x-axis :: underflow in y
            for (unsigned int xx=1;xx<numBinsX+1;++xx){
                double underflowContent = hist.second->GetBinContent(xx,numBinsY+1);

                int firstBin     = hist.second->GetBin(xx,1);
                int underflowBin = hist.second->GetBin(xx,0);
                hist.second->AddBinContent(firstBin,underflowContent); // add overflow to last bin
                hist.second->SetBinContent(underflowBin,0);            // set overflow to 0
            }
            // y-axis :: underflow in x
            for (unsigned int yy=1;yy<numBinsY;++yy){
                double underflowContent = hist.second->GetBinContent(0,yy);

                int firstBin     = hist.second->GetBin(1,yy);
                int underflowBin = hist.second->GetBin(0,yy);
                hist.second->AddBinContent(firstBin,underflowContent); // add overflow to last bin
                hist.second->SetBinContent(underflowBin,0);           // set overflow to 0
            }
        } // end 2D histogram underflow
    } // end else put underflow in first bin
}

// THE END
