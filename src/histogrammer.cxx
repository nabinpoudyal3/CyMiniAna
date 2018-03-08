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


    // loop over treenames (typically systematic uncertainties)
    for (const auto& treename : m_config->treeNames() )
        bookHists( m_name+treename );

    // weight systematics
    if (m_isMC && m_doSystWeights){
        for (const auto& syst : m_config->listOfWeightSystematics()){
            bookHists( m_name+syst );
        } // end weight systematics

        // vector weight systematics
        for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
            for (unsigned int el=0;el<syst.second;++el){
                std::string weightIndex = std::to_string(el);
                bookHists( m_name+weightIndex+"_"+syst.first );
            } // end components of vector
        } // end vector weight systematics
    } // end if MC and save weight systematics
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
        init_hist("ljet_tau1_"+name,    200,  0.0,    2.0);
        init_hist("ljet_tau2_"+name,    200,  0.0,    2.0);
        init_hist("ljet_tau3_"+name,    200,  0.0,    2.0);
        init_hist("ljet_tau21_"+name,   100,  0.0,    1.0);
        init_hist("ljet_tau32_"+name,   100,  0.0,    1.0);
        init_hist("ljet_subjet0_bdisc_"+name, 100, 0.0, 1.0);
        init_hist("ljet_subjet1_bdisc_"+name, 100, 0.0, 1.0);
        init_hist("ljet_subjet0_charge_"+name,1000, -5.0, 5.0);
        init_hist("ljet_subjet1_charge_"+name,1000, -5.0, 5.0);

        init_hist("ljet_subjet0_charge_bdisc_"+name,1000,-5.0,5.0, 100,0.0,1.0);  // charge vs bdisc (charge=x-axis)
        init_hist("ljet_subjet1_charge_bdisc_"+name,1000,-5.0,5.0, 100,0.0,1.0);  // charge vs bdisc (charge=x-axis)

        init_hist("ljet_pt_eta_"+name,    200,  0.0, 2000.0,  50, -2.5, 2.5);  // pt vs eta (pt=x-axis)
        init_hist("ljet_pt_SDmass_"+name, 200,  0.0, 2000.0,  50,  0, 500);    // pt vs SDmass (pt=x-axis)

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
    }

    if (m_useLeptons){
        init_hist("el_pt_"+name,  500, 0.0, 2000);
        init_hist("el_eta_"+name,  50, -2.5, 2.5);
        init_hist("el_phi_"+name,  64, -3.2, 3.2);
        init_hist("el_charge_"+name, 240, -1.2, 1.2);

        init_hist("mu_pt_"+name,  500, 0.0, 2000);
        init_hist("mu_eta_"+name,  50, -2.5, 2.5);
        init_hist("mu_phi_"+name,  64, -3.2, 3.2);
        init_hist("mu_charge_"+name, 240, -1.2, 1.2);
    }

    if (m_useNeutrinos){
        init_hist("nu_pt_"+name,  500, 0.0, 2000);
        init_hist("nu_eta_"+name,  50, -2.5, 2.5);
        init_hist("nu_phi_"+name,  64, -3.2, 3.2);
    }

    // kinematics
    init_hist("met_met_"+name, 500,  0.0,  500);
    init_hist("met_phi_"+name, 6.4, -3.2,  3.2);
    init_hist("ht_"+name,     5000,  0.0, 5000);

    //init_hist("dnn_"+name,  100, 0.0,   1.);

    // asymmetry variables
    init_hist("mttbar_"+name,  5000, 0.0, 5000);
    init_hist("pTttbar_"+name,  300, 0.0,  300);
    init_hist("yttbar_"+name,   100,-10.,   10);

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


void histogrammer::fill( Event& event ){
    /* Fill histograms -- fill histograms based on treename or systematic weights ("nominal" but different weight)
       This is the function to modify / inherit for analysis-specific purposes
    */
    std::string treeName = event.treeName();
    double event_weight  = event.nominal_weight();

    fill( m_name+treeName, event, event_weight );


    // if there are systematics stored as weights (e.g., b-tagging, pileup, etc.)
    // the following calls the fill() function with different event weights
    // to make histograms
    bool isNominal = m_config->isNominalTree( treeName );
    if (m_isMC && isNominal && m_doSystWeights){
        // weight systematics
        event_weight = 1.0;
        for (const auto& syst : m_config->listOfWeightSystematics()){
            event_weight = event.getSystEventWeight( syst );
            fill( m_name+syst, event, event_weight );
        } // end weight systematics

        // vector weight systematics
        event_weight = 1.0;
        for (const auto& syst : m_config->mapOfWeightVectorSystematics()){
            for (unsigned int el=0;el<syst.second;++el){
                event_weight = event.getSystEventWeight( syst.first, el );
                std::string weightIndex = std::to_string(el);

                fill( m_name+weightIndex+"_"+syst.first, event, event_weight );
            } // end components of vector
        } // end vector weight systematics
    } // end if nominal and doSystWeights

    return;
}


void histogrammer::fill( const std::string& name, Event& event, double event_weight){
    /* Fill histograms -- just use information from the event and fill histogram
       This is the function to modify / inherit for analysis-specific purposes
    */
    cma::DEBUG("HISTOGRAMMER : Fill histograms.");
    cma::DEBUG("HISTOGRAMMER : event weight = "+std::to_string(event_weight) );

    // physics information
    std::vector<Jet> jets = event.jets();
    std::vector<Ljet> ljets = event.ljets();
    std::vector<Muon> muons = event.muons();
    std::vector<Electron> electrons = event.electrons();
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

            fill("ljet_tau2_"+name,  ljet.tau2,  event_weight);
            fill("ljet_tau3_"+name,  ljet.tau3,  event_weight);
            fill("ljet_tau21_"+name, ljet.tau21, event_weight);
            fill("ljet_tau32_"+name, ljet.tau32, event_weight);
            fill("ljet_subjet0_bdisc_"+name, ljet.subjet0_bdisc, event_weight);
            fill("ljet_subjet1_bdisc_"+name, ljet.subjet1_bdisc, event_weight);
            fill("ljet_subjet0_charge_"+name,ljet.subjet0_charge,event_weight);
            fill("ljet_subjet1_charge_"+name,ljet.subjet1_charge,event_weight);

            fill("ljet_subjet0_charge_bdisc_"+name, ljet.subjet0_charge, ljet.subjet0_bdisc, event_weight);
            fill("ljet_subjet1_charge_bdisc_"+name, ljet.subjet1_charge, ljet.subjet1_bdisc, event_weight);

            fill("ljet_pt_eta_"+name,   ljet.p4.Pt(), ljet.p4.Eta(), event_weight);
            fill("ljet_pt_SDmass_"+name,ljet.p4.Pt(), ljet.softDropMass, event_weight);

            int charge(-999);
            if (m_useLeptons) {
                // only interested in these plots for lepton+jets channel (1 lepton reconstructed)
                if (electrons.size()>0)  charge = electrons.at(0).charge;
                else if (muons.size()>0) charge = muons.at(0).charge;

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
        } // end loop over ljets
    } // end if use ljets

    if (m_useLeptons){
        cma::DEBUG("HISTOGRAMMER : Fill leptons");
        for (const auto& el : electrons){
            if (!el.isGood) continue;
            fill("el_pt_"+name,  el.p4.Pt(),  event_weight);
            fill("el_eta_"+name, el.p4.Eta(), event_weight);
            fill("el_phi_"+name, el.p4.Phi(), event_weight);
            fill("el_charge_"+name, el.charge, event_weight);
        }

        for (const auto& mu : muons){
            if (!mu.isGood) continue;
            fill("mu_pt_"+name,  mu.p4.Pt(),  event_weight);
            fill("mu_eta_"+name, mu.p4.Eta(), event_weight);
            fill("mu_phi_"+name, mu.p4.Phi(), event_weight);
            fill("mu_charge_"+name, mu.charge, event_weight);
        }
    }

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

/*
    // DNN
    cma::DEBUG("HISTOGRAMMER : Fill DNN");
    fill("dnn_"+name, event.dnn(), event_weight); // N/A

    // Asymmetry
    cma::DEBUG("HISTOGRAMMER : Fill ttbar AC values");
    fill("mttbar_"+name,  (top.p4+antitop.p4).M(),  event_weight);
    fill("pTttbar_"+name, (top.p4+antitop.p4).Pt(), event_weight);
    fill("yttbar_"+name,  (top.p4+antitop.p4).Rapidity(), event_weight);
*/
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
