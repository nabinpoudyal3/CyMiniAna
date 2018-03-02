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
#include "Analysis/CyMiniAna/interface/histogrammer4ML.h"


histogrammer4ML::histogrammer4ML( configuration& cmaConfig, std::string name ) :
  histogrammer::histogrammer(cmaConfig,name),
  m_config(&cmaConfig),
  m_name(name){}

histogrammer4ML::~histogrammer4ML() {}


/**** INITIALIZE HISTOGRAMS ****/


void histogrammer4ML::initialize( TFile& outputFile ){
    /* Setup some values and book histograms */
    outputFile.cd();

    bookHists();

    return;
}

void histogrammer4ML::bookHists(){
    /* 
      Book histograms -- modify/inherit this function for analysis-specific hists 

      @param name   This is the string used to identify histograms for different systematics/event weights

      0 :: NONE = QCD (background)
      1 :: QB-Q = Signal AK8(QB) + AK4(Q)
      2 :: QQ-B = Signal AK8(W)  + AK4(B)
    */
    cma::DEBUG("HISTOGRAMMER : Init. histograms: "+m_name);

    for (const auto& target : m_targets){
        // AK8 mass vs DeltaR(AK8,AK4)
        histogrammer::init_hist("MassvDR_ljet_jet_"+target+"_"+m_name,  5000,0,5000,50,0,5);

        // AK8 N-subjettiness vs DeltaR(AK8,AK4)
//        histogrammer::init_hist("Tau21vDR_ljet_jet_"+target+"_"+m_name,  100,0,1,50,0,5);
//        histogrammer::init_hist("Tau32vDR_ljet_jet_"+target+"_"+m_name,  100,0,1,50,0,5);

        // Features
        for (unsigned int i=0;i<16;i++)
            histogrammer::init_hist("ljet_deepAK8-"+std::to_string(i)+"_"+target+"_"+m_name,  100, 0,1);

        histogrammer::init_hist("jet_bdisc_"+target+"_"+m_name,  200, -1,1);
        histogrammer::init_hist("jet_charge_"+target+"_"+m_name, 500, -5,5);
        histogrammer::init_hist("ljet_jet_deltaR_"+target+"_"+m_name,   50, 0,    5);
        histogrammer::init_hist("ljet_jet_m_"+target+"_"+m_name,      5000, 0, 5000);

/*      ORIGINAL : pre-DEEPAK8
        histogrammer::init_hist("ljet_SDmass-"+target+"_"+m_name,  500,  0.0,  500.0);
        histogrammer::init_hist("ljet_tau1-"+target+"_"+m_name,    200,  0.0,    2.0);
        histogrammer::init_hist("ljet_tau2-"+target+"_"+m_name,    200,  0.0,    2.0);
        histogrammer::init_hist("ljet_tau3-"+target+"_"+m_name,    200,  0.0,    2.0);
        histogrammer::init_hist("ljet_tau21-"+target+"_"+m_name,   100,  0.0,    1.0);
        histogrammer::init_hist("ljet_tau32-"+target+"_"+m_name,   100,  0.0,    1.0);
        histogrammer::init_hist("ljet_subjet0_bdisc-"+target+"_"+m_name, 100, 0.0, 1.0);
        histogrammer::init_hist("ljet_subjet1_bdisc-"+target+"_"+m_name, 100, 0.0, 1.0);
        histogrammer::init_hist("ljet_subjet0_pTrel-"+target+"_"+m_name, 100, 0.0, 1.0);
        histogrammer::init_hist("ljet_subjet1_pTrel-"+target+"_"+m_name, 100, 0.0, 1.0);
*/
    }

    return;
}


/**** FILL HISTOGRAMS ****/
void histogrammer4ML::fill( const std::map<std::string,double> features, double weight ){
    /* Fill histograms -- 
       Fill information from single top object (inputs to deep learning)
    */
    std::string target = std::to_string( int(features.at("target")) );

    cma::DEBUG("HISTOGRAMMER : Fill histograms: "+m_name+"; target = "+target);

    // AK8 N-subjettiness vs DeltaR(AK8,AK4)
    histogrammer::fill("MassvDR_ljet_jet_"+target+"_"+m_name,  features.at("ljet_jet_m"), features.at("ljet_jet_deltaR"), weight );
//    histogrammer::fill("Tau21vDR_ljet_jet_"+target+"_"+m_name, features.at("ljet_tau21"), features.at("ljet_jet_deltaR"), weight );
//    histogrammer::fill("Tau32vDR_ljet_jet_"+target+"_"+m_name, features.at("ljet_tau32"), features.at("ljet_jet_deltaR"), weight );

    // Features
    histogrammer::fill("jet_bdisc_"+target+"_"+m_name,  features.at("jet_bdisc"), weight);
    histogrammer::fill("jet_charge_"+target+"_"+m_name, features.at("jet_charge"), weight);

    cma::DEBUG("HISTOGRAMMER : ljets deepAK8 ");
    for (unsigned int i=0;i<16;i++){
        std::string idx = std::to_string(i);
        histogrammer::fill("ljet_deepAK8-"+idx+"_"+target+"_"+m_name, features.at("ljet_deepAK8_"+idx), weight);
    }

    cma::DEBUG("HISTOGRAMMER : ljets-jets ");
    histogrammer::fill("ljet_jet_deltaR_"+target+"_"+m_name, features.at("ljet_jet_deltaR"), weight);
    histogrammer::fill("ljet_jet_m_"+target+"_"+m_name,   features.at("ljet_jet_m"),   weight);

/*
    histogrammer::fill("ljet_SDmass-"+target+"_"+m_name, features.at("ljet_SDmass"), weight);
    histogrammer::fill("ljet_tau1-"+target+"_"+m_name,   features.at("ljet_tau1"),  weight);
    histogrammer::fill("ljet_tau2-"+target+"_"+m_name,   features.at("ljet_tau2"),  weight);
    histogrammer::fill("ljet_tau3-"+target+"_"+m_name,   features.at("ljet_tau3"),  weight);
    histogrammer::fill("ljet_tau21-"+target+"_"+m_name,  features.at("ljet_tau21"), weight);
    histogrammer::fill("ljet_tau32-"+target+"_"+m_name,  features.at("ljet_tau32"), weight);
    histogrammer::fill("ljet_subjet0_bdisc-"+target+"_"+m_name, features.at("ljet_subjet0_bdisc"), weight);
    histogrammer::fill("ljet_subjet0_pTrel-"+target+"_"+m_name, features.at("ljet_subjet0_pTrel"),  weight);
    histogrammer::fill("ljet_subjet1_bdisc-"+target+"_"+m_name, features.at("ljet_subjet1_bdisc"), weight);
    histogrammer::fill("ljet_subjet1_pTrel-"+target+"_"+m_name, features.at("ljet_subjet1_pTrel"), weight);
*/
    cma::DEBUG("HISTOGRAMMER : End histograms");

    return;
}

// THE END

