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

      0 :: Top      (lepton Q > 0)
      1 :: Anti-top (lepton Q < 0)
    */
    cma::DEBUG("HISTOGRAMMER : Init. histograms: "+m_name);

    for (const auto& target : m_targets){
        // Non-features just to compare consistency between top/anti-top
        histogrammer::init_hist("ljet_BEST_t-"+target+"_"+m_name,  100,  0.0,  100.0);
        histogrammer::init_hist("ljet_BEST_w-"+target+"_"+m_name,  100,  0.0,  100.0);
        histogrammer::init_hist("ljet_BEST_z-"+target+"_"+m_name,  100,  0.0,  100.0);
        histogrammer::init_hist("ljet_BEST_h-"+target+"_"+m_name,  100,  0.0,  100.0);
        histogrammer::init_hist("ljet_BEST_j-"+target+"_"+m_name,  100,  0.0,  100.0);
        histogrammer::init_hist("ljet_SDmass-"+target+"_"+m_name,  500,  0.0,  500.0);
        histogrammer::init_hist("ljet_tau1-"+target+"_"+m_name,    200,  0.0,    2.0);
        histogrammer::init_hist("ljet_tau2-"+target+"_"+m_name,    200,  0.0,    2.0);
        histogrammer::init_hist("ljet_tau3-"+target+"_"+m_name,    200,  0.0,    2.0);
        histogrammer::init_hist("ljet_tau21-"+target+"_"+m_name,   100,  0.0,    1.0);
        histogrammer::init_hist("ljet_tau32-"+target+"_"+m_name,   100,  0.0,    1.0);

        // Features
        histogrammer::init_hist("ljet_charge-"+target+"_"+m_name,  3,-1.5, 1.5);
        histogrammer::init_hist("ljet_subjet0_bdisc-"+target+"_"+m_name, 100, 0.0, 1.0);
        histogrammer::init_hist("ljet_subjet0_pTrel-"+target+"_"+m_name, 100, 0.0, 1.0);
        histogrammer::init_hist("ljet_subjet0_charge-"+target+"_"+m_name,  3,-1.5, 1.5);
        histogrammer::init_hist("ljet_subjet1_bdisc-"+target+"_"+m_name, 100, 0.0, 1.0);
        histogrammer::init_hist("ljet_subjet1_pTrel-"+target+"_"+m_name, 100, 0.0, 1.0);
        histogrammer::init_hist("ljet_subjet1_charge-"+target+"_"+m_name,  3,-1.5, 1.5);
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

    for (const auto& key : features) cma::INFO("HIST4ML : Feature "+key.first);

    histogrammer::fill("ljet_BEST_t-"+target+"_"+m_name, features.at("ljet_BEST_t"), weight);
    histogrammer::fill("ljet_BEST_w-"+target+"_"+m_name, features.at("ljet_BEST_w"), weight);
    histogrammer::fill("ljet_BEST_z-"+target+"_"+m_name, features.at("ljet_BEST_z"), weight);
    histogrammer::fill("ljet_BEST_h-"+target+"_"+m_name, features.at("ljet_BEST_h"), weight);
    histogrammer::fill("ljet_BEST_j-"+target+"_"+m_name, features.at("ljet_BEST_j"), weight);
    histogrammer::fill("ljet_SDmass-"+target+"_"+m_name, features.at("ljet_SDmass"), weight);
    histogrammer::fill("ljet_tau1-"+target+"_"+m_name,   features.at("ljet_tau1"),  weight);
    histogrammer::fill("ljet_tau2-"+target+"_"+m_name,   features.at("ljet_tau2"),  weight);
    histogrammer::fill("ljet_tau3-"+target+"_"+m_name,   features.at("ljet_tau3"),  weight);
    histogrammer::fill("ljet_tau21-"+target+"_"+m_name,  features.at("ljet_tau21"), weight);
    histogrammer::fill("ljet_tau32-"+target+"_"+m_name,  features.at("ljet_tau32"), weight);

    histogrammer::fill("ljet_charge-"+target+"_"+m_name, features.at("ljet_charge"), weight);
    histogrammer::fill("ljet_subjet0_bdisc-"+target+"_"+m_name,  features.at("ljet_subjet0_bdisc"),  weight);
    histogrammer::fill("ljet_subjet0_pTrel-"+target+"_"+m_name,  features.at("ljet_subjet0_pTrel"),  weight);
    histogrammer::fill("ljet_subjet0_charge-"+target+"_"+m_name, features.at("ljet_subjet0_charge"), weight);
    histogrammer::fill("ljet_subjet1_bdisc-"+target+"_"+m_name,  features.at("ljet_subjet1_bdisc"),  weight);
    histogrammer::fill("ljet_subjet1_pTrel-"+target+"_"+m_name,  features.at("ljet_subjet1_pTrel"),  weight);
    histogrammer::fill("ljet_subjet1_charge-"+target+"_"+m_name, features.at("ljet_subjet1_charge"), weight);

    cma::DEBUG("HISTOGRAMMER : End histograms");

    return;
}

// THE END

