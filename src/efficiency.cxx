/*
Created:         4 September 2016
Last Updated:    4 September 2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109

-----

Make TEfficiencies for measuring efficiencies
https://root.cern.ch/doc/master/classTEfficiency.html

*/
#include "Analysis/CyMiniAna/interface/efficiency.h"


efficiency::efficiency(configuration &cmaConfig) : 
  m_config(&cmaConfig){
   m_map_efficiencies.clear();
  }

efficiency::~efficiency() {}


void efficiency::bookEffs( TFile& outputFile ){
    /* Book efficiencies -- modify/inherit this function for analysis-specific efficiencies */
    m_names.resize(0); // append names to this to keep track of later
    outputFile.cd();

    init_eff("dilution_inclusive",    1, 0, 1);     // delta|y| reconstruction efficiency

    return;
}

// -- Existing efficiencies
void efficiency::init_eff(const TH1 &passed, const TH1 &total){
    /* Initialize efficiency -- existing efficiencies */
    std::string name(total.GetName());
    name+="_clone";
    m_map_efficiencies[name] = new TEfficiency(passed,total);

    return;
} 

// -- 1D efficiencies
void efficiency::init_eff( const std::string &name, const unsigned int nBins, const double x_min, const double x_max ){
    /* Initialize efficiency -- equal bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),nBins,x_min,x_max);
    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}

void efficiency::init_eff( const std::string &name, const unsigned int nBins, const double *xbins ){
    /* Initialize efficiency -- variable bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),nBins,xbins);
    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}

// -- 2D efficiencies
void efficiency::init_eff( const std::string &name, const unsigned int nBinsX, const double x_min, const double x_max,
                              const unsigned int nBinsY, const double y_min, const double y_max ){
    /* Initialize efficiency -- equal bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),
                                               nBinsX,x_min,x_max,nBinsY,y_min,y_max);
    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}

void efficiency::init_eff( const std::string &name, const unsigned int nBinsX, const double *xbins,
                              const unsigned int nBinsY, const double *ybins ){
    /* Initialize efficiency -- variable bins */
    m_map_efficiencies[name] = new TEfficiency((name).c_str(), (name).c_str(),
                                               nBinsX,xbins,nBinsY,ybins);
    m_map_efficiencies[name]->SetUseWeightedEvents();

    return;
}



void efficiency::fill( Event &event ){
    /* Fill efficiencies -- just use information from the event 
       This is the function to modify / inherit for analysis-specific purposes
       Example
       Fill an efficiency for jet trigger vs leading jet pT 
    */
    double reco_deltay  = 0.0;
    double weight       = event.nominal_weight();

    fill("dilution_inclusive", reco_deltay, true, weight);

    return;
}

void efficiency::fill( const std::string &name, const double &value, const bool &decision, const double &weight ){
    /* Fill efficiencies with values! */
    m_map_efficiencies.at(name)->FillWeighted(decision, weight, value);

    return;
}

void efficiency::fill( const std::string &name, 
                         const double &xvalue, const double &yvalue, const bool &decision, const double &weight ){
    /* Fill efficiencies with values! */
    m_map_efficiencies.at(name)->Fill(decision, weight, xvalue, yvalue);

    return;
}

// THE END
