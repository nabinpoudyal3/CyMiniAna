#include "diHiggs/CyMiniAna/interface/efficiencyTruth.h"

efficiencyTruth::efficiencyTruth(configuration &cmaConfig) : 
  efficiency::efficiency(cmaConfig){
    m_map_efficiencies.clear();
}


void efficiencyTruth::bookEffs(TFile& outputFile){
    /* Initialize truth efficiencies */
    double bins[] = {-3.,-0.4,0.,0.4,3.};
    unsigned int nbins = sizeof(bins)/sizeof(double)-1;
    init_eff("acceptance_deltay_inclu", nbins, bins);
}


void efficiencyTruth::fill(Event &event, bool decision){
    /* Fill truth efficiencies */
    double truth_deltay(0.0);
    double weight = (event.isValidRecoEntry() ? event.nominal_weight() : 1.0);

    fill("acceptance_deltay_inclu", truth_deltay, decision, weight);
}
