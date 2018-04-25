#ifndef NEUTRINORECO_H
#define NEUTRINORECO_H

#include <string>
#include <vector>
#include <cmath> 

#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/physicsObjects.h"


class neutrinoReco {
  public:
    neutrinoReco( configuration& cmaConfig );

    ~neutrinoReco();

    void setObjects(Lepton& lepton, MET& met);
    void setLepton(Lepton& lepton);
    void setMET(MET& met);
    Neutrino execute(float wmass=80.4);   // build the neutrino assuming W mass [GeV]

    std::vector<float> pzSolutions();

  protected:

    configuration *m_config;

    Neutrino m_nu;
    Lepton m_lepton;
    MET m_met;

    std::vector<float> m_pz_solutions;
};

#endif
