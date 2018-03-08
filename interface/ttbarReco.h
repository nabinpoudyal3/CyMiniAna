#ifndef TTBARRECO_H
#define TTBARRECO_H

#include <string>
#include <map>
#include <vector>

#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/physicsObjects.h"


class ttbarReco {
  public:
    ttbarReco( configuration& cmaConfig );

    ~ttbarReco();

    std::vector<Top> tops();
    void execute(const std::vector<Jet>& jets, const std::vector<Ljet>& ljets);
    bool isTopTagged(const Ljet& ljet);

  protected:

    configuration *m_config;

    std::vector<Top> m_ttbar;
    std::map<std::string,int> m_mapContainment;
    std::map<std::string,int> m_targetMap;

    std::vector<Jet> m_jets;
    std::vector<Ljet> m_ljets;
};

#endif
