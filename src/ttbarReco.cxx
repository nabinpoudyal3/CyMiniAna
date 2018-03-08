/*
Created:        --
Last Updated:   19 February 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Tool for building ttbar system
 - all-had:  Two top-tagged AK8
 - l+jets:   1 top-tagged AK8 + 1 leptonic top (AK4+lep+nu)
 - dilepton: 2 leptonic tops (AK4+lep+nu)
*/
#include "Analysis/CyMiniAna/interface/ttbarReco.h"


ttbarReco::ttbarReco( configuration& cmaConfig ) :
  m_config(&cmaConfig){
    m_mapContainment = m_config->mapOfPartonContainment();  // containment map (ints and strings)
    m_targetMap = m_config->mapOfTargetValues();
  }

ttbarReco::~ttbarReco() {}


std::vector<Top> ttbarReco::tops(){
    /* Return the ttbar system */
    return m_ttbar;
}


void ttbarReco::execute(const std::vector<Jet>& jets, const std::vector<Ljet>& ljets){
    /* Build top quarks system */
    m_ttbar.clear();

    m_ljets = ljets;

    bool isQCD(m_config->isQCD());
    bool isTtbar(m_config->isTtbar());

    HadTop top_cand;  // reconstructed top candidates

    cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(m_ljets.size())+" ak8 candidates");
    for (const auto& ljet : ljets){
        top_cand.jets.clear();
        top_cand.ljet = ljet.index;

        if (isTopTagged(ljet))
            m_ttbar.push_back( top_cand );
    } // end loop over ak8 candidates

    cma::DEBUG("TTBARRECO : Ttbar built ");

    return;
}


bool ttbarReco::isTopTagged(const Ljet& ljet){
    /* Check if ljet is top tagged */
    bool istagged = (std::abs(ljet.BEST_class-4)<1e-6);
    return istagged;
}

// THE END //

