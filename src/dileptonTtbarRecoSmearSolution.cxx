/*
Created:        --
Last Updated:   15 October  2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

 Methods for dileptonTtbarReco
  Class for dilepton ttbar reconstruction.
  Imported from
   https://gitlab.cern.ch/cms-desy-top/TopAnalysis/blob/master/
           Configuration/analysis/common/src/KinematicReconstruction_MeanSol.cc
  on 15 October 2017

*/
#include "cms-ttbarAC/CyMiniAna/interface/dileptonTtbarRecoMeanSolution.h"



dileptonTtbarRecoMeanSolution::dileptonTtbarRecoMeanSolution(const double& topm):
  sum_weight_(0),
  max_sum_weight_(0),
  m_mass_top(topm){}


dileptonTtbarRecoMeanSolution::~dileptonTtbarRecoMeanSolution() {}


void dileptonTtbarRecoMeanSolution::clear() {
    m_tops.clear();
    m_topbars.clear();
    m_ns.clear();
    m_nbars.clear();
    m_weight.clear();
    m_sum_weight     = 0;
    m_max_sum_weight = 0;

    return;
}



void dileptonTtbarRecoMeanSolution::Add(const std::vector<ttbarDilepton> ttSolution, 
                                        const double& weight, const double& mbl_weight){
    m_tops.push_back(ttSolution.top);
    m_topbars.push_back(ttSolution.topBar);
    m_ns.push_back(ttSolution.neutrino);
    m_nbars.push_back(ttSolution.neutrinoBar);
    m_weight.push_back(weight);

    m_sum_weight     += weight;
    m_max_sum_weight += mbl_weight;

    return;
}


void dileptonTtbarRecoMeanSolution::Add(const std::vector<ttbarDilepton> ttSolution, 
                                        const double& weight){
    m_tops.push_back(ttSolution.top);
    m_topbars.push_back(ttSolution.topBar);
    m_ns.push_back(ttSolution.neutrino);
    m_nbars.push_back(ttSolution.neutrinoBar);
    m_weight.push_back(weight);

    m_sum_weight     += weight;
    m_max_sum_weight += weight;

    return;
}


void dileptonTtbarRecoMeanSolution::getMeanVect(cmaBase& lv, 
                                                const std::vector<cmaBase>& vlv, 
                                                const double& mass) const {
    double px_sum(0);
    double py_sum(0);
    double pz_sum(0);
    double px(0);
    double py(0);
    double pz(0);

    for(unsigned int i=0,size=vlv.size();i<size;++i){
          px_sum += m_weight.at(i) * vlv.at(i).p4.Px();
          py_sum += m_weight.at(i) * vlv.at(i).p4.Py();
          pz_sum += m_weight.at(i) * vlv.at(i).p4.Pz();
    }

    px = px_sum/m_sum_weight;
    py = py_sum/m_sum_weight;
    pz = pz_sum/m_sum_weight;

    lv.p4.SetXYZM(px,py,pz,mass);

    return;
}



void dileptonTtbarRecoMeanSolution::getMeanSol(Top& top, Top& topbar, Neutrino& n, Neutrino& nbar) const {
    /* Get the average four-vector for tops/neutrinos
       Return the top, tobar, n, nbar back to the user that called this function
    */
    getMeanVect(top,   m_tops,   m_mass_top);
    getMeanVect(topbar,m_topbars,m_mass_top);
    getMeanVect(n,   m_ns,   0);
    getMeanVect(nbar,m_nbars,0);

    return;
}


double dileptonTtbarRecoMeanSolution::getSumWeight() const {
    return m_sum_weight; // for 1 weight
}


unsigned int dileptonTtbarRecoMeanSolution::getNsol() const {
    return m_tops.size();
}

// THE END