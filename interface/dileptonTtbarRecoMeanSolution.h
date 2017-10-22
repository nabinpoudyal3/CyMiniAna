#ifndef DILEPTONTTBARRECOMEANSOLUTION_H
#define DILEPTONTTBARRECOMEANSOLUTION_H

#include <vector>
#include <TLorentzVector.h>

#include "cms-ttbarAC/CyMiniAna/interface/tools.h"

class dileptonTtbarRecoMeanSolution {

public:

    dileptonTtbarRecoMeanSolution(const double& topm=172.5);
    ~dileptonTtbarRecoMeanSolution();  
    void add(const Top& top, const Top& topbar, const Neutrino& n, const Neutrino& nbar, const double& weight, const double& mbl_weight);
    void add(const Top& top, const Top& topbar, const Neutrino& n, const Neutrino& nbar, const double& weight);

    void getMeanVect(cmaBase& lv, const std::vector<cmaBase>& vlv, const double& mass) const;
    void getMeanSol(Top& top, Top& topbar, Neutrino& n, Neutrino& nbar) const;
    double getSumWeight() const;
    unsigned int getNsol() const;
    void clear();

private:

    std::vector<TLorentzVector> v_top_;
    std::vector<TLorentzVector> v_topbar_;
    std::vector<TLorentzVector> v_n_;
    std::vector<TLorentzVector> v_nbar_;

    std::vector<double> v_weight_;
    double sum_weight_;  
    double max_sum_weight_;

    const double mass_top_;
};

#endif