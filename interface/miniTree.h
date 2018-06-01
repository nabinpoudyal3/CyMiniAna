#ifndef MINITREE_H
#define MINITREE_H

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <memory>
#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/eventSelection.h"
#include "Analysis/CyMiniAna/interface/configuration.h"

class miniTree {
  public:
    // Default - so root can load based on a name;
    miniTree(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~miniTree();

    // Run once at the start of the job;
    virtual void initialize(TTree * t, TFile& outputFile, const std::string directory, const int cloneFactor=0);
    virtual void createBranches();
    virtual void disableBranches();

    bool branch_exists(const std::string& br);

    // Run for every event (in every systematic) that needs saving;
    virtual void saveEvent(Event &event, const std::vector<unsigned int>& evtsel_decisions=std::vector<unsigned int>());

    // Clear stuff;
    virtual void finalize();


  protected:

    TTree * m_ttree;
    TTree * m_oldTTree;
    configuration * m_config;

    std::vector<std::string> m_selections;
    std::vector<std::string> m_listOfBranches;

    // new branches defined here 
    float m_dnn;
    std::vector<float> m_BEST_t_j;
    std::vector<unsigned int> m_passSelection;

    unsigned int m_leptop_jet;
    unsigned int m_hadtop_ljet;

    std::vector<float> m_nu_pt;
    std::vector<float> m_nu_eta;
    std::vector<float> m_nu_phi;
};

#endif
