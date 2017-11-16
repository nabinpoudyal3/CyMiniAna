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

#include "cms-ttbarAC/CyMiniAna/interface/Event.h"
#include "cms-ttbarAC/CyMiniAna/interface/eventSelection.h"
#include "cms-ttbarAC/CyMiniAna/interface/configuration.h"

class miniTree {
  public:
    // Default - so root can load based on a name;
    miniTree(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~miniTree();

    // Run once at the start of the job;
    virtual void initialize(TTree * t, TFile& outputFile);

    // Run for every event (in every systematic) that needs saving;
    virtual void saveEvent(Event &event);

    // Clear stuff;
    virtual void finalize();


  protected:

    TTree * m_ttree;
    TTree * m_oldTTree;
    configuration * m_config;

    // new branches defined here 
    float m_dnn;
};

#endif
