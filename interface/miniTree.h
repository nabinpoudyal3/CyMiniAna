#ifndef MINITREE_H_
#define MINITREE_H_

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

#include "diHiggs/CyMiniAna/interface/Event.h"
#include "diHiggs/CyMiniAna/interface/eventSelection.h"
#include "diHiggs/CyMiniAna/interface/configuration.h"

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
    double m_dnn;
    double m_hme;
};

#endif
