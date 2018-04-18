#ifndef METADATATREE_H
#define METADATATREE_H

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
#include "Analysis/CyMiniAna/interface/miniTree.h"
#include "Analysis/CyMiniAna/interface/eventSelection.h"
#include "Analysis/CyMiniAna/interface/configuration.h"

class metadataTree : public miniTree {
  public:
    // Default
    metadataTree(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~metadataTree();

    // Run once at the start of the job;
    virtual void initialize(TTree * t, TFile& outputFile, const std::string directory, const bool recalculateMetadata=true);
    virtual void createBranches();

    bool branch_exists(const std::string& br);

    // Run for every event (in every systematic) that needs saving;
    virtual void saveMetaData(const Sample& smp);

  protected:

    TTree * m_ttree;
    TTree * m_oldTTree;
    configuration * m_config;

    // new branches defined here 
    std::string m_sampleName;
    float m_sumOfWeights;
    float m_xsection;
    float m_kfactor;
    unsigned int m_NEvents;
};

#endif
