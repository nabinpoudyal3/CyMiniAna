#ifndef EVENTSELECTION_H_
#define EVENTSELECTION_H_

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "diHiggs/CyMiniAna/interface/Event.h"
#include "diHiggs/CyMiniAna/interface/configuration.h"

class eventSelection{

  public:
    // constructor and destructor
    eventSelection(configuration &cmaConfig, const std::string &level="");
    virtual ~eventSelection();

    // Run once at the start of the job to setup the cuts
    virtual void initialize();                            // use cut file defined in configuration
    virtual void initialize(const std::string &cutsfile); // use any cut file desired
    virtual void identifySelection();

    // Run for every event (in every systematic) that needs saving
    virtual bool applySelection(Event &event, TH1D &cutflow, TH1D &cutflow_unweighted);

    // External access to information in this class
    virtual void getCutNames();
    virtual std::vector<std::string> cutNames();     // Return a vector of the cut names (for labeling bins in cutflow histograms)
    virtual unsigned int numberOfCuts();             // Return the number of cuts (for binning cutflow histograms)

  protected:

    // struct for holding information on a 'cut'
    //  ideally this could be extended so that cuts are parsed & written by code, not humans!
    struct Cut{
        std::string name;       // name of cut
        std::string comparison; // sign of cut (<,<=,>,>=,==,!=)
        float value;            // float value -- cast to int if needed
    };

    configuration* m_config;

    // cut information
    std::string m_level;
    std::string m_selection;
    std::string m_cutsfile;
    unsigned int m_numberOfCuts;
    std::vector<std::string> m_cutflowNames;
    std::vector<Cut> m_cuts;

    // booleans for each selection
    bool m_exampleSelection;
    bool m_example2Selection;
};

#endif


