#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

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

#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/physicsObjects.h"

class eventSelection{

  public:
    // constructor and destructor
    eventSelection(configuration &cmaConfig, const std::string &level="");
    virtual ~eventSelection();

    // Run once at the start of the job to setup the cuts
    virtual void initialize();                            // use cut file defined in configuration
    virtual void initialize(const std::string &cutsfile); // use any cut file desired
    virtual void identifySelection();

    // Run for every tree (before the event loop)
    void setCutflowHistograms(TH1D &cutflow, TH1D &cutflow_unweighted);

    // Run for every event (in every systematic) that needs saving
    virtual bool applySelection(const Event& event);

    // Selections put into functions (easily reference them in other cuts)
    bool allHadDNNSelection(double cutflow_bin);
    bool zeroLeptonSelection(double cutflow_bin);
    bool oneLeptonSelection(double cutflow_bin);
    bool ejetsSelection(double cutflow_bin);
    bool mujetsSelection(double cutflow_bin);
    bool twoLeptonSelection(double cutflow_bin);

    // Helper functions: Provide external access to information in this class
    void fillCutflows(double cutflow_bin);                                // fill cutflow histograms
    virtual void getCutNames();
    virtual std::vector<std::string> cutNames(){ return m_cutflowNames;}  // Return a vector of the cut names 
    virtual unsigned int numberOfCuts(){ return m_numberOfCuts;}          // Return the number of cuts


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
    std::string m_level;     // useful if you want to define one 'family' of selections 
                             // e.g., ABCD method for selection='qcd'; m_level = 'A', or something 
    std::string m_selection;
    std::string m_cutsfile;
    unsigned int m_numberOfCuts;
    std::vector<std::string> m_cutflowNames;
    std::vector<Cut> m_cuts;

    // cutflow histograms
    TH1D m_cutflow;
    TH1D m_cutflow_unw;

    // booleans for each selection
    bool m_dummySelection;
    bool m_isZeroLeptonAnalysis;
    bool m_isOneLeptonAnalysis;
    bool m_isTwoLeptonAnalysis;

    bool m_allHadDNNSelection;

    // physics information
    float m_nominal_weight;
    std::vector<Ljet> m_ljets;
    std::vector<Jet> m_jets;
    std::vector<Muon> m_muons;
    std::vector<Electron> m_electrons;
    std::vector<Neutrino> m_neutrinos;
    MET m_met;
    float m_ht;
    float m_st;

    unsigned int m_NElectrons;
    unsigned int m_NMuons;
    unsigned int m_NJets;
    unsigned int m_NLjets;
};

#endif
