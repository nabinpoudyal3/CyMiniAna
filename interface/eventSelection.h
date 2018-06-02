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
    virtual void initialize(const std::string& selection, const std::string& cutsfile);
    virtual void initialize(const std::string &cutsfile);

    virtual void finalize();

    virtual void identifySelection();

    // Run for every tree (before the event loop)
    void setCutflowHistograms(TFile& outputFile);

    // Run for every event (in every systematic) that needs saving
    virtual void setObjects(const Event& event);
    virtual bool applySelection();

    // -- Selections put into functions (easily reference them in other cuts)
    // DNN selection
    bool allHadDNNSelection(float& cutflow_bin);

    // Allhadronic selection
    bool zeroLeptonSelection(float& cutflow_bin);

    // Single lepton selections
    bool oneLeptonSelection(float& cutflow_bin, std::string selection="");
    bool ejetsSelection( float& cutflow_bin, std::string selection="");
    bool mujetsSelection(float& cutflow_bin, std::string selection="");
    bool cwolaSelection( float& cutflow_bin);
    bool afbSelection(   float& cutflow_bin);

    // Dilepton selections
    bool twoLeptonSelection(float& cutflow_bin);

    // Helper functions: Provide external access to information in this class
    void fillCutflows(float& cutflow_bin);                                // fill cutflow histograms
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
    TH1D* m_cutflow;
    TH1D* m_cutflow_unw;

    // booleans for each selection
    bool m_dummySelection;
    bool m_isZeroLeptonAnalysis;
    bool m_isOneLeptonAnalysis;
    bool m_isTwoLeptonAnalysis;

    bool m_isCWoLaAnalysis;
    bool m_allHadDNNSelection;
    bool m_isAFBAnalysis;

    // physics information
    bool m_valid;
    float m_nominal_weight;
    std::vector<Ljet> m_ljets;
    std::vector<Jet> m_jets;
    std::vector<Muon> m_muons;
    std::vector<Electron> m_electrons;
    std::vector<Lepton> m_leptons;
    std::vector<Neutrino> m_neutrinos;
    MET m_met;
    float m_ht;
    float m_st;

    std::vector<std::string> m_zeroLeptonTriggers;
    std::vector<std::string> m_ejetsTriggers;
    std::vector<std::string> m_mujetsTriggers;
    std::vector<std::string> m_dileptonTriggers;

    std::map<std::string,unsigned int> m_triggers;
    std::map<std::string,unsigned int> m_filters;

    unsigned int m_Nbtags;
    unsigned int m_NLeptons;
    unsigned int m_NElectrons;
    unsigned int m_NMuons;
    unsigned int m_NJets;
    unsigned int m_NLjets;

    Ttbar0L m_ttbar0L;
    Ttbar1L m_ttbar1L;
    Ttbar2L m_ttbar2L;
};

#endif
