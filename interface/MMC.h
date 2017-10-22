#ifndef MMC_H
#define MMC_H

//std lib
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"  
#include <time.h>

//root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TDirectory.h"

// local functions/classes
#include "cms-ttbarAC/CyMiniAna/interface/tools.h"
#include "cms-ttbarAC/CyMiniAna/interface/physicsObjects.h"
#include "cms-ttbarAC/CyMiniAna/interface/configuration.h"


typedef std::pair<float, float> EtaPhi;

class MMC {

  public:
    //constructor
    MMC( configuration &cmaConfig );
/*    MMC(TLorentzVector* mu1_lorentz, TLorentzVector* mu2_lorentz, TLorentzVector* b1jet_lorentz, TLorentzVector* b2jet_lorentz, 
    TLorentzVector* totjets_lorentz,TLorentzVector* met_lorentz, TLorentzVector* nu1_lorentz, TLorentzVector* nu2_lorentz,
    TLorentzVector* b_genp_lorentz, TLorentzVector* bbar_genp_lorentz, TLorentzVector* h2tohh_lorentz, int onshellMarker, bool simulation, bool PUsample_,
    int ievent, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist,
    int iterations, std::string RefPDFfile, bool useMET, int bjetrescaleAlgo, int metcorrection, int verbose_=0);
*/
    ~MMC();

    void initialize();
    bool execute( unsigned int entry, DileptonReco& dilep );
    TH1F getMMCh2();
    TH1F getMMCh2weight1();
    TH1F getMMCh2weight4();


  private:

    TVector2 metCorrection();
    bool bjetsCorrection();
    void initTree(TTree* mmctree);
    float genEtaGuass(float mean, float rms);
    float genPhiFlat();
    EtaPhi generatenu1_etaphi();
    float nu1pt_onshellW(EtaPhi nu1_etaphi, float wMass);

    bool getNuOffShellW(const TVector2& met, const Neutrino& nu_onshellW, 
                        Neutrino& nu_offshellW, const int control, const float hMass);

    bool checkSolution(TLorentzVector* jetslorentz,
                        TLorentzVector* mu1lorentz,
                        TLorentzVector* mu2lorentz,
                        TLorentzVector* nu1lorentz, int control, float hMass); 
    bool cutsCheck();
    void assignMuons(int control);  

    float onshellWMassRandomWalk(float x0, float step, float random);
    float onshellWMassRandomWalk(float x0, float step, float random, bool hist_);
    float onshellWMassPDF(float wmass);

    TH1F* getOnShellWMassPDF();
    TH1F* getOffShellWMassPDF();
    TH2F* getOnOffShellWMassPDF();
    TH1F* getOnShellNuPtPDF();
    TH1F* getBjetRescalec1PDF();
    TH1F* getBjetRescalec2PDF();
    TH2F* getBjetRescalec1c2PDF();

    float weightfromhist(TH1F* pdf, float x); 
    float weightfromhist(TH2F* pdf, float x, float y, bool whole=true); 
    float weightfromonshellnupt(float nupt); 


    configuration* m_config;
    TH1F m_MMC_h2Mass;
    TH1F m_MMC_h2Massweight1;
    TH1F m_MMC_h2Massweight4;
    TFile* m_file;
    bool m_isMC;
    TRandom3* m_generator;
    double m_pi;

    bool m_weightfromonshellnupt_func;
    bool m_weightfromonshellnupt_hist;
    bool m_weightfromonoffshellWmass_hist;
    bool m_weightfrombjetrescalec1c2_hist;

    int m_onshellMarker;
    bool m_puSample;
    int m_iterations;
    int m_metcorrection;
    int m_bjetrescale;
    float m_b1rescalefactor;
    float m_b2rescalefactor;
    float m_rescalec1;
    float m_rescalec2;
    float m_met_sigma;
    bool m_useMET;
    float m_WMass;
    float hmass_gen;

    Lepton m_lep_p;
    Lepton m_lep_n;
    TVector2 m_met;
    std::vector<Jet> m_jets;
    TLorentzVector m_bbar;
    TLorentzVector m_totjets;

    TLorentzVector* mmc_mu1_lorentz;
    TLorentzVector* mmc_mu2_lorentz;
    TLorentzVector* mmc_bjets_lorentz;
    TLorentzVector* mmc_b1jet_lorentz;
    TLorentzVector* mmc_b2jet_lorentz;
    TLorentzVector* mmc_totjets_lorentz;
    TLorentzVector* nu1_lorentz_true;
    TLorentzVector* nu2_lorentz_true;
    TLorentzVector* onshellW_lorentz_true;
    TLorentzVector* offshellW_lorentz_true;
    TVector2* mmcmet_vec2;
    TLorentzVector* b1_lorentz;
    TLorentzVector* b2_lorentz;
    TLorentzVector* htoWW_lorentz_true;
    TLorentzVector* htoBB_lorentz_true;
    TLorentzVector* h2tohh_lorentz_true;

    Lepton m_mu_onshellW;
    Lepton m_mu_offshellW;
    TLorentzVector* jets_lorentz;
    Neutrino nu_onshellW;
    Neutrino nu_offshellW;
    TLorentzVector* offshellW_lorentz;
    TLorentzVector* onshellW_lorentz;
    TLorentzVector* htoWW_lorentz;
    TLorentzVector* htoBB_lorentz;
    TLorentzVector* h2tohh_lorentz;

    TLorentzVector ideal_met_lorentz;
    TLorentzVector h2tohh_expect_lorentz;

    TH2F* m_h_onoffshellWmass;
    TH1F* m_h_onshellnupt;
    TH1F* m_h_bjetrescalec1;
    TH1F* m_h_bjetrescalec2;

    float b1rescalefactor_true;
    float b2rescalefactor_true;
    float rescalec1_true;
    float rescalec2_true;

    float eta_nuoffshellW_true;
    float phi_nuoffshellW_true;
    float pt_nuoffshellW_true;
    float px_nuoffshellW_true;
    float py_nuoffshellW_true;
    float eta_nuonshellW_true;
    float phi_nuonshellW_true;
    float pt_nuonshellW_true;
    float px_nuonshellW_true;
    float py_nuonshellW_true;
    float mass_offshellW_true;
    float mass_onshellW_true;
    float mass_htoWW_true;
    float pt_h2tohh_true;
    float mass_h2tohh_true;
    float mass_h2_expect;
};

#endif
