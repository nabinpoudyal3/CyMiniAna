#ifndef AMWT_H
#define AMWT_H

// Source: https://github.com/cms-ljmet/topDileptonMassAnalysis

#include <memory>
#include <vector>
#include <sstream>
#include <iostream>
#include <math.h>
#include <random>
#include <stdint.h>
#include <TH1.h>
#include <TF2.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <LHAPDF/LHAPDF.h>
#include "LHAPDF/GridPDF.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "cms-ttbarAC/CyMiniAna/interface/physicsObjects.h"
#include "cms-ttbarAC/CyMiniAna/interface/tools.h"
#include "cms-ttbarAC/CyMiniAna/interface/configuration.h"
#include "cms-ttbarAC/CyMiniAna/interface/MassSolver.h"

class AMWT{

  public:
    AMWT(configuration &cmaConfig);
    ~AMWT();

    void initialize();
    std::map<std::string,Top> findMass(DileptonReco& dilep);

    int massPoints() const   {return m_NMassPts;}
    int massRangeMin() const {return m_rangeLow;}
    int massRangeMax() const {return m_rangeHigh;}
    int NJetSmearIterations() const {return m_NJetSmear;}
    TH1F * sampleWeightHisto() const;

  private:

    double get_dalitz_prob( const Top& top) const;
    double get_weight(const Top& t, const Top& tbar, const double top_iter) const;
    double get_weight_ES(Neutrino& nu_1, Neutrino& nu_2) const;

    // smear Jet MET
    void smear_JetMET(const std::vector<Jet>& orig_jets, const TVector2& orig_met,
        std::vector<Jet>& smear_jets, TVector2& smear_met,
        const TLorentzVector& lep_sum, bool smear) const;
    // -- original function:
    void smear_JetMET(const std::vector<Jet>& orig_jets, const TVector2& orig_met,
        std::vector<Jet>& smear_jets, TVector2& smear_met,
        const TLorentzVector& lep_sum, int old) const;
    // -- updated function: 
    void smear_JetMET(const std::vector<Jet>& orig_jets, const TVector2& orig_met,
        std::vector<Jet>& smear_jets, TVector2& smear_met, const TLorentzVector& lep_sum) const;

    configuration* m_config;
    float m_sqrt_s;

    float m_topQuarkMass;
    float m_bQuarkMass;
    float m_WMass;
    bool m_doGaussian;

    int m_NMassPts;
    int m_rangeLow, m_rangeHigh;
    int m_NJetSmear;
    float m_max_weight;         // max weight for reconstruction
    float m_max_weight_ES;      // max weight for reconstruction (EventShape)
    Top m_top;
    Top m_antitop;

    MassSolver* m_massSolver;
    TRandom3*   m_rand3;
    TF2* m_eventShape;

    // original resolution smearing implementation
    JetResolution* m_ptResol;
    JetResolution* m_etaResol;
    JetResolution* m_phiResol;

    // new implementation of jet resolution smearing
    LHAPDF::PDF* m_pdf;
    JME::JetResolution m_resolution_pt;
    JME::JetResolution m_resolution_phi;
    JME::JetResolutionScaleFactor m_resolution_sf;

    std::unique_ptr<JME::JetResolution> m_resolution_from_file_pt;
    std::unique_ptr<JME::JetResolution> m_resolution_from_file_phi;
    std::unique_ptr<JME::JetResolutionScaleFactor> m_scale_factor_from_file;
};

#endif
