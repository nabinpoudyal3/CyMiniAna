#ifndef DILEPTONTTBARRECOUTILS_H
#define DILEPTONTTBARRECOUTILS_H

#include <vector>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TF1.h>
#include <TF2.h>
#include <TRandom3.h>

#include "cms-ttbarAC/CyMiniAna/interface/physicsObjects.h"
#include "cms-ttbarAC/CyMiniAna/interface/configuration.h"
#include "cms-ttbarAC/CyMiniAna/interface/tools.h"


class dileptonTtbarRecoUtils{

  public:

    // constructor with arguments for masses [GeV]
    dileptonTtbarRecoUtils(const std::map<std::string,double>& truth_masses={}, const double& beamEnergy=13000.);
    ~dileptonTtbarRecoUtils();
    void fDelete() const;

    void setTopMass( double topmass );
    void setAntiTopMass( double antitopmass );
    void setWplus( double wplus );
    void setWminus( double wminus );

    void execute(const DileptonReco& ttSystem);
    void execute();
    int getNsol() const;
    std::vector<TtbarDilepton> getTtSol() const;

    /* Assign truth information */
    void setTruthInfo(const Top& top, const Top& antiTop,
                      const Neutrino& neutrino, const Neutrino& antiNeutrino);

    void sortBy(std::string ch);  // re-arrange the vector based on a different metric

    void angle_rot(const double& alpha, const double& e, const CmaBase& inJet, CmaBase& jet_sm) const;


  private:

    void setdR();
    void setdN();

    void swapTopSol(TtbarDilepton& sol1, TtbarDilepton& sol2) const;
    void sortTopSol(std::vector<TtbarDilepton>& v) const;

    void doAll();
    void topRec(const double& px_neutrino);

    void findCoeff(double* const koeficienty);

    void quartic_equation(const double& h0, const double& h1, const double& h2, const double& h3, const double& h4, std::vector<double>& v) const;
    void cubic_equation(const double& a, const double& b, const double& c, const double& d, std::vector<double> &v) const;
    void quadratic_equation(const double& a, const double& b, const double& c, std::vector<double>& v) const;
    void linear_equation(const double& a, const double& b, std::vector<double>& v) const;

    double landau2D(const double& x, const double& y) const;

    int sign(const long double& ld) const;
    double sqr(const double& x) const;
    void swap(double& realone, double& realtwo) const;

    // Member variables
    double m_mass_top=172.5;
    double m_mass_topbar=172.5;
    double m_mass_b=4.8;
    double m_mass_bbar=4.8;
    double m_mass_Wp=80.4;
    double m_mass_Wm=80.4;
    double m_mass_al=0;
    double m_mass_l=0;
    double m_mass_v=0;
    double m_mass_av=0;

    int m_NSol;
    double coeffs_[5];
    std::vector<double> vect_pxv_;
    std::vector<TtbarDilepton> m_ttSol;
    TRandom3* m_r3;                /// Random number generation

    Lepton m_al;
    Lepton m_l;
    Jet m_b;
    Jet m_bbar;
    LepTop m_top;
    LepTop m_topbar;
    Neutrino m_neutrino;
    Neutrino m_neutrinobar;
    TLorentzVector m_w;
    TLorentzVector m_wbar;
    TLorentzVector m_tt;

    Top m_true_top;
    Top m_true_topBar;
    Neutrino m_true_neutrino;
    Neutrino m_true_neutrinoBar;

    double m_px_miss;
    double m_py_miss;

    double m_mt;
    double m_mtbar;
    double m_mb;
    double m_mbbar;
    double m_mw;
    double m_mwbar;
    double m_ml;
    double m_mal;
    double m_mv;
    double m_mav;

    float m_beamEnergy;

    double a1_,a2_,a3_,a4_;
    double b1_,b2_,b3_,b4_;
    double c22_,c21_,c20_,c11_,c10_,c00_;
    double d22_,d21_,d20_,d11_,d10_,d00_;
    double d0_,d1_,d2_;
    double c0_,c1_,c2_;
};

#endif

