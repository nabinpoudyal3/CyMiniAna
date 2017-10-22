/**
Original Author Lars Sonnenschein <sonne@lpnhep.in2p3.fr>

Minor modifcations by Dan Boline <ddboline@fnal.gov>

Further minor modifications by Aram Avetisyan <aram.avetisyan@cern.ch>
Debug by Thomas Speer.


Source: https://github.com/cms-ljmet/topDileptonMassAnalysis
 */
#ifndef MASSSOLVER_CLASS_H
#define MASSSOLVER_CLASS_H


#include "TLorentzVector.h"
#include "TVector2.h"
#include "TMath.h"
#include <iostream>
#include <vector>

#include "cms-ttbarAC/CyMiniAna/interface/physicsObjects.h"
#include "cms-ttbarAC/CyMiniAna/interface/tools.h"
#include "cms-ttbarAC/CyMiniAna/interface/configuration.h"



class MassSolver {
  public:
  
    MassSolver(configuration* config);
    ~MassSolver();

    bool solve( const TVector2& met, const Jet& bq1, const Jet& bq2,
                Lepton& lep1, Lepton& lep2, double mt, 
                std::vector<Neutrino>& nu1, std::vector<Neutrino>& nu2 );

    void solve(const TVector2& met, const Jet& bq1, const Jet& bq2,
               const Lepton& lep1, const Lepton& lep2, double mt,
               std::vector<double> *pnux, std::vector<double> *pnuy, std::vector<double> *pnuz,
               std::vector<double> *pnubx, std::vector<double> *pnuby, std::vector<double> *pnubz,
               std::vector<double> *cd_diff, int& cubic_single_root_cmplx);


    void quartic(std::vector<double> const & poly, std::vector<double> *pnuy, int& cubic_single_root_cmplx);

    void cubic(std::vector<double> const & poly, std::vector<double> *pnuy);

    void quadratic(std::vector<double> const & poly, std::vector<double> *pnuy);

    int algebraic_pz(double* b, double* lp, double mt, double mb, double mlp, 
                     double pnux, double pnuy, double* pnuz);

    double evalterm1(std::vector<double> *a1, double pnux, double pnuy);

    double evalterm2(std::vector<double> *a2, double pnux, double pnuy);
  
  private: 
  
    configuration *m_config;

    static constexpr double epsilon = 1e-6; //numerical precision
    bool flagProb;
    bool debug;

    double m_topQuarkMass;
    double m_bQuarkMass;
    double m_WMass;

    float m_threshold;
    double m_Pi;

    inline double sqr(double const & x) const {return x*x;}

    inline double sign(double const & a) const {return (a < 0) ? -1 : (a > 0) ? 1 : 0;}

    inline double quad(double const & x) const {return x*x*x*x;}
}; 

#endif
