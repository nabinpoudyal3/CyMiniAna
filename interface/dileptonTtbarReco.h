#ifndef DILEPTONTTBARRECO_H
#define DILEPTONTTBARRECO_H

#include <vector>
#include <map>
#include <TLorentzVector.h>
#include <TMath.h>

#include "diHiggs/CyMiniAna/interface/tools.h"
#include "diHiggs/CyMiniAna/interface/physicsObjects.h"
#include "diHiggs/CyMiniAna/interface/configuration.h"
#include "diHiggs/CyMiniAna/interface/dileptonTtbarRecoSolution.h"
#include "diHiggs/CyMiniAna/interface/dileptonTtbarRecoUtils.h"



class dileptonTtbarReco{

  public:

    dileptonTtbarReco(configuration& cmaConfig,   const configuration::Era era, 
                      const int minNumberOfBtags, const bool preferBtags, 
                      const bool massLoop=false);
    ~dileptonTtbarReco();

    int getNSol() const;
    TtbarDilepton getSol() const;
    std::vector<TtbarDilepton> getSols() const;

    void loadData();

    /// Retrieve all solutions valid for setup of kinematic reconstruction
    std::map<std::string,Top> execute(const DileptonReco& ttSystem);
    dileptonTtbarRecoSolution buildTtbar(const TtbarDilepton& ttSystem);


    /* NOT USED */
//     void kinReco(const Lepton& leptonMinus,   const Lepton& leptonPlus, 
//                  const std::vector<Jet> jets, const std::vector<double> btags, const TVector2 met);
//     void kinRecoMassLoop(const Lepton& leptonMinus,   const Lepton& leptonPlus, 
//                          const std::vector<Jet> jets, const std::vector<double> btags, const TVector2 met);

  private:

    /// Calculate solution for specific lepton, antilepton and pair of jets
    std::vector<TtbarDilepton> getSolutions(const TtbarDilepton& singleTtbarSystem, const int numberOfBtags);

    /// Calculate solution using smearing
    bool getSmearedSolutions(std::vector<TtbarDilepton>&, std::vector<double>&, const DileptonReco&) const;
    void averageSmearedSolutions(const std::vector<TtbarDilepton>& smeared_solutions,
                                 const std::vector<double>& smeared_weights,
                                 TtbarDilepton& solution);

    // FIXME: temporary helper variables for cleanup
    void setSolutions();
    void setSolutions(std::vector<TtbarDilepton> sols);


    // Member variables
    configuration *m_config;        /// CyMiniAna configuration
    TRandom3* m_r3;                 /// Random number generation
    const configuration::Era m_era; /// Analysis era
    const int m_minNumberOfBtags;   /// Minimum number of b-tags required for solutions (0, 1, 2)
    const bool m_preferBtags;       /// Prefer solutions with b-tags (2 tags if existing, else 1 tag if existing, else 0 tags)
    const bool m_massLoop;          /// Whether to run mass loop for top mass, instead of smearings according to uncertainties
    std::string m_btag_wp;
    int m_rangeLow, m_rangeHigh;

    int m_NSol;
    TtbarDilepton m_sol;
    std::vector<TtbarDilepton> m_sols;

    double m_topMass;

    TH1* m_h_wmass;        // W mass

    TH1* m_h_jetAngleRes;  // jet resolution
    TH1* m_h_jetEres;

    TH1* m_h_lepAngleRes;  // lepton resolution
    TH1* m_h_lepEres;

    TH1* m_h_mbl_w;        // mbl

    /* NOT USED */
//     void inputNoJetMerging(std::vector<int>& b1_id, 
//                            std::vector<int>& b2_id, 
//                            std::vector<int>& nb_tag,
//                            const std::vector<double>& btags) const;
//     void inputNoJetMerging(std::vector<int>& b1_id, 
//                            std::vector<int>& b2_id, 
//                            std::vector<int>& nb_tag,
//                            const std::vector<double>& btags,
//                            std::vector<double> btag_ww) const;
};

#endif


