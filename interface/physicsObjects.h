#ifndef PHYSICSOBJECTS_H_
#define PHYSICSOBJECTS_H_

/* 
   Physics objects to be used in analyses
   This structure allows the Event class
   and other classes to access these objects
   without circular inclusion (which breaks!)
*/
#include "TLorentzVector.h"
#include <map>
#include <string>

// base object (consistent reference to TLorentzVector)
struct cmaBase {
    TLorentzVector p4;
};

struct Tjet : cmaBase{
    // extra track jet attributes
    float mv2c10;
    float mv2c20;
    std::map<std::string, char> isbtagged;
    float jvt;
    double charge;
    int numConstituents;
    int true_flavor;
};

struct Jet : cmaBase{
    // extra jet attributes
    float cMVAv2;
    float mv2c10;
    float mv2c20;
    std::map<std::string, bool> isbtagged;
    float jvt;
    int true_flavor;
    double rho;      // jet energy density, 1 value per event (attaching to each jet for convenience)
    int index;       // index in vector of jets
};

struct Ljet : cmaBase{
    // extra ljet attributes
    int isGood;
    float Split23;
    float tau1_wta;
    float tau2_wta;
    float tau3_wta;
    float tau21_wta;
    float tau32_wta;
    double charge;
};

struct Lepton : cmaBase{
    // extra lepton attributes
    int charge;
    bool isElectron;
    bool isMuon;
    float Iso;
    int index;       // index in vector of leptons
};

struct Neutrino : cmaBase{
    // neutrino attributes
};


struct Top {
    // Define a top quark
    TLorentzVector p4;

    // leptonically-decaying top quark aspects
    Lepton lepton;
    Jet jet;
    Neutrino neutrino;
    float weight;       // reconstruction weight
    float weight_ES;    // reconstruction weight (EventShape)
    float weight_tt;    // reconstruction weight (ttbar XS method)

    void set_p4_lep(){
        p4 = lepton.p4 + jet.p4 + neutrino.p4;
    }

    // Resolved hadronically-decaying top quark
    std::vector<Jet> jets;  // contains all associated jets
    Jet bjet;               // if b-jet is identified
    Jet q1;                 // if quark from W is identified
    Jet q2;                 // if quark from W is identified

    void set_p4_had(){
        p4 = TLorentzVector(0,0,0,0);
        for (const auto& jet : jets)
            p4 += jet.p4;
    }
};

// ttbar system
struct DileptonReco {
    // struct of information needed for neutrino reconstruction (AMWT)
    Lepton lepton_pos;      // positively charged lepton
    Lepton lepton_neg;      // negatively charged lepton
    TVector2 met;           // MET (stored as TVector2 for convenience)
    std::vector<Jet> jets;  // jets in event
    std::vector<Jet> bjets; // two 'b'-jets in ttbar decay
    Jet bJet;
    Jet bbarJet;
    size_t bJet_index, bbarJet_index;
};

struct TtbarDilepton : DileptonReco {
    /* dilepton ttbar system for dileptonTtbarReco setup
       formerly: Struct_KinematicReconstruction
       incorporated some attributes of the "TopSolution" struct
    */
    Neutrino neutrino;
    Neutrino neutrinoBar;

    TLorentzVector Wplus;
    TLorentzVector Wminus;

    Top top;
    Top topBar;
    TLorentzVector ttbar;

    double recMtop;        // top mass used in the reconstruction (needed if doing massLoop)
    double weight;         // weight of the solution
    int ntags;             // number of b-tags

    double dR;
    double dN;
    double x1;
    double x2;
    double mtt;            // invariant mass of the ttbar system
    bool isNoSmearSol;     // if smearing, this bool identifies solutions found w/o smearing

    /// Enumeration for all defined weights of a kinematic reconstruction solution
    enum WeightType{defaultForMethod, neutrinoEnergy, averagedSumSmearings_mlb, undefinedWeight};
    std::map<WeightType,double> mapOfWeights;

    void set_ttbar(){
        ttbar = top.p4 + topBar.p4;
    }
    void set_W(){
        Wplus  = lepton_pos.p4 + neutrino.p4;
        Wminus = lepton_neg.p4 + neutrinoBar.p4;
    }
};

#endif

