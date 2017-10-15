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
    std::map<std::string, char> isbtagged;
    float jvt;
    int true_flavor;
    double rho;      // jet energy density, 1 value per event (attaching to each jet for convenience)
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
};

struct Neutrino : cmaBase{
    // neutrino attributes
};


struct Top {
    // Define a leptonically-decaying top quark
    TLorentzVector p4;

    Lepton lepton;
    Jet jet;
    Neutrino neutrino;
    float weight;       // reconstruction weight
    float weight_ES;    // reconstruction weight (EventShape)

    void set_p4(){
        p4 = lepton.p4 + jet.p4 + neutrino.p4;
    }
};


struct DileptonReco {
    // struct of information needed for neutrino reconstruction (AMWT)
    Lepton lepton_pos;      // positively charged lepton
    Lepton lepton_neg;      // negatively charged lepton
    TVector2 met;           // MET (stored as TVector2 for convenience)
    std::vector<Jet> jets;  // two 'b'-jets in ttbar decay
};

#endif