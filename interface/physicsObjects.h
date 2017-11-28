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
struct CmaBase {
    TLorentzVector p4;
};


// Struct for particle flow jets
// -- common to all types of jets
struct Jet : CmaBase{
    float cMVAv2;
    float CSVv2;
    float CvsL;
    float CvsB;
    std::map<std::string, bool> isbtagged;

    float partonFlavour;
    float hadronFlavour;
    float neutralMultiplicity;
    float neutralHadronEnergyFrac;
    float neutralEmEnergyFrac;
    float chargedHadronEnergyFrac;
    float chargedEmEnergyFrac;
    float chargedMultiplicity;

    float jecFactor0;
    float jetArea;
    float ptResolution;
    float smearedPt;
    std::vector<int> keys;

    float charge;
    int index;    // index in vector of jets
};

struct Ljet : Jet{
    // extra ljet attributes
    int isGood;
    float tau1_CHS;
    float tau2_CHS;
    float tau3_CHS;
    float tau21_CHS;
    float tau32_CHS;
    float softDropMass_CHS;
    float vSubjetIndex0;
    float vSubjetIndex1;
};

struct Tjet : Jet{
    // Track Jets not used in CMS -- here as a placeholder
    int numConstituents;
};



// Extra lepton attributes
struct Lepton : CmaBase{
    // common to electrons and muons
    int charge;
    bool isElectron;
    bool isMuon;
    int index;       // index in vector of leptons

    float key;
    float miniIso;
    float loose;
    float medium;
    float tight;
};

struct Electron : Lepton{
    // extra electron attributes
    Lepton::isElectron = true;
    Lepton::isMuon     = false;

    float iso03;
    float iso03db;
    float SCEta;
    float SCPhi;
    float vidVeto;
    float vidLoose;
    float vidMedium;
    float vidTight;
    float vidHEEP;
    float vidVetonoiso;
    float vidLoosenoiso;
    float vidMediumnoiso;
    float vidTightnoiso;
    float vidHEEPnoiso;
    float vidMvaGPvalue;
    float vidMvaGPcateg;
    float vidMvaHZZvalue;
    float vidMvaHZZcateg;
    int veto_NoIsoID;
    int loose_NoIsoID;
    int medium_NoIsoID;
    int tight_NoIsoID;
    int isoVeto;
    int isoLoose;
    int isoMedium;
    int isoTight;
    int vetoID;
};
struct Muon : Lepton{
    // extra muon attributes
    Lepton::isElectron = false;
    Lepton::isMuon     = true;

    float iso04;
    float soft;
    float medium2016;
    float hightPt;
};

struct Neutrino : CmaBase{
    // extra neutrino attributes
};



// Top Quarks
struct Top : CmaBase{
    // only used for truth tops (parton-level) right now (no extra information, yet)
};

struct LepTop : CmaBase{
    // Leptonically-decaying top quark aspects
    Lepton lepton;
    Jet jet;
    Neutrino neutrino;
    float weight;        // reconstruction weight
    float weight_ES;     // reconstruction weight (EventShape)
    float weight_tt;     // reconstruction weight (ttbar XS method)

    void set_p4(){
        p4 = lepton.p4 + jet.p4 + neutrino.p4;
    }
};

struct HadTop : CmaBase{
    // Resolved hadronically-decaying top quark
    std::vector<Jet> jets; // contains all associated jets
    Jet bjet; // if b-jet is identified
    Jet q1;   // if quark from W is identified
    Jet q2;   // if quark from W is identified

    void set_p4(){
        p4 = TLorentzVector(0,0,0,0);
        for (const auto& jet : jets)
            p4 += jet.p4;
    }
};

struct BoostedTop : CmaBase{
    // Boosted hadronically-decaying top quark
    Ljet jet;
    float dnn;

    void set_p4(){
        p4 = jet.p4;
    }
};



// ************************************ //
// ttbar System

//struct TtbarOneLepton{};
//struct TtbarAllHad{};

struct DileptonReco {
    // struct of information needed for neutrino reconstruction (AMWT)
    Lepton lepton_pos;// positively charged lepton
    Lepton lepton_neg;// negatively charged lepton
    TVector2 met;// MET (stored as TVector2 for convenience)
    std::vector<Jet> jets;// jets in event
    std::vector<Jet> bjets;// two 'b'-jets in ttbar decay
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

    LepTop top;
    LepTop topBar;
    TLorentzVector ttbar;

    double recMtop;     // top mass used in the reconstruction (needed if doing massLoop)
    double weight;      // weight of the solution
    int ntags;          // number of b-tags

    double dR;
    double dN;
    double x1;
    double x2;
    double mtt;         // invariant mass of the ttbar system
    bool isNoSmearSol;  // if smearing, this bool identifies solutions found w/o smearing

    // Enumeration for all defined weights of a kinematic reconstruction solution
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

