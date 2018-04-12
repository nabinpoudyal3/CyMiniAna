/*
Created:        --
Last Updated:   21 March 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Tool for building ttbar system
 - all-had:  Two top-tagged AK8
 - l+jets:   1 top-tagged AK8 + 1 leptonic top (AK4+lep+nu)
 - dilepton: 2 leptonic tops (AK4+lep+nu)
*/
#include "Analysis/CyMiniAna/interface/ttbarReco.h"


ttbarReco::ttbarReco( configuration& cmaConfig ) :
  m_config(&cmaConfig){
    m_mapContainment = m_config->mapOfPartonContainment();  // containment map (ints and strings)
    m_targetMap = m_config->mapOfTargetValues();

    m_ttbar0L = {};
    m_ttbar1L = {};
    m_ttbar2L = {};
  }

ttbarReco::~ttbarReco() {}



// all-hadronic
void ttbarReco::execute(std::vector<Ljet>& ljets){
    /* Build top quarks system 
       - 2 AK8 jets
         > highest BEST_t scores
    */
    m_ttbar0L = {};

    cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(ljets.size())+" ak8 candidates");

    if (ljets.size()>2){
        // pick the 2 AK8 candidates with highest BEST_t scores
        unsigned int firstAK8(0);
        unsigned int secondAK8(1);

        float firstBEST_t(-999.);
        float secondBEST_t(-999.);

        for (const auto& ljet : ljets){
            if (ljet.BEST_t>firstBEST_t){
                firstAK8    = ljet.index;
                firstBEST_t = ljet.BEST_t;
            }
            else if (ljet.BEST_t>secondBEST_t){
                secondAK8    = ljet.index;
                secondBEST_t = ljet.BEST_t;
            }
        } // end loop over ak8 candidates

        // sort by pT
        Ljet first  = ljets.at(firstAK8);
        Ljet second = ljets.at(secondAK8);

        if (first.p4.Pt()>second.p4.Pt()){
            m_ttbar0L.ljets.push_back( first );
            m_ttbar0L.ljets.push_back( second );
        }
        else{
            m_ttbar0L.ljets.push_back( second );
            m_ttbar0L.ljets.push_back( first );
        }
    }
    else
        m_ttbar0L.ljets = ljets;

    return;
}


// single lepton
void ttbarReco::execute(std::vector<Lepton>& leptons, std::vector<Jet>& jets, std::vector<Ljet>& ljets){
    /* Build top quarks system 
       - lepton
       - AK4 near lepton (2D cut)
         > highest pT
       - AK8 away from lepton
         > Most 'top-like' = highest BEST_t score
    */
    m_ttbar1L = {};

    // Setup lepton (only 1 in the single lepton analysis)
    cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(leptons.size())+" leptons");
    Lepton lep;
    if (leptons.size()>0)
        lep = leptons.at(0);
    else
        lep.p4.SetPtEtaPhiE(0,0,0,0);   // this event will fail the selection anway, use dummy value


    // Get Jets (only if lepton exists)
    int ak4candidate(-1);
    int ak8candidate(-1);
    if (lep.p4.Pt()>10){
        // -- Setup AK4 jet (highest pT jet within [0.3,pi/2] of lepton)
        cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(jets.size())+" ak4 candidates");
        float ak4_pt(0);

        for (auto& jet : jets){
            float dr    = jet.p4.DeltaR(lep.p4);
            float ptrel = jet.p4.Perp( lep.p4.Vect() );
            //if (0.3<dr && dr<M_HALF_PI && jet.p4.Pt() > ak4_pt){
            if ( (dr>0.4 || ptrel>25) && jet.p4.Pt() > ak4_pt){        // 2D cut instead of DeltaR window
                ak4_pt = jet.p4.Pt();
                ak4candidate = jet.index;
            }
        } // end loop over ak4 candidates


        // -- Setup AK8 jet (highest BEST_t jet farther away than pi/2 from lepton)
        cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(ljets.size())+" ak8 candidates");
        float BEST_t(-999.);            // between 0 and 1 for real jets

        for (const auto& ljet : ljets){
            float dr = ljet.p4.DeltaR(lep.p4);
            if (dr>M_HALF_PI && ljet.BEST_t>BEST_t){
                BEST_t = ljet.BEST_t;
                ak8candidate = ljet.index;
            }
        } // end loop over ak8 candidates
    } // end if lepton has sufficient pT


    // Define objects in the struct -- no neutrino yet
    m_ttbar1L.lepton = lep;

    Jet dummy_jet;
    dummy_jet.isGood = false;
    if (ak4candidate>=0) m_ttbar1L.jet = jets.at(ak4candidate);
    else m_ttbar1L.jet = dummy_jet;

    Ljet dummy_ljet;
    dummy_ljet.isGood = false;
    if (ak8candidate>=0) m_ttbar1L.ljet = ljets.at(ak8candidate);
    else m_ttbar1L.ljet = dummy_ljet;

    cma::DEBUG("TTBARRECO : Ttbar built ");

    return;
}


// dilepton
void ttbarReco::execute(std::vector<Electron>& electrons, std::vector<Muon>& muons, std::vector<Jet>& jets){
    return;
}

// THE END //
