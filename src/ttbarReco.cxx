/*
Created:        21 March 2018
Last Updated:   16 May   2018

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

    if (m_config->isTwoLeptonAnalysis())
        m_dileptonTtbar = new dileptonTtbarReco(cmaConfig, configuration::run2_13tev_2016_25ns, 2, true);
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
void ttbarReco::execute(std::vector<Lepton>& leptons, std::vector<Neutrino>& nu, std::vector<Jet>& jets, std::vector<Ljet>& ljets){
    /* Build top quarks system 
       - lepton
       - AK4 near lepton (2D cut)
         > highest pT
         > Ref: https://github.com/UHH2/UHH2/blob/master/common/src/Utils.cxx#L34
       - AK8 away from lepton
         > Most 'top-like' = highest BEST_t score
    */
    m_ttbar1L = {};

    // Setup lepton (only 1 in the single lepton analysis)
    cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(leptons.size())+" leptons");
    Lepton lep;
    if (leptons.size()>0)
        lep = leptons.at(0);
    else{
        lep.p4.SetPtEtaPhiE(0,0,0,0);   // this event will fail the selection anway, use dummy value
        lep.isGood = false;
    }

    // Get Jets (only if lepton exists)
    int ak4candidate(-1);   // index in jets that corresponds to AK4 from leptonic top
    int ak8candidate(-1);   // index in ljets that corresponds to AK8 (hadronic top)

    if (leptons.size()>0){
        // -- Setup AK4 jet : 2D Cut
        cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(jets.size())+" ak4 candidates");
        float ak4_pt(0);

        for (auto& jet : jets){
            float jpt   = jet.p4.Pt();
            float dr    = jet.p4.DeltaR(lep.p4);              // DeltaR( lepton,AK4 )
            float ptrel = cma::ptrel( lep.p4,jet.p4 );        // pTrel(  lepton,AK4 )

            // Two possible scenarios -- only consider jets with DeltaR<PI/2:
            // > jet within (0.4<DeltaR<PI/2)
            // > jet closer than 0.4 but ptrel>25
            // Choose highest pT option

            if (dr<M_HALF_PI) { 
                if (dr<0.4 && ptrel>25 && jpt>ak4_pt){
                    ak4_pt = jpt;
                    ak4candidate = jet.index;
                }
                else if (dr>0.4 && jpt>ak4_pt){
                    ak4_pt = jpt;
                    ak4candidate = jet.index;
                }
            } // end possible AK4 candidate
        } // end loop over ak4 


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


    // Define objects in the struct
    m_ttbar1L.neutrino = nu.at(0);
    m_ttbar1L.lepton   = lep;

    Jet dummy_jet;
    dummy_jet.isGood = false;
    m_ttbar1L.jet = (ak4candidate>=0) ? jets.at(ak4candidate) : dummy_jet;

    Ljet dummy_ljet;
    dummy_ljet.isGood = false;
    m_ttbar1L.ljet = (ak8candidate>=0) ? ljets.at(ak8candidate) : dummy_ljet;

    cma::DEBUG("TTBARRECO : Ttbar built ");

    return;
}


// dilepton
void ttbarReco::execute(std::vector<Electron>& electrons, std::vector<Muon>& muons, std::vector<Jet>& jets){
    /* Dilepton ttbar reconstruction */
    m_dileptonTtbar->execute(m_dilepton);
    return;
}

// THE END //
