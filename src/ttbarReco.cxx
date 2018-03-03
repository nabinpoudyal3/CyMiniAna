/*
Created:        --
Last Updated:   19 February 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Tool for performing deep learning tasks
*/
#include "Analysis/CyMiniAna/interface/ttbarReco.h"


ttbarReco::ttbarReco( configuration& cmaConfig ) :
  m_config(&cmaConfig){
    m_mapContainment = m_config->mapOfPartonContainment();  // containment map (ints and strings)
    m_targetMap = m_config->mapOfTargetValues();
  }

ttbarReco::~ttbarReco() {}


std::vector<Top> ttbarReco::tops(){
    /* Return the ttbar system */
    return m_ttbar;
}


void ttbarReco::execute(const std::vector<Jet>& jets, const std::vector<Ljet>& ljets){
    /* Build top quarks system
       Testing semi-resolved tagger; only interested in AK8(QB/W)+AK4
        0 :: NONE = QCD (background)
        1 :: QB-Q = Signal AK8(QB) + AK4(Q)
        2 :: QQ-B = Signal AK8(W)  + AK4(B)
        3 :: FULL = Signal AK8
        4 :: Other = placeholder

        -- if (fully contained): top_cand.target = 'full'; no extra AK4 needed
        -- if (q/b-only) top_cand.target = 4; need two extra AK4 -> this is effectively the 'fully resolved' case
    */
    m_ttbar.clear();

    m_jets  = jets;
    m_ljets = ljets;


    // Reconstruct ttbar, define containment
    int FULL  = m_mapContainment.at("FULL");
    int QB    = m_mapContainment.at("BQ");
    int W     = m_mapContainment.at("W");
    int QONLY = m_mapContainment.at("QONLY");
    int BONLY = m_mapContainment.at("BONLY");

    bool isQCD(m_config->isQCD());
    bool isTtbar(m_config->isTtbar());

    unsigned int target(4);
    if (isQCD) target = 0;

    HadTop top_cand;  // reconstructed top candidates

    cma::DEBUG("TTBARRECO : building ttbar with "+std::to_string(m_ljets.size())+" ak8 candidates");
    for (const auto& ljet : ljets){

        // Overlap Removal (subtract AK4 if matched to AK8 subjet)
        std::vector<Jet> ak4candidates;
        overlapRemoval(ljet,ak4candidates,isTtbar);  // overlap removal

        int absLjetCt = std::abs(ljet.containment);

        cma::DEBUG("TTBARRECO : Access AK8 jet matched to parton "+std::to_string(ljet.matchId)+"; containment = "+std::to_string(absLjetCt));

        // reset struct
        top_cand.jets.clear();
        top_cand.ljet = ljet.index;

        if (isTtbar && (absLjetCt==QB || absLjetCt==W)){
            // find the AK4 jet that creates "FULL" with AK8
            cma::DEBUG("TTBARRECO : -- AK8 that is QB or W");
            target = (absLjetCt==QB) ? m_targetMap["QB"] : m_targetMap["W"];
            for (const auto& jet : ak4candidates){

                // ensure jets are matched to the same parton
                if (jet.matchId!=ljet.matchId) continue;
                // Don't use an AK4 jet that matches to the same parton as the AK8
                if (std::find(ljet.truth_partons.begin(), ljet.truth_partons.end(), jet.truth_partons.at(0)) != ljet.truth_partons.end())
                    continue;

                int total_containment = jet.containment+ljet.containment;

                if (total_containment==FULL || total_containment==-FULL){
                    cma::DEBUG("TTBARRECO : AK8+AK4 top in ttbar");
                    top_cand.jets.push_back(jet.index);
                    top_cand.target = target;
                    m_ttbar.push_back( top_cand );
                }
            } // end loop over AK4
        } // end if ljet is QB or W
        else if (isQCD){
            // Just assign AK8+AK4 jets as background (target 0)
            for (const auto& jet : ak4candidates){
                top_cand.jets.clear();
                top_cand.ljet   = ljet.index;
                top_cand.target = target;
                top_cand.jets.push_back(jet.index);

                m_ttbar.push_back( top_cand );
            }
        } // end if QCD
    } // end loop over ak8 candidates

    cma::DEBUG("TTBARRECO : Ttbar built ");

    return;
}




void ttbarReco::overlapRemoval(const Ljet& ak8, std::vector<Jet>& new_objects, const bool match_truth){
    /* 
      Remove AK4 jets that overlap with AK8 candidate 
      'match_truth' for checking the AK4 and AK8 match the same parton
    */
    new_objects.clear();

    for (const auto& jet : m_jets){

        // check if AK4 is DeltaR matched to AK8 candidate subjets
        bool subjet_matched(false);

        bool subjet_match(false);
        if (ak8.subjets.size()<1)  // match to AK8 if there are no subjets
            subjet_match = cma::deltaRMatch(ak8.p4, jet.p4, ak8.radius);
        else{
            for (const auto& sj : ak8.subjets)
                subjet_match = cma::deltaRMatch(sj.p4, jet.p4, sj.radius);
        }

        if (match_truth && ak8.matchId == jet.matchId && subjet_match){
           subjet_matched = true;
           break;
        }
        else if (!match_truth && subjet_match){
           subjet_matched = true;
           break;
        }
        // end if DeltaR matched

        if (!subjet_matched) new_objects.push_back( jet );
    }

    return;
}


// THE END //

