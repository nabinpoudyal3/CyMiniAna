/*
Created:        25 January   2018
Last Updated:   25 January   2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109

-----

Match reconstructed objects to truth objects

*/
#include "Analysis/CyMiniAna/interface/truthMatching.h"


truthMatching::truthMatching(configuration &cmaConfig) : 
  m_config(&cmaConfig){
  }

truthMatching::~truthMatching() {}

void truthMatching::initialize(){
    m_truth_tops.clear();
    m_truth_partons.clear();
    return;
}


void truthMatching::setTruthPartons(const std::vector<Parton> truth_partons){
    /* Set truth partons */
    m_truth_partons = truth_partons;
    return;
}


void truthMatching::setTruthTops(const std::vector<TruthTop> truth_tops){
    /* Set truth tops */
    m_truth_tops = truth_tops;
    return;
}


void truthMatching::matchJetToTruthTop(Jet& jet){
    /* Match reconstructed or truth jets to partons
       Currently setup to process partons from top quarks (qqb)
    */
    jet.matchId     = -1;
    jet.containment = 0;         // initialize containment
    jet.truth_partons.clear();

    for (unsigned int t_idx=0, size=m_truth_tops.size(); t_idx<size; t_idx++){
        auto truthtop = m_truth_tops.at(t_idx);
        if (!truthtop.isHadronic) continue;         // only want hadronically-decaying tops

        Parton bottomQ = m_truth_partons.at( truthtop.bottom );
        Parton wdecay1 = m_truth_partons.at( truthtop.Wdecays.at(0) );
        Parton wdecay2 = m_truth_partons.at( truthtop.Wdecays.at(1) );

        parton_match(bottomQ,jet);
        parton_match(wdecay1,jet);
        parton_match(wdecay2,jet);

        // if the jet is matched to a truth top, exit
        if (jet.containment!=0){
            cma::DEBUG("TRUTHMATCHING : Jet deltaR bottomQ = "+std::to_string(jet.p4.DeltaR(bottomQ.p4)));
            cma::DEBUG("TRUTHMATCHING : Jet deltaR wdecay1 = "+std::to_string(jet.p4.DeltaR(wdecay1.p4)));
            cma::DEBUG("TRUTHMATCHING : Jet deltaR wdecay2 = "+std::to_string(jet.p4.DeltaR(wdecay2.p4)));

            jet.matchId = t_idx;
            break;
        }
    }

    return;
}


void truthMatching::parton_match(const Parton& p, Jet& r, double dR){
    /* Match parton to physics object 
       @param p    truth parton
       @param r    reconstructed object
       @param dR   distance measure to use in DeltaR (default = -1)
    */
    bool match = (dR>0) ? cma::deltaRMatch( p.p4, r.p4, dR ) : cma::deltaRMatch( p.p4, r.p4, r.radius );

    // if matched update parameters of the object
    // check matches -- use map in header to avoid errors misremembering the integer values
    if (match){
        cma::DEBUG("TRUTHMATCH : parton_match() "+std::to_string(p.index)+" with containment "+std::to_string(p.containment));
        r.truth_partons.push_back(p.index);
        r.containment += p.containment;
    }

    return;
}


void truthMatching::matchJetToTruthJet(Jet& jet, const std::vector<Jet>& truth_jets){
    /* Match reco jet to truth jet */
    float truthDR( jet.radius );

    for (const auto& truth_jet : truth_jets){
        float thisDR = truth_jet.p4.DeltaR( jet.p4 );
        if ( thisDR < truthDR ){
            jet.truth_jet = truth_jet.index;
            truthDR = thisDR;
        }
    } // end loop over truth jets

    return;
}

// THE END

