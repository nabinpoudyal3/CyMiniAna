/*
Created:        --
Last Updated:   25 April 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Tool for reconstructing the neutrino
- 1-lepton: Use W-mass constraint
*/
#include "Analysis/CyMiniAna/interface/neutrinoReco.h"


neutrinoReco::neutrinoReco( configuration& cmaConfig ) :
  m_config(&cmaConfig){
    m_nu = {};
    m_lepton = {};
    m_met = {};
  }

neutrinoReco::~neutrinoReco() {}


void neutrinoReco::setObjects(Lepton& lepton, MET& met){
    /* Set lepton and MET from one function call */
    setLepton(lepton);
    setMET(met);

    return;
}

void neutrinoReco::setLepton(Lepton& lepton){
    /* Set the lepton (easily test different methods without passing this repeatedly) */
    m_lepton = lepton;
    return;
}

void neutrinoReco::setMET(MET& met){
    /* Set the MET (easily test different methods without passing this repeatedly) */
    m_met = met;
    return;
}


Neutrino neutrinoReco::execute(float wmass){
    /* Build the neutrino
       - Use the W-mass constraint
         See AN2015-107-v9 (Equation 3)
         > http://cms.cern.ch/iCMS/user/noteinfo?cmsnoteid=CMS%20AN-2015/107
         For imaginary solutions, take the real part
         For multiple real solutions, choose the smallest one
    */
    m_pz_solutions.clear();          // keep track of pz solutions (in case you want them all later)
    m_nu.p4.SetPtEtaPhiM(m_met.p4.Pt(),0.,m_met.p4.Phi(),0.);

    cma::DEBUG("NEUTRINORECO : Reconstructing the neutrino");

    float lepPt = m_lepton.p4.Pt();
    float nuPt  = m_nu.p4.Pt();
    float metPx = m_met.p4.Px();
    float metPy = m_met.p4.Py();

    float mu = 0.5 * pow(wmass,2) + (lepPt * nuPt);

    float A = -1. * pow(lepPt,2);
    float B = mu * m_lepton.p4.Pz();
    float C = pow(mu,2) - pow(m_lepton.p4.E(),2) * pow(nuPt,2);

    float discriminant = pow(B,2) - A*C;      // radicand

    if (discriminant<0) {
        // Imaginary! Take the real part of the solution for pz
        float pz  = -B/A;
        float nuE = sqrt( pow(metPx,2) + pow(metPy,2) + pow(pz,2));
        m_nu.p4.SetPxPyPzE( metPx, metPy, pz, nuE );
        m_pz_solutions.push_back(pz);
    }
    else {
        discriminant = sqrt(discriminant);
        float pz1 = (-B-discriminant) / A;
        float pz2 = (-B+discriminant) / A;

        float pz  = (pz1<pz2) ? pz1 : pz2;   // choose the smallest solution
        float nuE = sqrt( pow(metPx,2) + pow(metPy,2) + pow(pz,2));
        m_nu.p4.SetPxPyPzE( metPx, metPy, pz, nuE );

        m_pz_solutions.push_back(pz1);
        m_pz_solutions.push_back(pz2);
    }

    return m_nu;
}


std::vector<float> neutrinoReco::pzSolutions(){
    /* Return all the solutions to the user -- then they can use whichever they prefer */
    return m_pz_solutions;
}


// THE END //
