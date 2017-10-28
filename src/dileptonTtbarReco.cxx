/*
Created:        --
Last Updated:   22 October  2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

 Methods for dileptonTtbarReco
  Class for dilepton ttbar reconstruction.
  Imported from
   https://gitlab.cern.ch/cms-desy-top/TopAnalysis/blob/master/
           Configuration/analysis/common/src/KinematicReconstruction.cc

 Main interface to the dilepton reconstruction

*/
#include "diHiggs/CyMiniAna/interface/dileptonTtbarReco.h"



dileptonTtbarReco::dileptonTtbarReco(configuration& cmaConfig,
                                     const configuration::Era era, const int minNumberOfBtags, 
                                     const bool preferBtags,       const bool massLoop):
  m_config(&cmaConfig),
  m_era(era),
  m_minNumberOfBtags(minNumberOfBtags),
  m_preferBtags(preferBtags),
  m_massLoop(massLoop),
  m_NSol(0),
  m_h_wmass(0),
  m_h_jetAngleRes(0),
  m_h_jetEres(0),
  m_h_lepAngleRes(0),
  m_h_lepEres(0),
  m_h_mbl_w(0){
    if(m_minNumberOfBtags<0 || m_minNumberOfBtags>2){
        cma::ERROR("DILEPTON : Minimum number of b-tags needs to be within [0,2]. Value is: "+std::to_string(m_minNumberOfBtags));
        std::exit(96);
    }

    m_r3 = new TRandom3();
    m_btag_wp   = m_config->jet_btagWkpt();  // 0.244; FIXME:Warning , hardcoded value of b-tag working point. 
    m_rangeLow  = m_config->massMin();       // 173;
    m_rangeHigh = m_config->massMax();       // 300;
    m_topMass   = m_config->topQuarkMass();  // 172.5

    // Read all histograms for smearing objects
    if (!m_massLoop) loadData();
  }

dileptonTtbarReco::~dileptonTtbarReco() {
    delete m_r3;
}


std::map<std::string,Top> dileptonTtbarReco::execute(const DileptonReco& ttSystem){
    /* Pass information from event to this class and return the solution */
    TtbarDilepton ttbarSystem;
    ttbarSystem.lepton_pos = ttSystem.lepton_pos;
    ttbarSystem.lepton_neg = ttSystem.lepton_neg;
    ttbarSystem.met        = ttSystem.met;
    ttbarSystem.jets       = ttSystem.jets;
    ttbarSystem.bjets      = ttSystem.bjets;

    dileptonTtbarRecoSolution ttDileptonSolution = buildTtbar(ttbarSystem);
    // convert this solution to a map for use in CyMiniAna

    std::map<std::string,Top> ttbarDileptonSolution;  // solution for the Event class
    if (ttDileptonSolution.numberOfSolutions()<1)
        return ttbarDileptonSolution;

    TtbarDilepton soln = ttDileptonSolution.solution(TtbarDilepton::defaultForMethod,0); // take the first solution
    cma::DEBUG("DILEPTON : Solution obtained, now save");

    ttbarDileptonSolution["top"]     = soln.top;
    ttbarDileptonSolution["antitop"] = soln.topBar;

    cma::DEBUG("DILEPTON : Solution obtained now exit");

    return ttbarDileptonSolution;
}


dileptonTtbarRecoSolution dileptonTtbarReco::buildTtbar(const TtbarDilepton& ttSystem){
    /* Internal function to solve for the ttbar dilepton system */
    dileptonTtbarRecoSolution result;

    // Check if minimum number of objects for any solution is present
    const int numberOfJets(ttSystem.jets.size());
    const int numberOfBtags(ttSystem.bjets.size());
    if(numberOfJets<2 || numberOfBtags<m_minNumberOfBtags) {
        cma::ERROR("DILEPTON : Not enough jets or not enough b-tags! ");
        return result;
    }

    Lepton lepton     = ttSystem.lepton_neg;
    Lepton antiLepton = ttSystem.lepton_pos;

    // -- Find solutions with 2 b-tagged jets
    if (numberOfBtags>=2){
        for (int i_index=0; i_index<numberOfBtags-1; ++i_index) {
            for (int j_index=i_index+1; j_index<numberOfBtags; ++j_index) {
                cma::DEBUG("DILEPTON : Check indices "+std::to_string(i_index)+", "+std::to_string(j_index));

                Jet jet1 = ttSystem.bjets.at(i_index);
                Jet jet2 = ttSystem.bjets.at(j_index);

                TtbarDilepton tt_tmp1;
                tt_tmp1.lepton_pos = antiLepton;
                tt_tmp1.lepton_neg = lepton;
                tt_tmp1.bJet    = jet1;
                tt_tmp1.bbarJet = jet2;
                tt_tmp1.met     = ttSystem.met;

                TtbarDilepton tt_tmp2;           // same as tt_tmp1 with the jets switched
                tt_tmp2.lepton_pos = antiLepton;
                tt_tmp2.lepton_neg = lepton;
                tt_tmp2.bJet    = jet2;
                tt_tmp2.bbarJet = jet1;
                tt_tmp2.met     = ttSystem.met;

                const std::vector<TtbarDilepton> solution1( getSolutions(tt_tmp1,2) );
                const std::vector<TtbarDilepton> solution2( getSolutions(tt_tmp2,2) );
                cma::DEBUG("DILEPTON : Add solutions");
                result.Add(solution1);
                result.Add(solution2);
                cma::DEBUG("DILEPTON : Added solutions");
            }
        }
    }

    cma::DEBUG("DILEPTON : end 2 b-tags");
    if(m_preferBtags && result.numberOfSolutions(2)){
        cma::DEBUG("DILEPTON : Prefer b-tags and number of solutions");
        return result;
    }
    if(m_minNumberOfBtags > 1) {
        // require at least 2 b-tags, but not enough solutions; exit.
        cma::DEBUG("DILEPTON : Prefer b-tags > 1");
        return result;
    }

    cma::DEBUG("DILEPTON : Continuing to 1 b-tag solutions");


    // Collect non b-tagged jets for reconstruction
    const int numberOfNonBtags(numberOfJets - numberOfBtags);
    if(numberOfNonBtags < 1) 
        return result;

    if (numberOfBtags==1){
        Jet bJet = ttSystem.bjets.at(0);

        // -- Find solutions with 1 b-tagged jet
        for(const auto& jet : ttSystem.jets) {
            if (jet.isbtagged.at(m_btag_wp)) continue;  // skip the b-tagged jet

            TtbarDilepton tt_tmp1;
            tt_tmp1.lepton_pos = antiLepton;
            tt_tmp1.lepton_neg = lepton;
            tt_tmp1.bJet    = bJet;
            tt_tmp1.bbarJet = jet;
            tt_tmp1.met     = ttSystem.met;

            TtbarDilepton tt_tmp2;           // same as tt_tmp1 with the jets switched
            tt_tmp2.lepton_pos = antiLepton;
            tt_tmp2.lepton_neg = lepton;
            tt_tmp2.bJet    = jet;
            tt_tmp2.bbarJet = bJet;
            tt_tmp2.met     = ttSystem.met;

            const std::vector<TtbarDilepton> solution1( getSolutions(tt_tmp1,1) );
            const std::vector<TtbarDilepton> solution2( getSolutions(tt_tmp2,1) );
            result.Add(solution1);
            result.Add(solution2);
        }
    }

    if(m_preferBtags && result.numberOfSolutions(1)) 
        return result;

    if(m_minNumberOfBtags > 0) 
        // require >0 b-tags, but the 1,2 b-tag solutions don't exist!
        return result;

    // -- Find solutions with 0 b-tagged jets
    if (numberOfBtags==0){
        for (int i_index=0; i_index<numberOfJets-1; ++i_index) {
            for (int j_index=i_index+1; i_index<numberOfJets; ++j_index) {
                Jet jet1 = ttSystem.jets.at(i_index);
                Jet jet2 = ttSystem.jets.at(j_index);

                TtbarDilepton tt_tmp1;
                tt_tmp1.lepton_pos = antiLepton;
                tt_tmp1.lepton_neg = lepton;
                tt_tmp1.bJet    = jet1;
                tt_tmp1.bbarJet = jet2;
                tt_tmp1.met     = ttSystem.met;

                TtbarDilepton tt_tmp2;           // same as tt_tmp1 with the jets switched
                tt_tmp2.lepton_pos = antiLepton;
                tt_tmp2.lepton_neg = lepton;
                tt_tmp2.bJet    = jet2;
                tt_tmp2.bbarJet = jet1;
                tt_tmp2.met     = ttSystem.met;

                const std::vector<TtbarDilepton> solution1( getSolutions(tt_tmp1,0) );
                const std::vector<TtbarDilepton> solution2( getSolutions(tt_tmp2,0) );
                result.Add(solution1);
                result.Add(solution2);
            }
        }
    }

    return result;
}



std::vector<TtbarDilepton> dileptonTtbarReco::getSolutions(const TtbarDilepton& singleTtbarSystem, const int numberOfBtags) {
    /* Obtain the solutions from Utils*/
    std::vector<TtbarDilepton> result;

    // assign objects to the (simple) struct of information
    DileptonReco ttSystem;
    ttSystem.lepton_pos = singleTtbarSystem.lepton_pos;
    ttSystem.lepton_neg = singleTtbarSystem.lepton_neg;
    ttSystem.bJet       = singleTtbarSystem.bJet;
    ttSystem.bbarJet    = singleTtbarSystem.bbarJet;
    ttSystem.met        = singleTtbarSystem.met;

    if(m_massLoop){
        // resolve the ttbar kinematics assuming different values for top mass
        double neutrinoWeightMax(0.);

        dileptonTtbarRecoUtils tp_m;
        for(double iTopMass=m_rangeLow; iTopMass<m_rangeHigh; iTopMass+=1){

            tp_m.setTopMass( iTopMass );
            tp_m.setAntiTopMass( iTopMass );
            tp_m.execute(ttSystem);

            if(tp_m.getNsol() < 1)
                continue;

            // Save the solution! Only take the solution with largest weight (first entry)
            TtbarDilepton tp_m_soln = tp_m.getTtSol().at(0);

            const double& weight = tp_m_soln.weight;
            if(weight <= neutrinoWeightMax)
                continue;

            neutrinoWeightMax = weight;

            std::map<TtbarDilepton::WeightType, double> tmp_weight;
            tmp_weight[TtbarDilepton::defaultForMethod] = weight;
            tmp_weight[TtbarDilepton::neutrinoEnergy]   = weight;

            tp_m_soln.top.weight_tt    = weight;
            tp_m_soln.topBar.weight_tt = weight;
            tp_m_soln.top.weight       = -999.;
            tp_m_soln.topBar.weight    = -999.;
            tp_m_soln.top.weight_ES    = -999.;
            tp_m_soln.topBar.weight_ES = -999.;

            tp_m_soln.mapOfWeights = tmp_weight;
            tp_m_soln.isNoSmearSol = false;
            tp_m_soln.recMtop = iTopMass;
            tp_m_soln.ntags   = numberOfBtags;

            result.push_back(tp_m_soln);
        }
        tp_m.fDelete();
    }
    else{
        // solve for the ttbar kinematics after smearing inputs
        std::vector<TtbarDilepton> smeared_solutions;
        std::vector<double> smeared_weights;
        bool smearSolution = getSmearedSolutions(smeared_solutions,smeared_weights,ttSystem);

        // Check if there is also a solution without smearing
        dileptonTtbarRecoUtils tp_NOsm;
        tp_NOsm.execute(ttSystem);
        bool noSmearSolution = (tp_NOsm.getNsol()>0);

        std::map<TtbarDilepton::WeightType, double> tmp_weight;

        if (smearSolution && !noSmearSolution){
            // Solution with smearing
            double weight = std::accumulate(smeared_weights.begin(),smeared_weights.end(),0.0);

            tmp_weight[TtbarDilepton::defaultForMethod]         = weight;
            tmp_weight[TtbarDilepton::averagedSumSmearings_mlb] = weight;

            // from the collection of smeared solutions, average the four-vectors
            // for the top/anti-top/neutrino/anti-neutrino
            TtbarDilepton tp_m_soln;
            averageSmearedSolutions(smeared_solutions,smeared_weights,tp_m_soln);

            // Add the extra attributes to this solution
            tp_m_soln.top.weight_tt    = weight;
            tp_m_soln.topBar.weight_tt = weight;

            tp_m_soln.mapOfWeights = tmp_weight;
            tp_m_soln.isNoSmearSol = false;
            tp_m_soln.recMtop = m_topMass;
            tp_m_soln.ntags   = numberOfBtags;

            // keep the original objects (even though smeared inputs achieved the solution!)
            tp_m_soln.lepton_pos = singleTtbarSystem.lepton_pos;
            tp_m_soln.lepton_neg = singleTtbarSystem.lepton_neg;
            tp_m_soln.bJet       = singleTtbarSystem.bJet;
            tp_m_soln.bbarJet    = singleTtbarSystem.bbarJet;
            tp_m_soln.met        = singleTtbarSystem.met;

            // save the solution
            result.push_back(tp_m_soln);
        }
        else if (noSmearSolution) {
            // Solution without smearing (even if smeared solution exists!)
            TtbarDilepton tp_m_soln = tp_NOsm.getTtSol().at(0);

            const double weight = tp_m_soln.weight;
            tmp_weight[TtbarDilepton::defaultForMethod]         = weight;
            tmp_weight[TtbarDilepton::averagedSumSmearings_mlb] = weight;

            tp_m_soln.top.weight_tt    = weight;
            tp_m_soln.topBar.weight_tt = weight;

            tp_m_soln.mapOfWeights = tmp_weight;
            tp_m_soln.isNoSmearSol = true;
            tp_m_soln.recMtop = m_topMass;
            tp_m_soln.ntags   = numberOfBtags;

            result.push_back(tp_m_soln);
        }

        tp_NOsm.fDelete(); // clean-up memory
    }

    cma::DEBUG("DILEPTON : End getSolutions");

    return result;
}


bool dileptonTtbarReco::getSmearedSolutions(std::vector<TtbarDilepton>& smeared_solutions,
                                            std::vector<double>& smeared_weights,
                                            const DileptonReco& ttSystem) const {
    /* Determine solution to ttbar reconstruction after smearing the inputs 
    
    ++++ Don't know what this is for, removing it ++++
    if((antiLepton.p4+jet1.p4).M()>180. || (lepton.p4+jet2.p4).M()>180.)
        // protection against large masses?
        return false;
    */
    const unsigned int seed = cma::setRandomNumberSeeds(ttSystem.lepton_neg,ttSystem.lepton_neg,ttSystem.bJet,ttSystem.bbarJet);
    m_r3->SetSeed(seed);  // Set random number generator seeds based on this arrangement

    smeared_solutions.clear();

    bool isSolved(false);         // true if at least 1 smeared iteration produces a solution
    double err(0.001);            // tolerance for angle rotation
    dileptonTtbarRecoUtils tp_sm; // object responsible for finding the ttbar solution

    Lepton lepton     = ttSystem.lepton_neg;
    Lepton antiLepton = ttSystem.lepton_pos;
    Jet jet1     = ttSystem.bJet;
    Jet jet2     = ttSystem.bbarJet;
    TVector2 met = ttSystem.met;
    TLorentzVector met_tv;
    met_tv.SetPxPyPzE( met.Px(), met.Py(), 0., 0.);

    TVector3 vX_reco = -jet1.p4.Vect() - jet2.p4.Vect() - lepton.p4.Vect() - antiLepton.p4.Vect() - met_tv.Vect();

    for(int sm=0,size=m_config->NJetSmear(); sm<size; ++sm){
        Jet b_sm;
        Jet bbar_sm;
        Lepton l_sm;
        Lepton al_sm;
        TVector2 met_sm;

        // jets energy smearing values
        double fB    = m_h_jetEres->GetRandom();
        double xB    = sqrt((fB*fB*pow(jet1.p4.E(),2)-jet1.p4.M2())/(pow(jet1.p4.P(),2)));
        double fBbar = m_h_jetEres->GetRandom();
        double xBbar = sqrt((fBbar*fBbar*pow(jet2.p4.E(),2)-jet2.p4.M2())/(pow(jet2.p4.P(),2)));

        // b-jet smearing
        b_sm.p4.SetXYZT(jet1.p4.Px()*xB,jet1.p4.Py()*xB,jet1.p4.Pz()*xB,jet1.p4.E()*fB);
        tp_sm.angle_rot(m_h_jetAngleRes->GetRandom(),err,b_sm,b_sm);

        // bbar-jet smearing
        bbar_sm.p4.SetXYZT(jet2.p4.Px()*xBbar,jet2.p4.Py()*xBbar,jet2.p4.Pz()*xBbar,jet2.p4.E()*fBbar);
        tp_sm.angle_rot(m_h_jetAngleRes->GetRandom(),err,bbar_sm,bbar_sm);


        // leptons energy smearing values
        double fL  = m_h_lepEres->GetRandom();
        double xL  = sqrt((fL*fL*pow(lepton.p4.E(),2)-lepton.p4.M2())/(pow(lepton.p4.P(),2)));
        double faL = m_h_lepEres->GetRandom();
        double xaL = sqrt((faL*faL*pow(antiLepton.p4.E(),2)-antiLepton.p4.M2())/( pow(antiLepton.p4.P(),2)));

        // lepton smearing
        l_sm.p4.SetXYZT(lepton.p4.Px()*xL,lepton.p4.Py()*xL,lepton.p4.Pz()*xL,lepton.p4.E()*fL);
        tp_sm.angle_rot(m_h_lepAngleRes->GetRandom(),err,l_sm,l_sm);

        // anti-lepton smearing
        al_sm.p4.SetXYZT(antiLepton.p4.Px()*xaL,antiLepton.p4.Py()*xaL,antiLepton.p4.Pz()*xaL,antiLepton.p4.E()*faL);
        tp_sm.angle_rot(m_h_lepAngleRes->GetRandom(),err,al_sm,al_sm);

        // Get smeared MET
        TVector3 metV3_sm = -b_sm.p4.Vect()-bbar_sm.p4.Vect()-l_sm.p4.Vect()-al_sm.p4.Vect()-vX_reco;
        met_sm.SetX(metV3_sm.Px());
        met_sm.SetY(metV3_sm.Py());


        // Get new solution! (random masses for the W)
        double wplus_sm  = m_h_wmass->GetRandom();
        double wminus_sm = m_h_wmass->GetRandom();
        tp_sm.setWplus( wplus_sm );
        tp_sm.setWminus( wminus_sm );

        DileptonReco sm_ttSystem;
        sm_ttSystem.lepton_pos = al_sm;
        sm_ttSystem.lepton_neg = l_sm;
        sm_ttSystem.bJet       = b_sm;
        sm_ttSystem.bbarJet    = bbar_sm;
        sm_ttSystem.met        = met_sm;

        tp_sm.execute(sm_ttSystem);

        if(tp_sm.getNsol()>0) {
            isSolved = true;    // at least 1 smeared solution gives a result
            double mbl_weight = m_h_mbl_w->GetBinContent(m_h_mbl_w->FindBin((al_sm.p4+b_sm.p4).M())) * 
                                                         m_h_mbl_w->GetBinContent(m_h_mbl_w->FindBin((l_sm.p4+bbar_sm.p4).M())) / 
                                                         100000000; // where does this come from?
            smeared_solutions.push_back( tp_sm.getTtSol().at(0) ); // only take the first solution (highest weight)
            smeared_weights.push_back( mbl_weight );
        }
    }

//    tp_sm.fDelete(); // clean-up memory

    return isSolved;
}


void dileptonTtbarReco::averageSmearedSolutions(const std::vector<TtbarDilepton>& smeared_solutions,
                                                const std::vector<double>& smeared_weights,
                                                TtbarDilepton& solution) {
    /* Average the four-vectors for tops & neutrinos */
    double px_sum_top(0);    double px_sum_antitop(0);
    double py_sum_top(0);    double py_sum_antitop(0);
    double pz_sum_top(0);    double pz_sum_antitop(0);
    double px_sum_nu(0);     double px_sum_antinu(0);
    double py_sum_nu(0);     double py_sum_antinu(0);
    double pz_sum_nu(0);     double pz_sum_antinu(0);

    double sumOfWeights(0.0);

    for(int i=0,size=smeared_weights.size();i<size;++i) {

          double w(smeared_weights.at(i));
          TtbarDilepton thisSolution = smeared_solutions.at(i);

          px_sum_top     += w * thisSolution.top.p4.Px();
          py_sum_top     += w * thisSolution.top.p4.Py();
          pz_sum_top     += w * thisSolution.top.p4.Pz();

          px_sum_antitop += w * thisSolution.topBar.p4.Px();
          py_sum_antitop += w * thisSolution.topBar.p4.Py();
          pz_sum_antitop += w * thisSolution.topBar.p4.Pz();

          px_sum_nu      += w * thisSolution.neutrino.p4.Px();
          py_sum_nu      += w * thisSolution.neutrino.p4.Py();
          pz_sum_nu      += w * thisSolution.neutrino.p4.Pz();

          px_sum_antinu  += w * thisSolution.neutrinoBar.p4.Px();
          py_sum_antinu  += w * thisSolution.neutrinoBar.p4.Py();
          pz_sum_antinu  += w * thisSolution.neutrinoBar.p4.Pz();

          sumOfWeights += w;
    }

    double px_top     = px_sum_top/sumOfWeights;
    double py_top     = py_sum_top/sumOfWeights;
    double pz_top     = pz_sum_top/sumOfWeights;
    double px_antitop = px_sum_antitop/sumOfWeights;
    double py_antitop = py_sum_antitop/sumOfWeights;
    double pz_antitop = pz_sum_antitop/sumOfWeights;
    double px_nu      = px_sum_nu/sumOfWeights;
    double py_nu      = py_sum_nu/sumOfWeights;
    double pz_nu      = pz_sum_nu/sumOfWeights;
    double px_antinu  = px_sum_antinu/sumOfWeights;
    double py_antinu  = py_sum_antinu/sumOfWeights;
    double pz_antinu  = pz_sum_antinu/sumOfWeights;

    Neutrino nu;
    nu.p4.SetXYZM(px_nu,py_nu,pz_nu,0.0);

    Neutrino antiNu;
    antiNu.p4.SetXYZM(px_antinu,py_antinu,pz_antinu,0.0);

    // redeclare the other attributes of the struct (to keep information available)
    Lepton lep  = smeared_solutions.at(0).lepton_neg; // same for all solutions
    Lepton alep = smeared_solutions.at(0).lepton_pos; 

    TLorentzVector wplus  = alep.p4 + nu.p4;
    TLorentzVector wminus = lep.p4  + antiNu.p4;

    Top top;
    top.lepton   = alep;
    top.neutrino = nu;
    top.jet      = smeared_solutions.at(0).bJet;
    top.p4.SetXYZM(px_top,py_top,pz_top,m_topMass);

    Top antiTop;
    antiTop.lepton   = lep;
    antiTop.neutrino = antiNu;
    antiTop.jet      = smeared_solutions.at(0).bbarJet;
    antiTop.p4.SetXYZM(px_antitop,py_antitop,pz_antitop,m_topMass);

    TLorentzVector m_tt = top.p4 + antiTop.p4;

    solution.top         = top;
    solution.topBar      = antiTop;
    solution.Wplus       = wplus;
    solution.Wminus      = wminus;
    solution.neutrino    = nu;
    solution.neutrinoBar = antiNu;

    solution.x1     = (top.p4.E()+antiTop.p4.E()+top.p4.Pz()+antiTop.p4.Pz())/(m_config->beamEnergy());
    solution.x2     = (top.p4.E()+antiTop.p4.E()-top.p4.Pz()-antiTop.p4.Pz())/(m_config->beamEnergy());
    solution.mtt    = m_tt.M();
    solution.weight = 1./m_tt.M();

    return;
}


void dileptonTtbarReco::setSolutions() {
    setSolutions(m_sols);
}


void dileptonTtbarReco::setSolutions(std::vector<TtbarDilepton> sols) {
    m_NSol = (int)(sols.size());

    if(m_NSol>0){
        std::nth_element(begin(sols), begin(sols), end(sols),
                         [](const TtbarDilepton& a, const TtbarDilepton& b){
                            return b.ntags < a.ntags || (b.ntags == a.ntags && b.weight < a.weight);
                         });
        m_sol = sols[0];
    }

    return;
}


int dileptonTtbarReco::getNSol() const {
    return m_NSol;
}

TtbarDilepton dileptonTtbarReco::getSol() const {
    return m_sol;
}

std::vector<TtbarDilepton> dileptonTtbarReco::getSols() const {
    return m_sols;
}



void dileptonTtbarReco::loadData() {
    /* Load data for smearing 
         jet & lepton resolutions
         mbl & W mass
    */
    std::string data_path = m_config->getAbsolutePath()+"/config";

    if (m_era == configuration::run2_13tev_25ns || m_era == configuration::run2_13tev_2015_25ns)
        data_path += "/KinReco_input_2015_Run2CD_25ns_76X.root";
    else if (m_era == configuration::run2_13tev_2016_25ns)
        data_path += "/KinReco_input_2016_Run2BtoH_25ns_80X_v036_dilep_sel_25July2017.root";
    else{
        cma::ERROR("DILEPTON : loadData() : Era is not supported");
        std::exit(357);
    }

    TFile dataFile(data_path.c_str());

    //jet angle resolution
    m_h_jetAngleRes = (TH1F*)dataFile.Get("KinReco_d_angle_jet_step7");
    m_h_jetAngleRes->SetDirectory(0);

    //jet energy resolution
    m_h_jetEres = (TH1F*)dataFile.Get("KinReco_fE_jet_step7");
    m_h_jetEres->SetDirectory(0);

    //lep angle resolution
    m_h_lepAngleRes = (TH1F*)dataFile.Get("KinReco_d_angle_lep_step7");
    m_h_lepAngleRes->SetDirectory(0);

    //lep energy resolution
    m_h_lepEres = (TH1F*)dataFile.Get("KinReco_fE_lep_step7");
    m_h_lepEres->SetDirectory(0);

    //mbl mass
    m_h_mbl_w = (TH1F*)dataFile.Get("KinReco_mbl_true_step0");
    m_h_mbl_w->SetDirectory(0);

    // W mass
    m_h_wmass = (TH1F*)dataFile.Get("KinReco_W_mass_step0");
    m_h_wmass->SetDirectory(0);

    dataFile.Close();

    cma::DEBUG("DILEPTON : Found all histograms needed for smearing");

    return;
}




/**************************************************/
/**** NOT USED, BUT LEFT HERE FOR COMPLETENESS ****/
/**************************************************/

// void dileptonTtbarReco::inputNoJetMerging(std::vector<int>& b1_id, 
//                                           std::vector<int>& b2_id, 
//                                           std::vector<int>& nb_tag,
//                                           const std::vector<double>& btags) const {
//     std::vector<double> dummy;
//     inputNoJetMerging(b1_id, b2_id, nb_tag, btags, dummy);
// 
//     return;
// }
// 
// void dileptonTtbarReco::inputNoJetMerging(std::vector<int>& b1_id, 
//                                           std::vector<int>& b2_id, 
//                                           std::vector<int>& nb_tag,
//                                           const std::vector<double>& btags,
//                                           std::vector<double> btag_ww) const {
//     /* */
//     for (int i=0; i<(int)btags.size(); ++i){
//         for (int j=0; j<(int)btags.size(); ++j){
//             double wi = btags.at(i);
//             double wj = btags.at(j);
//             if(i==j || (wi<m_btag_wp && wj<m_btag_wp)) continue;
//             btag_ww.push_back(wi + wj);
// 
//             if(wi>m_btag_wp && wj>m_btag_wp) nb_tag.push_back(2);
//             else nb_tag.push_back(1);
// 
//             b1_id.push_back(i);
//             b2_id.push_back(j);
//         }
//     }
// 
//     return;
// }
// 
// 
// void dileptonTtbarReco::kinReco(const Lepton& leptonMinus,   const Lepton& leptonPlus, 
//                                 const std::vector<Jet> jets, const std::vector<double> btags, const TVector2 met) {
//     m_sols.clear();
// 
//     // Jets selection
//     std::vector<int> b1_id;
//     std::vector<int> b2_id;
//     std::vector<int> nb_tag;
//     std::vector<Jet> new_jets(jets);
//     std::vector<double> new_btags(btags);
//     this->inputNoJetMerging(b1_id, b2_id, nb_tag, btags);
// 
//     if(b1_id.size() < 2)
//         return; 
// 
//     dileptonTtbarRecoMeanSolution meanSolution(m_config->topQuarkMass());
// 
//     for(int ib=0; ib<(int)b1_id.size(); ++ib){
//         const int bjetIndex     = b1_id.at(ib);
//         const int antiBjetIndex = b2_id.at(ib);
//         const int numberOfBtags = nb_tag.at(ib);
// 
//         const Jet jet1 = new_jets.at(bjetIndex);
//         const Jet jet2 = new_jets.at(antiBjetIndex);
// 
//         const bool hasSolution(this->solutionSmearing(meanSolution,leptonMinus,leptonPlus,jet1,jet2,met));
// 
//         if(hasSolution){
//             meanSolution.getMeanSol(m_sol.top, m_sol.topBar, m_sol.neutrino, m_sol.neutrinoBar);
//             m_sol.weight = meanSolution.getSumWeight();
//             m_sol.Wplus  = m_sol.lp + m_sol.neutrino;
//             m_sol.Wminus = m_sol.lm + m_sol.neutrinoBar;
//             m_sol.ttbar  = m_sol.top + m_sol.topBar;
//             m_sol.jetB_index    = bjetIndex;
//             m_sol.jetBbar_index = antiBjetIndex;
//             m_sol.ntags = numberOfBtags;
//             m_sols.push_back(m_sol);
//         }
//         meanSolution.clear();
//     }
// 
//     this->setSolutions();
// 
//     return;
// }
// 
// void dileptonTtbarReco::kinRecoMassLoop(const Lepton& leptonMinus, const Lepton& leptonPlus, 
//                                         const std::vector<Jet> jets, const std::vector<double> btags, 
//                                         const TVector2 met) {
//     /* Reconstruction with optional mass loop */
//     std::vector<TtbarDilepton> vect_sol;
// 
//     const TLorentzVector leptonPlus_tlv  = leptonPlus.p4;
//     const TLorentzVector leptonMinus_tlv = leptonMinus.p4;
//     const TLorentzVector met_tlv;
//     met_tlv.SetPxPyPzE(met.Px(),met.Py(),0,0);
// 
//     std::vector<TLorentzVector> jets_tlv;
//     for (const auto& jet : jets) {
//         jets_tlv.push_back(jet.p4);
//     }
// 
//     std::vector<int> b1_id;
//     std::vector<int> b2_id;
//     std::vector<double> btag_ww;
//     std::vector<int> nb_tag;
// 
//     inputNoJetMerging(b1_id, b2_id, nb_tag, btags, btag_ww);
// 
//     if(b1_id.size()<2) {
//         m_NSol=0;
//         return;
//     }
// 
//     for(int i=0; i<(int)btag_ww.size()-1; ++i) {
//         if(btag_ww[i]>=btag_ww[i+1]) continue;
// 
//         double aux = btag_ww[i];
//         btag_ww[i] = btag_ww[i+1];
//         btag_ww[i+1] = aux;
// 
//         int aix = b1_id[i];
//         b1_id[i] = b1_id[i+1];
//         b1_id[i+1] = aix;
// 
//         aix = b2_id[i];
//         b2_id[i] = b2_id[i+1];
//         b2_id[i+1] = aix;
// 
//         aix = nb_tag[i];
//         nb_tag[i] = nb_tag[i+1];
//         nb_tag[i+1] = aix;
// 
//         i=-1;
//     }
// 
//     // jets loop
//     m_NSol=0;
//     for(int ib=0; ib<(int)b1_id.size(); ++ib) {
// 
//         int j1=b1_id[ib];
//         int j2=b2_id[ib];
// 
//         const TLorentzVector b_temp    = jets_tlv.at(j1);
//         const TLorentzVector bbar_temp = jets_tlv.at(j2);
//         const TLorentzVector l_temp    = leptonMinus_tlv;
//         const TLorentzVector al_temp   = leptonPlus_tlv;
//         const TLorentzVector met_temp  = met_tlv;
// 
//         // mass scan
//         double vw_max = 0.;
//         if(m_massLoop){
//            for(double iTopMass = 100.; iTopMass < 300.5; iTopMass += 1.){
// 
//                 dileptonTtbarRecoUtils tp_m(iTopMass, iTopMass);
//                 tp_m.execute(al_temp, l_temp, b_temp, bbar_temp, met_temp.Px(), met_temp.Py());
// 
//                 if(tp_m.getNsol()<1) continue;
//                 if(!(tp_m.getTtSol()->at(0).weight>vw_max)) continue;
// 
//                 m_NSol++;
// 
//                 vw_max=tp_m.getTtSol()->at(0).weight;
//                 m_sol.jetB    = b_temp;
//                 m_sol.jetBbar = bbar_temp;
//                 m_sol.lm  = leptonMinus_tlv;
//                 m_sol.lp  = leptonPlus_tlv;
//                 m_sol.met = met_temp;
//                 m_sol.neutrino    = tp_m.getTtSol()->at(0).neutrino;
//                 m_sol.neutrinoBar = tp_m.getTtSol()->at(0).neutrinobar;
//                 m_sol.weight = tp_m.getTtSol()->at(0).weight;
//                 m_sol.Wplus  = m_sol.lp + m_sol.neutrino;
//                 m_sol.Wminus = m_sol.lm + m_sol.neutrinoBar;
//                 m_sol.top    = m_sol.Wplus + m_sol.jetB;
//                 m_sol.topBar = m_sol.Wminus + m_sol.jetBbar;
//                 m_sol.ttbar  = m_sol.top + m_sol.topBar;
//                 m_sol.jetB_index    = j1;
//                 m_sol.jetBbar_index = j2;
//                 m_sol.ntags = nb_tag[ib];
//             }
//         }
// 
//         if(vw_max>0){
//             vect_sol.push_back(m_sol);
//         }
//     }
// 
//     setSolutions(vect_sol);
// 
//     return;
// }

// THE END

