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

*/
#include "cms-ttbarAC/CyMiniAna/interface/dileptonTtbarReco.h"



dileptonTtbarReco::dileptonTtbarReco(configuration& cmaConfig,
                                     const Era::Era era,     const int minNumberOfBtags, 
                                     const bool preferBtags, const bool massLoop):
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

    // Read all histograms for smearing objects
    if (!m_massLoop) loadData();
  }

dileptonTtbarReco::~dileptonTtbarReco() {
    delete m_r3;
}


dileptonTtbarRecoSolution dileptonTtbarReco::execute(const ttbarDilepton& ttSystem){
    /* Pass event information to this class and solve for the reconstructed system */
    dileptonTtbarRecoSolution result;

    // Check if minimum number of objects for any solution is present
    const int numberOfJets(ttSystem.jets.size());
    const int numberOfBtags(ttSystem.bjets.size());
    if(numberOfJets<2 || numberOfBtags<m_minNumberOfBtags) 
        return result;

    Lepton lepton     = ttSystem.lepton_neg;
    Lepton antiLepton = ttSystem.lepton_pos;

    // -- Find solutions with 2 b-tagged jets
    if (numberOfBtags>=2){
        for (unsigned int i_index=0; i_index<numberOfBtags-1; ++i_index) {
            for (unsigned int j_index=i_index+1; i_index<numberOfBtags; ++j_index) {
                Jet jet1 = ttSystem.bjets.at(i_index);
                Jet jet2 = ttSystem.bjets.at(j_index);

                ttbarDilepton tt_tmp1;
                tt_tmp1.lepton_pos = antiLepton;
                tt_tmp1.lepton_neg = lepton;
                tt_tmp1.bJet    = jet1;
                tt_tmp1.bbarJet = jet2;
                tt_tmp1.met     = ttSystem.met;

                ttbarDilepton tt_tmp2;           // same as tt_tmp1 with the jets switched
                tt_tmp2.lepton_pos = antiLepton;
                tt_tmp2.lepton_neg = lepton;
                tt_tmp2.bJet    = jet2;
                tt_tmp2.bbarJet = jet1;
                tt_tmp2.met     = ttSystem.met;

                const std::vector<ttbarDilepton> solution1( getSolutions(tt_tmp1,2) );
                const std::vector<ttbarDilepton> solution2( getSolutions(tt_tmp2,2) );
                result.Add(solution1);
                result.Add(solution2);
            }
        }
    }

    if(m_preferBtags && result.numberOfSolutions(2)) 
        return result;

    if(m_minNumberOfBtags > 1) 
        // require at least 2 b-tags, but not enough solutions; exit.
        return result;


    // Collect non b-tagged jets for reconstruction
    const int numberOfNonBtags(numberOfJets - numberOfBtags);
    if(numberOfNonBtags < 1) 
        return result;

    if (numberOfBtags==1){
        Jet bJet = ttSystem.bjets.at(0);

        // -- Find solutions with 1 b-tagged jet
        for(const int jet : ttSystem.jets) {
            if (jet.isbtagged.at(m_btag_wp)) continue;  // skip the b-tagged jet

            ttbarDilepton tt_tmp1;
            tt_tmp1.lepton_pos = antiLepton;
            tt_tmp1.lepton_neg = lepton;
            tt_tmp1.bJet    = bJet;
            tt_tmp1.bbarJet = jet;
            tt_tmp1.met     = ttSystem.met;

            ttbarDilepton tt_tmp2;           // same as tt_tmp1 with the jets switched
            tt_tmp2.lepton_pos = antiLepton;
            tt_tmp2.lepton_neg = lepton;
            tt_tmp2.bJet    = jet;
            tt_tmp2.bbarJet = bJet;
            tt_tmp2.met     = ttSystem.met;

            const std::vector<ttbarDilepton> solution1( getSolutions(tt_tmp1,1) );
            const std::vector<ttbarDilepton> solution2( getSolutions(tt_tmp2,1) );
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
        for (unsigned int i_index=0; i_index<numberOfJets-1; ++i_index) {
            for (unsigned int j_index=i_index+1; i_index<numberOfJets; ++j_index) {
                Jet jet1 = ttSystem.jets.at(i_index);
                Jet jet2 = ttSystem.jets.at(j_index);

                ttbarDilepton tt_tmp1;
                tt_tmp1.lepton_pos = antiLepton;
                tt_tmp1.lepton_neg = lepton;
                tt_tmp1.bJet    = jet1;
                tt_tmp1.bbarJet = jet2;
                tt_tmp1.met     = ttSystem.met;

                ttbarDilepton tt_tmp2;           // same as tt_tmp1 with the jets switched
                tt_tmp2.lepton_pos = antiLepton;
                tt_tmp2.lepton_neg = lepton;
                tt_tmp2.bJet    = jet2;
                tt_tmp2.bbarJet = jet1;
                tt_tmp2.met     = ttSystem.met;

                const std::vector<ttbarDilepton> solution1( getSolutions(tt_tmp1,0) );
                const std::vector<ttbarDilepton> solution2( getSolutions(tt_tmp2,0) );
                result.Add(solution1);
                result.Add(solution2);
            }
        }
    }

    return result;
}



std::vector<ttbarDilepton> dileptonTtbarReco::getSolutions(ttbarDilepton& singleTtbarSystem, const int numberOfBtags);
    /* Obtain the solutions from Utils*/
    std::vector<dileptonTtbarRecoSolution> result;

    // assign objects to the (simple) struct of information
    const DileptonReco ttSystem;
    ttSystem.lepton_pos = singleTtbarSystem.lepton_pos;
    ttSystem.lepton_neg = singleTtbarSystem.lepton_neg;
    ttSystem.bJet       = singleTtbarSystem.bJet;
    ttSystem.bbarJet    = singleTtbarSystem.bbarJet;
    ttSystem.met        = singleTtbarSystem.met;

    if(m_massLoop){
        // resolve the ttbar kinematics assuming different values for top mass
        double neutrinoWeightMax(0.);

        for(double iTopMass=m_rangeLow; iTopMass<m_rangeHigh; iTopMass+=1){
            dileptonTtbarRecoUtils tp_m(iTopMass, iTopMass);
            tp_m.execute(ttSystem);

            if(tp_m.getNsol() < 1)
                continue;

            ttbarDilepton tp_m_soln = tp_m.getTtSol();
            const double& weight    = tp_m_soln.at(0).weight;
            if(weight <= neutrinoWeightMax)
                continue;

            neutrinoWeightMax = weight;
            std::map<ttbarDilepton::WeightType, double> tmp_weight;
            tmp_weight[ttbarDilepton::defaultForMethod] = weight;
            tmp_weight[ttbarDilepton::neutrinoEnergy]   = weight;

            const Neutrino neutrino     = tp_m_soln.at(0).neutrino;
            const Neutrino antiNeutrino = tp_m_soln.at(0).neutrinobar;

            Top top;
            top.lepton   = ttSystem.lepton_pos;
            top.neutrino = neutrino;
            top.jet      = ttSystem.bJet;
            top.set_p4_lep();

            Top antiTop;
            antiTop.lepton   = ttSystem.lepton_neg;
            antiTop.neutrino = antiNeutrino;
            antiTop.jet      = ttSystem.bbarJet;
            antiTop.set_p4_lep();

            tp_m_soln.at(0).mapOfWeights = tmp_weight;
            tp_m_soln.at(0).top    = top;
            tp_m_soln.at(0).topBar = antiTop;
            tp_m_soln.at(0).set_ttbar();
            tp_m_soln.at(0).set_W();

            dileptonTtbarRecoSolution this_soln;
            this_soln.Add( tp_m_soln.at(0) );

            result.push_back(this_soln);
        }
    }
    else{
        // solve for the ttbar kinematics after smearing inputs
        std::vector<ttbarDilepton> smeared_solutions,
        std::vector<double> smeared_weights,
        const bool smearSolution = solutionSmearing(smeared_solutions,smeared_weights,ttSystem);

        // Check if there is also a solution without smearing
        dileptonTtbarRecoUtils tp_NOsm();
        tp_NOsm.execute(ttSystem);
        const bool noSmearSolution = (tp_NOsm.getNsol()>0);

        if(smearSolution && !noSmearSolution){
            // Solution with smearing
            const double& weight = std::accumulate(smeared_weights.begin(),smeared_weights.end(),0.0);

            std::map<ttbarDilepton::WeightType, double> tmp_weight;
            tmp_weight[ttbarDilepton::defaultForMethod]         = weight;
            tmp_weight[ttbarDilepton::averagedSumSmearings_mlb] = weight;

            // save the solution
            // solution.mapOfWeights = tmp_weight;
            dileptonTtbarRecoSolution this_soln;
            this_soln.Add( smeared_solutions );

            result.push_back(this_soln);
        }
        else if (!smearSolution && noSmearSolution) {
            // Solution without smearing
            dileptonTtbarRecoSolution this_soln;
            this_soln.Add( tp_NOsm );

            result.push_back(this_soln);
        }
    }

    return result;
}


bool dileptonTtbarReco::getSmearedSolutions(std::vector<ttbarDilepton>& smeared_solutions,
                                            std::vector<double>& smeared_weights,
                                            DileptonReco& ttSystem) const {
    /* Determine solution to ttbar reconstruction after smearing the inputs 
    
    ++++ Don't know what this is for, removing it ++++
    if((antiLepton.p4+jet1.p4).M()>180. || (lepton.p4+jet2.p4).M()>180.)
        // protection against large masses?
        return false;
    */
    const unsigned int seed = cma::setRandomNumberSeeds(lepton,antiLepton,jet1,jet2);
    m_r3->SetSeed(seed);  // Set random number generator seeds based on this arrangement

    smeared_solutions.clear();

    bool isHaveSol(false);        // true if at least 1 smeared iteration produces a solution
    double err(0.001);            // tolerance for angle rotation
    dileptonTtbarRecoUtils util;  // object responsible for finding the ttbar solution

    Lepton lepton     = ttSystem.lepton_neg;
    Lepton antiLepton = ttSystem.lepton_pos;
    Jet jet1 = ttSystem.bJet;
    Jet jet2 = ttSystem.bbarJet;
    TVector2 met = ttSystem.met;

    TVector3 vX_reco = -jet1.p4.Vect() - jet2.p4.Vect() - lepton.p4.Vect() - antiLepton.p4.Vect() - met.Vect();

    for(int sm=0,size=m_config->NJetSmear(); sm<size; ++sm){
        TLorentzVector b_sm;
        TLorentzVector bbar_sm;
        TLorentzVector met_sm;
        TLorentzVector l_sm;
        TLorentzVector al_sm;

        // jets energy smearing
        double fB    = m_h_jetEres->GetRandom();
        double xB    = sqrt((fB*fB*pow(jet1.p4.E(),2)-jet1.p4.M2())/(pow(jet1.p4.P(),2)));
        double fBbar = m_h_jetEres->GetRandom();
        double xBbar = sqrt((fBbar*fBbar*pow(jet2.p4.E(),2)-jet2.p4.M2())/(pow(jet2.p4.P(),2)));

        // b-jet angle smearing
        b_sm.SetXYZT(jet1.p4.Px()*xB,jet1.p4.Py()*xB,jet1.p4.Pz()*xB,jet1.p4.E()*fB);
        util.angle_rot(m_h_jetAngleRes->GetRandom(),err,b_sm,b_sm);

        // bbar-jet angle smearing
        bbar_sm.SetXYZT(jet2.p4.Px()*xBbar,jet2.p4.Py()*xBbar,jet2.p4.Pz()*xBbar,jet2.p4.E()*fBbar);
        util.angle_rot(m_h_jetAngleRes->GetRandom(),err,bbar_sm,bbar_sm);


        // leptons energy smearing
        double fL  = m_h_lepEres->GetRandom();
        double xL  = sqrt((fL*fL*pow(lepton.p4.E(),2)-lepton.p4.M2())/(pow(lepton.p4.P(),2)));
        double faL = m_h_lepEres->GetRandom();
        double xaL = sqrt((faL*faL*pow(antiLepton.p4.E(),2)-antiLepton.p4.M2())/( pow(antiLepton.p4.P(),2)));

        // lepton angle smearing
        l_sm.SetXYZT(lepton.p4.Px()*xL,lepton.p4.Py()*xL,lepton.p4.Pz()*xL,lepton.p4.E()*fL);
        util.angle_rot(m_h_lepAngleRes->GetRandom(),err,l_sm,l_sm);

        // anti lepton angle smearing
        al_sm.SetXYZT(al_sm.Px()*xaL,al_sm.Py()*xaL,al_sm.Pz()*xaL,al_sm.E()*faL);
        util.angle_rot(m_h_lepAngleRes->GetRandom(),err,al_sm,al_sm);

        // Get smeared MET
        TVector3 metV3_sm = -b_sm.Vect()-bbar_sm.Vect()-l_sm.Vect()-al_sm.Vect()-vX_reco;
        met_sm.SetXYZM(metV3_sm.Px(),metV3_sm.Py(),0,0);


        // Get new solution!
        dileptonTtbarRecoUtils tp_sm(m_h_wmass->GetRandom(),m_h_wmass->GetRandom());
        tp_sm.execute(al_sm, l_sm, b_sm, bbar_sm, met_sm);

        if(tp_sm.getNsol()>0) {
            isHaveSol = true;    // at least 1 smeared solution gives a result
            double mbl_weight = m_h_mbl_w->GetBinContent(m_h_mbl_w->FindBin((al_sm+b_sm).M())) * 
                                                         m_h_mbl_w->GetBinContent(m_h_mbl_w->FindBin((l_sm+bbar_sm).M())) / 
                                                         100000000; // where does this come from?
            smeared_solutions.push_back( tp_sm.getTtSol()->at(0) );
            smeared_weights.push_back( mbl_weight );
            // use `mbl_weight` for "sumWeights", meanVect, 
        }
    }

    return isHaveSol;
}


void dileptonTtbarReco::setSolutions() {
    setSolutions(m_sols);
}


void dileptonTtbarReco::setSolutions(std::vector<ttbarDilepton> sols) {
    m_NSol = (int)(sols.size());

    if(m_NSol>0){
        std::nth_element(begin(sols), begin(sols), end(sols),
                         [](const ttbarDilepton& a, const ttbarDilepton& b){
                            return b.ntags < a.ntags || (b.ntags == a.ntags && b.weight < a.weight);
                         });
        m_sol = sols[0];
    }

    return;
}


int dileptonTtbarReco::getNSol() const {
    return m_NSol;
}

ttbarDilepton dileptonTtbarReco::getSol() const {
    return m_sol;
}

std::vector<ttbarDilepton> dileptonTtbarReco::getSols() const {
    return m_sols;
}



void dileptonTtbarReco::loadData() {
    /* Load data for smearing 
         jet & lepton resolutions
         mbl & W mass
    */
    std::string data_path = m_config->getAbsolutePath()+"/data";

    if (m_era == Era::run2_13tev_25ns || m_era == Era::run2_13tev_2015_25ns)
        data_path.Append("/KinReco_input_2015_Run2CD_25ns_76X.root");
    else if (m_era == Era::run2_13tev_2016_25ns)
        data_path.Append("/KinReco_input_2016_Run2BtoH_25ns_80X_v036_dilep_sel_25July2017.root");
    else{
        cma::ERROR("DILEPTON : loadData() : Era is not supported");
        std::exit(357);
    }

    TFile dataFile(data_path);

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
//     std::vector<ttbarDilepton> vect_sol;
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
