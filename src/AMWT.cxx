/*
Created:        --
Last Updated:   22 August 2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

@source: https://github.com/cms-ljmet/topDileptonMassAnalysis

AMWT class
    - Reconstruct the dilepton ttbar system using analytical matrix weighting technique

*/
#include "cms-ttbarAC/CyMiniAna/interface/AMWT.h"


// Constructor
AMWT::AMWT(configuration &cmaConfig) : 
  m_config(&cmaConfig),
  m_max_weight(-1e30),
  m_ptResol(nullptr),
  m_etaResol(nullptr),
  m_phiResol(nullptr){
    m_sqrt_s       = m_config->beamEnergy();  // 13000.;
    m_topQuarkMass = m_config->topQuarkMass();
    m_bQuarkMass   = m_config->bQuarkMass();
    m_WMass        = m_config->WMass();
    m_NJetSmear    = m_config->NJetSmear();   // 500
    m_doGaussian   = false; // not used: m_config->doGaussianResolution();  // false

    m_top     = {};
    m_antitop = {};

    // put in histogrammer
    m_NMassPts    = m_config->NMassPoints(); // 1
    m_rangeLow    = m_config->massMin();     // 173;
    m_rangeHigh   = m_rangeLow+m_NMassPts;   // 175

    // Mass Solver
    m_massSolver = new MassSolver(m_config);

    // Random numbers
    m_rand3 = new TRandom3();

    // EventShape calculation for neutrino "likelihood"
    m_eventShape = new TF2("landau2D","[0]*TMath::Landau(x,[1],[2],0)*TMath::Landau(y,[3],[4],0)",0,500,0,500);
    m_eventShape->SetParameters(30.7137,56.2880,23.0744,59.1015,24.9145);
}


AMWT::~AMWT() {
    delete m_massSolver;
    delete m_rand3;
    delete m_eventShape;
/*  // The following aren't used:
    delete m_pdf;
    delete m_ptResol;
    delete m_etaResol;
    delete m_phiResol; */
}



void AMWT::initialize(){
    /* Initialize many variables / objects */
    cma::DEBUG("AMWT : Initialize() ");

/*    // PDF initialization
    // https://alcaraz.web.cern.ch/alcaraz/PDF_reweight.txt (@ int rewtest())
    std::string pdfset("cteq6l1"); //_0000.dat");

    const char* lhapdf_dir = getenv("LHAPDF_DATA_PATH");
    std::string lhapdf_name = lhapdf_dir;
    lhapdf_name += "/"+pdfset;
    cma::INFO("AMWT : LHAPDF file "+lhapdf_name);

//    cma::check_file(lhapdf_name);
    //LHAPDF::setPDFPath(lhapdf_dir);
    //LHAPDF::initPDFSet(lhapdf_name);
    //LHAPDF::initPDFSet(pdfset);
    m_pdf = LHAPDF::mkPDF(pdfset, 0);
*/
    // Jet resolution functions
    //JRDatabase/textFiles/Spring16_25nsV6_MC/Spring16_25nsV6_MC_PtResolution_AK4PF.txt
    const char* jr_dir = getenv("JRDATABASE");  // /home/demarley/diHiggs/CMSSW_8_0_24_patch1/src/diHiggs/CyMiniAna
    std::string jrdatabase = jr_dir;
    std::string resolutionFilePt  = jrdatabase+"/textFiles/Spring16_25nsV6_MC/Spring16_25nsV6_MC_PtResolution_AK4PF.txt";
    std::string resolutionFilePhi = jrdatabase+"/textFiles/Spring16_25nsV6_MC/Spring16_25nsV6_MC_PhiResolution_AK4PF.txt";
    std::string scaleFactorFile   = jrdatabase+"/textFiles/Spring16_25nsV6_MC/Spring16_25nsV6_MC_SF_AK4PF.txt";

    m_resolution_from_file_pt.reset(  new JME::JetResolution(resolutionFilePt) );
    m_resolution_from_file_phi.reset( new JME::JetResolution(resolutionFilePhi) );   //std::unique_ptr<JME::JetResolution>(new JME::JetResolution(resolutionFilePhi)) );
    m_scale_factor_from_file.reset(   new JME::JetResolutionScaleFactor(scaleFactorFile) );   //std::unique_ptr<JME::JetResolution>(new JME::JetResolutionScaleFactor(scaleFactorFile)) );

    m_resolution_pt  = *m_resolution_from_file_pt;
    m_resolution_phi = *m_resolution_from_file_phi;
    m_resolution_sf  = *m_scale_factor_from_file;

/*
    // Original implementation:
    std::string jr_name(name+"/JR_Standalone/JetMETObjects/data/");
    std::string JR_Standalone_Path(jr_name);
    cma::INFO("JetResolution path "+JR_Standalone_Path);

    std::string ptFileName  = JR_Standalone_Path+"Spring10_PtResolution_AK5PF.txt";
    std::string etaFileName = JR_Standalone_Path+"Spring10_EtaResolution_AK5PF.txt";
    std::string phiFileName = JR_Standalone_Path+"Spring10_PhiResolution_AK5PF.txt";

    cma::check_file(ptFileName);
    cma::check_file(etaFileName);
    cma::check_file(phiFileName);

    cma::INFO("pT filename "+ptFileName);
    m_ptResol  = new JetResolution(ptFileName,m_doGaussian);

    cma::INFO("eta filename "+etaFileName);
    m_etaResol = new JetResolution(etaFileName,m_doGaussian);

    cma::INFO("phi filename "+phiFileName);
    m_phiResol = new JetResolution(phiFileName,m_doGaussian);
*/
    return;
}


std::map<std::string,Top> AMWT::findMass( DileptonReco& dilep) {
    /* Get the ttbar system */
    m_max_weight    = -1e30;                // reset to really large negative number (such that any solution will be selected)
    m_max_weight_ES = -1e30;                // reset to really large negative number (such that any solution will be selected)

    std::map<std::string,Top> ttbar;     // store Top candidates
    m_top     = {};
    m_antitop = {};
    m_top.weight        = m_max_weight;
    m_antitop.weight    = m_max_weight;
    m_top.weight_ES     = m_max_weight_ES;
    m_antitop.weight_ES = m_max_weight_ES;


    Lepton lepton_p = dilep.lepton_pos;
    Lepton lepton_n = dilep.lepton_neg;
    TVector2 met_sh = dilep.met;
    std::vector<Jet> jets_by_pt_sh = dilep.jets;


    // ** Begin JetMET smearing loop **
    TVector2 met_sm;
    std::vector<Jet> jets_by_pt_sm;
    double top_iter, sum_weight;

    // Loop over JetMET smearing
    bool jetmet_smear = (m_NJetSmear == 1) ? false : true;
    for (int iSM_JMT=0; iSM_JMT<m_NJetSmear; iSM_JMT++){
        cma::DEBUG("AMWT : Smearing iteration "+std::to_string(iSM_JMT)+" of "+std::to_string(m_NJetSmear));

        smear_JetMET(jets_by_pt_sh, met_sh, jets_by_pt_sm, met_sm,
                     (lepton_p.p4+lepton_n.p4), jetmet_smear);


        // Loop over mass points (only 1 for now)
        for (int i_t=0; i_t<m_NMassPts; i_t++){
            top_iter = m_rangeLow + i_t;

            for (unsigned int j_id=0; j_id<2; j_id++){

                // change combination of jets+lepton
                std::vector<Jet> tmp_jets;
                tmp_jets.resize(2);
                if (j_id==0){
                    tmp_jets[0] = jets_by_pt_sm[0];
                    tmp_jets[1] = jets_by_pt_sm[1];
                }
                else{
                    tmp_jets[0] = jets_by_pt_sm[1];
                    tmp_jets[1] = jets_by_pt_sm[0];
                }

                // vectors to hold neutrino solutions
                std::vector<Neutrino> nu1;
                std::vector<Neutrino> nu2;
                nu1.clear();
                nu2.clear();

                // Get the neutrino solutions
                m_massSolver->solve( met_sm, tmp_jets[0], tmp_jets[1], 
                                     lepton_p, lepton_n, top_iter, nu1, nu2 );
                cma::DEBUG("AMWT : Solutions for mtop "+std::to_string(top_iter)+" : "+std::to_string(nu1.size()));

                // For each solution, determine the weight ("likelihood" the reconstruct is correct)
                for (unsigned int ui=0, size=nu1.size(); ui<size; ui++){

                    // define 4-vector for top quarks
                    Top top;
                    top.lepton   = lepton_p;
                    top.neutrino = nu1[ui];
                    top.jet      = tmp_jets[0];
                    top.set_p4();

                    Top antitop;
                    antitop.lepton   = lepton_n;
                    antitop.neutrino = nu2[ui];
                    antitop.jet      = tmp_jets[1];
                    antitop.set_p4();

                    double s_weight = get_weight(top, antitop, top_iter);

                    if (m_max_weight < s_weight){
                        // the current solution is greater than the others!
                        // save the information from here
                        cma::DEBUG("AMWT :  - "+std::to_string(s_weight)+" :: "+std::to_string(m_max_weight));
                        cma::DEBUG("AMWT :  - "+std::to_string(top.p4.Pt()));
                        m_max_weight = s_weight;

                        m_top        = top;
                        m_top.weight = m_max_weight;
                        m_antitop        = antitop;
                        m_antitop.weight = m_max_weight;
                    }

                    // weight from EventShape here: https://github.com/cms-sw/cmssw/blob/master/TopQuarkAnalysis/TopKinFitter/src/TtFullLepKinSolver.cc
                    double s_weight_ES = get_weight_ES(nu2[ui],nu1[ui]);

                    if (m_max_weight_ES < s_weight_ES){
                        // the current solution is greater than the others!
                        // save the information from here
                        m_max_weight_ES = s_weight_ES;

                        m_top.weight_ES     = m_max_weight_ES;
                        m_antitop.weight_ES = m_max_weight_ES;
                    }
                } // end loop over neutrino solutions
            } // end loop over jet combinations
        } // End loop over mass values
    } // End JetMet smearing loop


/*
    for (unsigned int i = 0; i < vPtRes.size(); i++){
        delete vPtRes[i];
        delete vEtaRes[i];
        delete vPhiRes[i];
    }
*/

    ttbar["top"]     = m_top;
    ttbar["antitop"] = m_antitop;

    return ttbar;
}


double AMWT::get_weight_ES( Neutrino& nu_1, Neutrino& nu_2 ) const{
    /* Get the weight using EventShape parameterization
     * https://github.com/cms-sw/cmssw/blob/master/TopQuarkAnalysis/TopKinFitter/src/TtFullLepKinSolver.cc#L309-L312
     */
    double prob = m_eventShape->Eval(nu_1.p4.E(),nu_2.p4.E());
    return prob;
}


double AMWT::get_weight(const Top& t,  const Top& tbar, const double top_mass) const {
    /* Get the AMWT weight */
    double prob_dalitz(1.0);
    prob_dalitz *= get_dalitz_prob( t );
    prob_dalitz *= get_dalitz_prob( tbar );

    cma::DEBUG("AMWT : probability dalitz = "+std::to_string( prob_dalitz ));

    return prob_dalitz;
}



double AMWT::get_dalitz_prob(const Top& top) const{
    /* Get the Dalitz probability */
    double mte = top.lepton.p4.Dot( top.p4 );
    double mt2 = pow(top.p4.M(),2);

    cma::DEBUG("AMWT : MTE = "+std::to_string(mte));

    double mb2 = pow(m_bQuarkMass,2);
    double mw2 = pow(m_WMass,2);
    double mt2_mb2 = mt2 - mb2;

    double result = 4. * mte * ( mt2-mb2 - 2.*mte ) / 
                    ( pow(mt2_mb2,2) + mw2*(mt2+mb2) - 2.*pow(mw2,2) );

    return result;
}


void AMWT::smear_JetMET(const std::vector<Jet>& orig_jets, const TVector2& orig_met,
                        std::vector<Jet>& smear_jets, TVector2& smear_met,
                        const TLorentzVector& lep_sum, bool smear) const {
    if (smear){
        smear_JetMET(orig_jets, orig_met, smear_jets, smear_met, lep_sum);
    }
    else{
        smear_met  = orig_met;
        smear_jets = orig_jets;
    }

    return;
}

void AMWT::smear_JetMET(const std::vector<Jet>& orig_jets, const TVector2& orig_met,
                        std::vector<Jet>& smear_jets, TVector2& smear_met, const TLorentzVector& lep_sum) const {
    /* Implementation of jet smearing using JER from updated sources 
     *   Smear jet pt using a random gaussian variation
     *   https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L239-L250
     *   https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
     *   jet energy density: https://arxiv.org/pdf/1607.03663.pdf; section 4.1
     */
    smear_jets.clear();

    double sum_jpx(0);
    double sum_jpy(0);
    double sum_jpx_sm(0);
    double sum_jpy_sm(0);

    Jet v_temp_jet;
    TLorentzVector v_temp;

    for (unsigned int sui=0, size=orig_jets.size(); sui<size; sui++){ 
        // get the resolution
        JME::JetParameters jet_params;
        jet_params.setJetPt( orig_jets.at(sui).p4.Pt());
        jet_params.setJetEta(orig_jets.at(sui).p4.Eta());
        jet_params.setRho(   orig_jets.at(sui).rho);       // offset jet energy density, 1 value per event
        double res_pt  = m_resolution_pt.getResolution(jet_params);
        double res_phi = m_resolution_phi.getResolution(jet_params);
        double jer_sf  = m_resolution_sf.getScaleFactor(jet_params, Variation::NOMINAL);

        // turn resolution into gaussian & smear
        double sigma_pt  = res_pt  * std::sqrt(jer_sf * jer_sf - 1);
        double sigma_phi = res_phi * std::sqrt(jer_sf * jer_sf - 1);
        cma::DEBUG("AMWT : Gaussian width pT:  "+std::to_string(sigma_pt) );
        cma::DEBUG("AMWT : Gaussian width phi: "+std::to_string(sigma_phi) );

        double smearFactor_pt  = 1. + m_rand3->Gaus(0.0,sigma_pt);
        double smearFactor_phi = 1. + m_rand3->Gaus(0.0,sigma_phi);

        // build a copy of the original jet
        v_temp.SetPtEtaPhiE(orig_jets.at(sui).p4.Pt(),orig_jets.at(sui).p4.Eta(),orig_jets.at(sui).p4.Phi(),orig_jets.at(sui).p4.E());
        // smear energy and phi
        v_temp *= smearFactor_pt;                                      // smear the energy
        v_temp.SetPhi( orig_jets.at(sui).p4.Phi() * smearFactor_phi ); // smear phi

        sum_jpx += orig_jets.at(sui).p4.Px();  // original p_x values
        sum_jpy += orig_jets.at(sui).p4.Py();  // original p_y values

        sum_jpx_sm += v_temp.Px();             // smeared p_x values
        sum_jpy_sm += v_temp.Py();             // smeared p_y values

        v_temp_jet.p4 = v_temp;
        smear_jets.push_back(v_temp_jet);
    }

    double unclust_metx = orig_met.Px() + sum_jpx + lep_sum.Px();
    double unclust_mety = orig_met.Py() + sum_jpy + lep_sum.Py();

    // 10% resolution
    double unclust_metx_sm = unclust_metx * (1 + 0.1*m_rand3->Gaus());
    double unclust_mety_sm = unclust_mety * (1 + 0.1*m_rand3->Gaus());

    smear_met.Set(orig_met.Px() + sum_jpx - unclust_metx - sum_jpx_sm + unclust_metx_sm,
                  orig_met.Py() + sum_jpy - unclust_mety - sum_jpy_sm + unclust_mety_sm);

    return;
}


void AMWT::smear_JetMET(const std::vector<Jet>& orig_jets, const TVector2& orig_met,
                        std::vector<Jet>& smear_jets, TVector2& smear_met, const TLorentzVector& lep_sum, int old) const {
    /* Original implementation from @source (https://github.com/cms-ljmet/topDileptonMassAnalysis) 
     * @param old    Used to distinguish this function from the new implementation (this function is from the original code)
     */
    std::vector<TF1*> vPtRes;
    std::vector<TF1*> vEtaRes;
    std::vector<TF1*> vPhiRes;

    if (m_NJetSmear != 1 && old>0){
        for (unsigned ijr = 0; ijr < orig_jets.size(); ijr++){

            TF1* fPtResol  = (TF1*) m_ptResol->resolutionEtaPt(orig_jets.at(ijr).p4.Eta(),  orig_jets.at(ijr).p4.Pt())->Clone();
            TF1* fEtaResol = (TF1*) m_etaResol->resolutionEtaPt(orig_jets.at(ijr).p4.Eta(), orig_jets.at(ijr).p4.Pt())->Clone();
            TF1* fPhiResol = (TF1*) m_phiResol->resolutionEtaPt(orig_jets.at(ijr).p4.Eta(), orig_jets.at(ijr).p4.Pt())->Clone();

            vPtRes.push_back(fPtResol);
            vEtaRes.push_back(fEtaResol);
            vPhiRes.push_back(fPhiResol);
        }
    }

    smear_jets.clear();

    double sum_jpx(0);
    double sum_jpy(0);

    double sum_jpx_sm(0);
    double sum_jpy_sm(0);

    double Pt_sm;
    double Eta_sm;
    double Phi_sm;

    Jet v_temp_jet;
    TLorentzVector v_temp;
    Long64_t iseed = (Long64_t) 10e10;

    // loop over un-smeared jets
    for (unsigned int sui=0, size=orig_jets.size(); sui<size; sui++){

        // get smeared values for pT, eta, and phi
        gRandom->SetSeed(m_rand3->Integer(iseed));
        Pt_sm  = orig_jets.at(sui).p4.Pt()  * vPtRes.at(sui)->GetRandom();    // pT
        gRandom->SetSeed(m_rand3->Integer(iseed));
        Eta_sm = orig_jets.at(sui).p4.Eta() + vEtaRes.at(sui)->GetRandom();   // eta
        gRandom->SetSeed(m_rand3->Integer(iseed));
        Phi_sm = orig_jets.at(sui).p4.Phi() + vPhiRes.at(sui)->GetRandom();   // phi

        v_temp.SetPtEtaPhiM(Pt_sm, Eta_sm, Phi_sm, orig_jets.at(sui).p4.M()); // smeared TLorentzVector
   
        sum_jpx += orig_jets.at(sui).p4.Px();  // original p_x values
        sum_jpy += orig_jets.at(sui).p4.Py();  // original p_y values

        sum_jpx_sm += v_temp.Px();             // smeared p_x values
        sum_jpy_sm += v_temp.Py();             // smeared p_y values

        v_temp_jet.p4 = v_temp;                // smeared "Jet"
        smear_jets.push_back(v_temp_jet);
    }

    double unclust_metx = orig_met.Px() + sum_jpx + lep_sum.Px();
    double unclust_mety = orig_met.Py() + sum_jpy + lep_sum.Py();

    // 10% resolution
    double unclust_metx_sm = unclust_metx * (1 + 0.1*m_rand3->Gaus());
    double unclust_mety_sm = unclust_mety * (1 + 0.1*m_rand3->Gaus());

    smear_met.Set(orig_met.Px() + sum_jpx - unclust_metx - sum_jpx_sm + unclust_metx_sm, 
                  orig_met.Py() + sum_jpy - unclust_mety - sum_jpy_sm + unclust_mety_sm);

    return;
}

// THE END //
