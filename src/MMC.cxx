/*
Created:       4 June      2015
Last Updated:  8 September 2017

Tao Huang
tao.huang@cernSPAMNOT.ch
Texas A&M University

Dan Marey
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Use to run MMC algorithm, aka Heavy Mass Estimator (HME).
-- Updated for CyMiniAna framework

Included from github (https://github.com/tahuang1991/delphes/blob/master/diHiggs/)
on 7 September 2017
  SHA = bb42b6158bf70678637e9b555d8f5e771783c37b
*/
#include "diHiggs/CyMiniAna/interface/MMC.h"



MMC::MMC( configuration &cmaConfig ) :
  m_config(&cmaConfig),
  m_weightfromonshellnupt_func(false),
  m_weightfromonshellnupt_hist(false),
  m_weightfromonoffshellWmass_hist(false),
  m_weightfrombjetrescalec1c2_hist(false),
  m_onshellMarker(-1),
  m_puSample(false),
  m_iterations(0),
  m_useMET(false),
  m_WMass(80.3){
    m_isMC = m_config->isMC();
    m_pi   = TMath::Pi();
    m_generator = new TRandom3();
  }


MMC::~MMC(){
  m_file->Close();
  delete m_generator;
}



void MMC::initialize(){
    /* Set some global parameters 
     * read from : https://github.com/tahuang1991/delphes/blob/master/diHiggs/parametersconfig.txt
     */
    m_weightfromonshellnupt_func     = false;
    m_weightfromonshellnupt_hist     = true;
    m_weightfromonoffshellWmass_hist = true;

    m_iterations = 1000000;
    m_useMET     = true;
    m_file       = TFile::Open( (m_config->getAbsolutePath()+"/config/MMCRefPDF.ROOT").c_str(),"READ"); 
    m_puSample   = false;
    m_met_sigma  = m_puSample ? 25.2 : 14.8;  // PU0=14.8; PU40=25.2

    m_b1rescalefactor = 1;
    m_b2rescalefactor = 1;
    m_rescalec1 = 1;
    m_rescalec2 = 1;

    // -1.ideal case, 
    //  0. no correction; 
    //  1. simple rescale, 
    //  2. elaborate rescale, 
    //  4. simple rescale from bjet simple rescale. 
    //  5: elaborate rescale from bjet elaborate rescale 
    m_metcorrection = 5;
    m_bjetrescale   = 2;

    m_h_onoffshellWmass = getOnOffShellWMassPDF(); 
    m_h_onshellnupt     = getOnShellNuPtPDF(); 
    m_h_bjetrescalec1   = getBjetRescalec1PDF(); 
    m_h_bjetrescalec2   = getBjetRescalec2PDF(); 

    m_h_onshellnupt->Scale(    1.0/m_h_onshellnupt->GetBinContent(m_h_onshellnupt->GetMaximumBin()));
    m_h_onoffshellWmass->Scale(1.0/m_h_onoffshellWmass->GetBinContent(m_h_onoffshellWmass->GetMaximumBin()));
    m_h_bjetrescalec1->Scale(  1.0/m_h_bjetrescalec1->GetBinContent(m_h_bjetrescalec1->GetMaximumBin()));
    m_h_bjetrescalec2->Scale(  1.0/m_h_bjetrescalec2->GetBinContent(m_h_bjetrescalec2->GetMaximumBin()));

    return;
}



bool MMC::execute( unsigned int entry, DileptonReco& dilep ){
    /* Should not include any gen level information here in final version 
     * method called to run MMC method for each case
     * control 0 : take muon from onshellW  as muon from onshell W and nu_offshellW_eta = some_eta+deltaeta
     * control 1 : take muon from onshellW  as muon from onshell W and nu_offshellW_eta = some_eta-deltaeta
     * control 2 : take muon from offshellW as muon from onshell W and nu_offshellW_eta = some_eta+deltaeta
     * control 3 : take muon from offshellW as muon from onshell W and nu_offshellW_eta = some_eta-deltaeta
     */
    cma::DEBUG("MMC : execute() in MMC class ");
    m_onshellMarker = -1;    // actually determined from MC information (data = -1)

    // Access information from DileptonReco struct
    m_lep_p = dilep.lepton_pos;
    m_lep_n = dilep.lepton_neg;
    m_met   = dilep.met;
    m_jets  = dilep.jets;
    m_bbar  = m_jets.at(0).p4 + m_jets.at(1).p4; // leading two are the 'b'-jets

    m_totjets = {};          // total TLorentzVector of jets (only 2 available right now)
    for (const auto& jet : m_jets)
        m_totjets += jet.p4;


    // TLV to store relevant information
    Neutrino nu_onshellW;
    Neutrino nu_offshellW;
    TLorentzVector offshellW_lorentz;
    TLorentzVector onshellW_lorentz;
    TLorentzVector htoWW_lorentz;
    TLorentzVector htoBB_lorentz(m_bbar);
    TLorentzVector h2tohh_lorentz;
    TVector2 met_vec2(m_met.Px(), m_met.Py());

    // initialize some values
    float weight  = 1.0;
    float weight1 = 1.0;
    float weight2 = 1.0;
    float weight3 = 1.0;
    float weight4 = 1.0;

    float eta_gen   = 0;
    float phi_gen   = 0;
    float metpx_gen = 1;
    float metpy_gen = 1;
    bool validrun   = false;

    float nu_onshellW_pt(0);
    float wmass_gen(m_WMass); // initial value
    float step(0.0);
    float random01(0.0);

    int seed = time(NULL);
    m_generator->SetSeed(seed+entry);

    // histograms of HME
    std::string histname("MMC_h2Mass_"+std::to_string(entry));
    std::string histweight1name("MMC_h2Mass_weight1_"+std::to_string(entry));
    std::string histweight4name("MMC_h2Mass_weight4_"+std::to_string(entry));

    m_MMC_h2Mass        = TH1F(histname.c_str(),histname.c_str(), 3800, 200, 4000);
    m_MMC_h2Massweight1 = TH1F(histweight1name.c_str(),histweight1name.c_str(), 3800, 200, 4000);
    m_MMC_h2Massweight4 = TH1F(histweight4name.c_str(),histweight4name.c_str(), 3800, 200, 4000);

    // loop over iterations
    for (int i=0; i<m_iterations; i++){
        cma::DEBUG("MMC : execute() iteration = "+std::to_string(i)+" of "+std::to_string(m_iterations));

        eta_gen   = m_generator->Uniform(-6,6); 
        phi_gen   = m_generator->Uniform(-m_pi,m_pi);
        hmass_gen = m_generator->Gaus(125.03,0.3);
        metpx_gen = (m_metcorrection>3) ? m_generator->Gaus(0.0,m_met_sigma) : 0;
        metpy_gen = (m_metcorrection>3) ? m_generator->Gaus(0.0,m_met_sigma) : 0;

        // Generate onshell Wmass
        step      = m_generator->Uniform(-4,4);
        random01  = m_generator->Uniform(0,1);
        wmass_gen = onshellWMassRandomWalk(wmass_gen, step, random01, true); // use histogram

        // b-jet re-scaling
        if (m_bjetrescale == 1){
            float scale = 125./m_bbar.M();
            m_b1rescalefactor = scale;
            m_b2rescalefactor = scale;
            m_rescalec1       = scale; 
            m_rescalec2       = scale;
        }
        if (m_bjetrescale == 2){
            m_rescalec1 = m_h_bjetrescalec1->GetRandom();
            if (!bjetsCorrection()) continue;
        }
        if (m_bjetrescale>0){
            htoBB_lorentz = m_b1rescalefactor*(m_jets.at(0).p4)+m_b2rescalefactor*(m_jets.at(1).p4);
            if (fabs(htoBB_lorentz.M()-125.)>1){
                cma::ERROR("MMC : The htobb mass is not 'close' 125!");
                cma::ERROR("MMC :  htobb_mass     = "+std::to_string(htoBB_lorentz.M()));
                cma::ERROR("MMC :  bjetrescale b1 = "+std::to_string(m_b1rescalefactor));
                cma::ERROR("MMC :              b2 = "+std::to_string(m_b2rescalefactor));
                continue;
            }
        }

        // MET corrections
        if ((m_metcorrection-3)>0 && ((m_metcorrection-3)==m_bjetrescale || m_metcorrection==m_bjetrescale))
            met_vec2 = metCorrection();
        else if ((m_metcorrection-3)==1 or m_metcorrection==1){
            float met_corr_factor = 125./m_bbar.M();
            m_b1rescalefactor = met_corr_factor;
            m_b2rescalefactor = met_corr_factor;
            met_vec2 = metCorrection(); 
        }
        else if ((m_metcorrection-3)==2 or m_metcorrection==2){
            m_rescalec1 = m_h_bjetrescalec1->GetRandom();
            if (!bjetsCorrection()) continue;
            met_vec2 = metCorrection(); 
        }

        TVector2 met_gen = TVector2(metpx_gen, metpy_gen);
        if (m_metcorrection>3 && m_useMET) 
            met_vec2 += met_gen;

        cma::DEBUG("MMC : met_vec2 x,y = "+std::to_string(met_vec2.Px())+", "+std::to_string(met_vec2.Py()) );

        // check for solutions
        std::vector<int> solutions;   // count num of solutions case
        for (int j=0; j<4; j++){
            assignMuons(j/2);  // modify mu_on/offshellW

            nu_onshellW_pt = nu1pt_onshellW(std::make_pair(eta_gen, phi_gen), wmass_gen); 
            nu_onshellW.p4.SetPtEtaPhiM(nu_onshellW_pt, eta_gen, phi_gen,0);

            bool solution = getNuOffShellW(met_vec2,nu_onshellW,nu_offshellW,j%2, hmass_gen);

            if (solution) solutions.push_back(j);
        }

        // no go through solutions and do stuff!
        weight = 1.0/solutions.size(); // change weight if we consider possibility factor like matrix elements
        for (const auto& j : solutions){
            assignMuons(j/2);   // modify mu_on/offshellW

            nu_onshellW_pt = nu1pt_onshellW(std::make_pair(eta_gen,phi_gen), wmass_gen); 
            nu_onshellW.p4.SetPtEtaPhiM(nu_onshellW_pt, eta_gen, phi_gen,0);

            bool dummy = getNuOffShellW(met_vec2,nu_onshellW,nu_offshellW, j%2, hmass_gen);

            onshellW_lorentz  = m_mu_onshellW.p4  + nu_onshellW.p4;
            offshellW_lorentz = m_mu_offshellW.p4 + nu_offshellW.p4;
            htoWW_lorentz     = onshellW_lorentz  + offshellW_lorentz;
            h2tohh_lorentz    = htoWW_lorentz     + htoBB_lorentz;

            if (h2tohh_lorentz.M()<245 or h2tohh_lorentz.M()>3800){
                cma::ERROR("MMC : h2tohh mass to small or large");
                continue;
            }

            cma::DEBUG("MMC : hmass_gen           "+std::to_string(hmass_gen));
            cma::DEBUG("MMC : Higgs mass from MMC "+std::to_string(htoWW_lorentz.M()));

            if (m_weightfromonshellnupt_func)
                weight1 = weightfromonshellnupt(nu_onshellW_pt); 
            else if (m_weightfromonshellnupt_hist)
                weight1 = weightfromhist(m_h_onshellnupt, nu_onshellW_pt); 

            if (m_weightfromonoffshellWmass_hist){
                weight2 = weightfromhist(m_h_onoffshellWmass, wmass_gen, offshellW_lorentz.M()); 
                weight3 = weightfromhist(m_h_onoffshellWmass, wmass_gen, offshellW_lorentz.M(), false);
            }
            if (m_weightfrombjetrescalec1c2_hist) weight4 = weightfromhist(m_h_bjetrescalec2, m_rescalec2);

            weight1 *= weight;
            weight2 *= weight1;
            weight3 *= weight1;
            weight4 *= weight1;

            float h2tohh_Mass = h2tohh_lorentz.M();
            m_MMC_h2Mass.Fill(h2tohh_Mass, weight);
            m_MMC_h2Massweight1.Fill(h2tohh_Mass, weight1);
            m_MMC_h2Massweight4.Fill(h2tohh_Mass, weight4);

            cma::DEBUG("MMC : h2tohh  mass           "+std::to_string(h2tohh_Mass));
            cma::DEBUG("MMC : h2tohh  weight         "+std::to_string(weight));
            cma::DEBUG("MMC : h2tohh  weight1        "+std::to_string(weight1));
            cma::DEBUG("MMC : h2tohh  nu_onshellW_pt "+std::to_string(nu_onshellW_pt));

            validrun = true;
        } // end controls loop,(0,1,2,3)
    } // end of iteration

    cma::DEBUG("MMC : Valid run =    "+std::to_string(validrun));

    return validrun;
}



void MMC::assignMuons(int control){
    /* method called to assign muons lorenz vector 
     * execute() control/2==0, namely control=0 here, we have correct muon TLV pair
     */
    if (m_isMC && m_onshellMarker>=0){
        if (m_onshellMarker == 1){
            if (control == 0){
                m_mu_onshellW  = m_lep_p;
                m_mu_offshellW = m_lep_n; 
            }
            else if (control == 1){
                m_mu_onshellW  = m_lep_n;
                m_mu_offshellW = m_lep_p;
            }
        }
        else if (m_onshellMarker == 2){
            if (control == 0){
                m_mu_onshellW  = m_lep_n;
                m_mu_offshellW = m_lep_p;
            }
            else if (control == 1){
                m_mu_onshellW  = m_lep_p;
                m_mu_offshellW = m_lep_n;
            }
        }
    }
    else if ((m_isMC && m_onshellMarker<0) || !m_isMC ) {
        // real case, assign them randomly
        double rand = m_generator->Uniform(0,1);
        if (rand > 0.5){
            m_mu_onshellW  = m_lep_p;
            m_mu_offshellW = m_lep_n;
        }
        else{
            m_mu_onshellW  = m_lep_n;
            m_mu_offshellW = m_lep_p;
        }
    }

    return;
}



float MMC::onshellWMassRandomWalk(float x0, float step, float random){
    /* use random walk to generate random onshellW mass accroding to wmass pdf */
    float value(0.0);
    float xmin(50);
    float xmax(90);

    float x1 = x0+step;
    while (x1 > xmax || x1 < xmin){
        if (x1 > xmax) x1 = x1-xmax+xmin;
        if (x1 < xmin) x1 = xmax-(xmin-x1);
    }

    // transition probability
    float w = onshellWMassPDF(x1)/onshellWMassPDF(x0);
    if (w >= 1.00 || (w < 1.00 && random < w) ) 
        value = x1;
    else 
        value = x0;

    return value;
}  



float MMC::onshellWMassRandomWalk(float x0, float step, float random, bool hist_){
    /* use random walk to generate random onshellW mass accroding to wmass pdf */
    float value(0.0);
    float xmin(50);
    float xmax(90);

    // periodic boundary codition
    while (x0 > xmax || x0 < xmin){
        if (x0 > xmax) x0 = x0-xmax+xmin;
        if (x0 < xmin) x0 = xmax-(xmin-x0);
    }

    float x1 = x0+step;
    while (x1 > xmax || x1 < xmin){
        if (x1 > xmax) x1 = x1-xmax+xmin;
        if (x1 < xmin) x1 = xmax-(xmin-x1);
    }

    //find
    int binx0_1,binx0_2;
    int binx1_1,binx1_2;
    double bincent0_1,bincont0_1;// center and content
    double bincent1_1,bincont1_1;

    TH1F* hist = getOnShellWMassPDF();
    binx0_1 = hist->FindBin(x0);
    binx1_1 = hist->FindBin(x1);

    if ((float)hist->GetBinCenter(binx0_1) < x0){
        binx0_2 = binx0_1+1;
    }
    else {
        binx0_2 = binx0_1;
        binx0_1 = binx0_1-1;
    }

    if ((float)hist->GetBinCenter(binx1_1) < x1){
        binx1_2 = binx1_1+1;
    }
    else {
        binx1_2 = binx1_1;
        binx1_1 = binx1_1-1;
    }

    bincent0_1 = hist->GetBinCenter(binx0_1);
    bincont0_1 = hist->GetBinContent(binx0_1);
    bincent1_1 = hist->GetBinCenter(binx1_1);
    bincont1_1 = hist->GetBinContent(binx1_1);
    double w0  = (x0-bincent0_1)*(bincont0_1-hist->GetBinContent(binx0_2))/(bincent0_1-hist->GetBinCenter(binx0_2))+bincont0_1;
    double w1  = (x1-bincent1_1)*(bincont1_1-hist->GetBinContent(binx1_2))/(bincent1_1-hist->GetBinCenter(binx1_2))+bincont1_1;

    // transition probability
    double w = w1/w0;

    if (w >= 1.00 || (w<1.00 && random<(float)w) ) 
        value = x1;
    else 
        value = x0;

    return value;
}  



float MMC::weightfromhist(TH1F* hist, float x){
    /* weight solution by a histogram: hist should be scaled */
    float weight = 0.0;
    int bin1 = hist->FindBin(x);
    //first make sure that x is within range
    if (bin1 == 0 || bin1 == hist->GetNbinsX()+1) 
        return 0;

    weight = hist->Interpolate(x);

    return weight;
}



float MMC::weightfromhist(TH2F* hist, float x, float y, bool whole){
    /* weight solution by a 2d histogram : hist should be scaled */
    float value(0.0);
    float weight(0.0);
    int bin1 = hist->GetXaxis()->FindBin(x);
    int bin2 = hist->GetYaxis()->FindBin(y);

    //first make sure that x is within range
    if (bin1 == 0 || bin1 == hist->GetNbinsX()+1) return 0;
    if (bin2 == 0 || bin2 == hist->GetNbinsY()+1) return 0;

    weight = hist->GetBinContent(bin1, bin2);
    if (whole)
        value = weight;
    else
        value = (weight<0.1) ? 0 : weight/( hist->Integral(bin1,bin1,0,hist->GetNbinsY()+1) );

    return value;
}



float MMC::weightfromonshellnupt(float nupt){
    /* weight solution by nu_pt */
    float weight = 0.0;
    float max = 170;
    if (nupt<0 || nupt>125) return 0.0;

    weight = -16.925+12.4066*nupt-0.2884*std::pow(nupt,2)+0.00203*std::pow(nupt,3)+7.695e-7*std::pow(nupt,4)
             -7.2191e-8*std::pow(nupt,5)+2.499e-10*std::pow(nupt,6);
    if (weight < 0 && nupt<5) return 0.0;
    if (weight < 0)
        cma::ERROR("MMC :  nu pT "+std::to_string(nupt)+" weight "+std::to_string(weight));
    weight = weight/max;

    return weight;
}



float MMC::nu1pt_onshellW(EtaPhi nu1_etaphi, float wMass){
    /* method called to calculate pt of nuetrinos from on-shell W decay */
    float deltaeta = nu1_etaphi.first  - m_mu_onshellW.p4.Eta();
    float deltaphi = nu1_etaphi.second - m_mu_onshellW.p4.Phi();
    float nu1_pt   = wMass*wMass/(2*m_mu_onshellW.p4.Pt()*(cosh(deltaeta)-cos(deltaphi)));

    cma::DEBUG("MMC : wMass                 = "+std::to_string(wMass) );
    cma::DEBUG("MMC : m_mu_onshellW.p4.Pt() = "+std::to_string(m_mu_onshellW.p4.Pt()) );
    cma::DEBUG("MMC : cosh(deltaeta)        = "+std::to_string(cosh(deltaeta)) );
    cma::DEBUG("MMC : cos(deltaphi)         = "+std::to_string(cos(deltaphi)) );
    cma::DEBUG("MMC : nu1_pt                = "+std::to_string(nu1_pt) );

    return nu1_pt;
}



bool MMC::getNuOffShellW(const TVector2& met, 
                         const Neutrino& nu_onshellW, Neutrino& nu_offshellW,
                         const int control, const float hMass){
    /* method called to calculate lorentzvector of second nueutrino, which is from offshell W 
     * return true if we can get nu_offshellW
     * Use MET
     */
    TLorentzVector tmplorentz = m_mu_onshellW.p4 + m_mu_offshellW.p4 + nu_onshellW.p4;
    TLorentzVector tmp2lorentz;  // fake one massless lorentzvector with same pz and E
    tmp2lorentz.SetPxPyPzE(sqrt(pow(tmplorentz.Pt(),2)+pow(tmplorentz.M(),2)),0,
                           tmplorentz.Pz(),tmplorentz.E());

    cma::DEBUG("MMC : tmplorentz = "+std::to_string(m_mu_onshellW.p4.Pt()) );

    float nu_tmp_px;
    float nu_tmp_py;
    float nu_tmp_pt;

    if (m_useMET){
        nu_tmp_px = met.Px()-nu_onshellW.p4.Px();
        nu_tmp_py = met.Py()-nu_onshellW.p4.Py();
    }
    else{
        nu_tmp_px = -m_totjets.Px()-m_mu_onshellW.p4.Px()-m_mu_offshellW.p4.Px()-nu_onshellW.p4.Px();
        nu_tmp_py = -m_totjets.Py()-m_mu_onshellW.p4.Py()-m_mu_offshellW.p4.Py()-nu_onshellW.p4.Py();
    }
    TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);
    nu_tmp_pt = nu_pxpy.Mod();

    // cosh(nu_offshellW_eta-tmp2lorentz_eta)
    float chdeltaeta = (m_useMET) ? (pow(hMass,2)+2*(nu_pxpy.Px()*tmplorentz.Px()+nu_pxpy.Py()*tmplorentz.Py())-pow(tmplorentz.M(),2))/(2*tmp2lorentz.Pt()*nu_tmp_pt) :
                                    (pow(hMass,2)+pow(m_totjets.Pt(),2)-pow(tmplorentz.M(),2)-pow(tmplorentz.Pt(),2)-pow(nu_tmp_pt,2))/(2*tmp2lorentz.Pt()*nu_tmp_pt);

    if (chdeltaeta < 1.0) {
        cma::DEBUG("MMC : chdeltaeta < 1.0; returning false");
        nu_offshellW.p4.SetPtEtaPhiM(0,0,0,0);
        return false;
    }

    float nu_tmp_phi = nu_pxpy.Phi_mpi_pi(nu_pxpy.Phi());
    float deltaeta   = acosh(chdeltaeta);
    float nu_tmp_eta = (control==1) ? (tmp2lorentz.Eta()-deltaeta) : (tmp2lorentz.Eta()+deltaeta);//control = j%2 
    // should check whether deltaeta > 1

    // from simulation: |nu_offshellW_Eta|<6
    if (fabs(nu_tmp_eta) > 7) {
        cma::DEBUG("MMC : fabs(nu_tmp_eta) > 7; returning false");
        nu_offshellW.p4.SetPtEtaPhiM(0,0,0,0);
        return false;
    }
    nu_offshellW.p4.SetPtEtaPhiM(nu_tmp_pt, nu_tmp_eta, nu_tmp_phi, 0);

    cma::DEBUG("MMC : Set Higgs Mass = "+std::to_string(hMass));
    cma::DEBUG("MMC : MMC higgs mass = "+std::to_string((tmplorentz+nu_offshellW.p4).M()));

    return true; 
}



bool MMC::bjetsCorrection(){
    /* bjets correction, based on c1, calculate c2 here
     * use m_rescalec1, m_rescalec2 to correct bjets
     * c1rescale taken from pdf
     */
    bool jets_ordered_by_pt(false);
    TLorentzVector b1lorentz;
    TLorentzVector b2lorentz;
    if (m_jets.at(0).p4.Pt() > m_jets.at(1).p4.Pt()){
        b1lorentz = m_jets.at(0).p4;
        b2lorentz = m_jets.at(1).p4;
        jets_ordered_by_pt = true;
    }
    else{
        cma::WARNING("MMC : b1jet is not jet with larger pt ");
        b1lorentz = m_jets.at(0).p4;
        b2lorentz = m_jets.at(1).p4;
    }

    // x1*c2*c2+x2*c2+x3=0, slove for c2
    float x1 = b2lorentz.M2();
    float x2 = 2*m_rescalec1*(b1lorentz*b2lorentz);
    float x3 = m_rescalec1*m_rescalec1*b1lorentz.M2()-125*125;

    if (x2<0)
        cma::ERROR("MMC : Bjets lorentzvector dot productor less than 0 ");

    if ((x2*x2-4*x1*x3)<0 or x1==0)
        return false;

    m_rescalec2 = (-x2+std::sqrt(x2*x2-4*x1*x3))/(2*x1);

    if (jets_ordered_by_pt){
        m_b1rescalefactor = m_rescalec1;
        m_b2rescalefactor = m_rescalec2;
    }
    else{
        cma::WARNING("MMC : Wired b1jet is not jet with larger pt ");
        m_b2rescalefactor = m_rescalec1;
        m_b1rescalefactor = m_rescalec2;
    }

    return true;
}



TVector2 MMC::metCorrection(){
    /* Correct the MET */
    Jet jet1 = m_jets.at(0);
    Jet jet2 = m_jets.at(1);
    float metpx_correction = m_met.Px() - (m_b1rescalefactor-1)*jet1.p4.Px()
                                        - (m_b2rescalefactor-1)*jet2.p4.Px();
    float metpy_correction = m_met.Py() - (m_b1rescalefactor-1)*jet1.p4.Py()
                                        - (m_b2rescalefactor-1)*jet2.p4.Py();
    TVector2 met_vec2 = TVector2(metpx_correction, metpy_correction);

    return met_vec2;
}



TH1F MMC::getMMCh2(){
    /* Return histogram */
    return m_MMC_h2Mass;
}

TH1F MMC::getMMCh2weight1(){
    /* Return histogram */
    return m_MMC_h2Massweight1;
}

TH1F MMC::getMMCh2weight4(){
    /* Return histogram */
    return m_MMC_h2Massweight4;
}


EtaPhi MMC::generatenu1_etaphi(){
    /* method called to generate a pair (eta,phi) for nuetrino1 */
	float mean = 0;
	float rms  = 1.403;

	float eta = genEtaGuass(mean,rms);
	float phi = genPhiFlat();

    return std::make_pair(eta,phi);
}



float MMC::genEtaGuass(float mean, float rms){
    /* method called to generate eta from Gauss distribution */
    return m_generator->Gaus(mean,rms);
}



float MMC::genPhiFlat(){
    /* method called to generate phi from Flat distribution */
    return m_generator->Uniform(-m_pi,m_pi);
}



TH1F* MMC::getOnShellWMassPDF(){
    /* method called to readout TH1F onshellWmasspdf from root file */
    TH1F* onshellWmasspdf = (TH1F*)m_file->Get("onshellWmasspdf");
    return onshellWmasspdf;
}



TH1F* MMC::getOffShellWMassPDF(){
    /* method called to readout TH1F offshellWmasspdf from root file */
    TH1F* offshellWmasspdf = (TH1F*)m_file->Get("offshellWmasspdf");
    return offshellWmasspdf;
}



TH2F* MMC::getOnOffShellWMassPDF(){
    /* method called to readout TH1F onoffshellWmasspdf from root file */
    TH2F* onoffshellWmasspdf = (TH2F*)m_file->Get("onoffshellWmasspdf");
    return onoffshellWmasspdf;
}



TH1F* MMC::getOnShellNuPtPDF(){
    /* method called to readout TH1F onshellnuptpdf from root file */
    TH1F* onshellnuptpdf = (TH1F*)m_file->Get("onshellnuptpdf");
    return onshellnuptpdf;
}



TH1F* MMC::getBjetRescalec1PDF(){
        /* method called to readout TH1F recobjetrescalec1pdfPU40 from root file */
    TH1F* bjetrescalec1pdf;
    bjetrescalec1pdf = m_puSample ? (TH1F*)m_file->Get("recobjetrescalec1pdfPU40") : 
                                    (TH1F*)m_file->Get("bjetrescalec1dR4pdf");
    return bjetrescalec1pdf;
}



TH1F* MMC::getBjetRescalec2PDF(){
    /* method called to readout TH1F recobjetrescalec2pdfPU40 from root file */
    TH1F* bjetrescalec2pdf;
    bjetrescalec2pdf = m_puSample ? (TH1F*)m_file->Get("recobjetrescalec2pdfPU40") :
                                    (TH1F*)m_file->Get("bjetrescalec2dR4pdf");
    return bjetrescalec2pdf;
}



TH2F* MMC::getBjetRescalec1c2PDF(){
    /* method called to readout TH1F bjetrescalec1c2pdf from root file */
    TH2F* bjetrescalec1c2pdf = (TH2F*)m_file->Get("bjetrescalec1c2pdf");
    return bjetrescalec1c2pdf;
}



float MMC::onshellWMassPDF(float mass){
    /* method to describe onshellW mass Probability density function */
    float p0 = 7.87161e-03;
    float p1 = 1.69085;
    float p2 = 603.474 ;
    float p = exp(mass*p0+p1) + p2*exp(-0.5*((mass-80.1)*0.5)*((mass-80.1)*0.5));

    return p;
}


bool MMC::checkSolution(TLorentzVector* jetslorentz,
                        TLorentzVector* mu1lorentz,
                        TLorentzVector* mu2lorentz,
                        TLorentzVector* nu1lorentz, int control, float hMass){
    /* Method called to check whether the solution in this case exist or not 
     * ** Not used now, may be helpful later 
     */
    TLorentzVector* tmplorentz = new TLorentzVector(mu1lorentz->Px()+mu2lorentz->Px()+nu1lorentz->Px(),
                      mu1lorentz->Py()+mu2lorentz->Py()+nu1lorentz->Py(),
                      mu1lorentz->Pz()+mu2lorentz->Pz()+nu1lorentz->Pz(),
                      mu1lorentz->Energy()+mu2lorentz->Energy()+nu1lorentz->Energy());

    float nu_tmp_px;
    float nu_tmp_py;
    float nu_tmp_pt;

    nu_tmp_px = -jetslorentz->Px()-mu1lorentz->Px()-mu2lorentz->Px()-nu1lorentz->Px();
    nu_tmp_py = -jetslorentz->Py()-mu1lorentz->Py()-mu2lorentz->Py()-nu1lorentz->Py();
    TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);

    nu_tmp_pt = nu_pxpy.Mod();

    float chdeltaeta;//cosh(nu2_eta-tmp2lorenz_eta)
    TLorentzVector* tmp2lorentz = new TLorentzVector(sqrt(pow(tmplorentz->Pt(),2)+pow(tmplorentz->M(),2)),0,tmplorentz->Pz(),tmplorentz->Energy());// construct massless lorentzvector with same pz and E as tmplorentzvector

    chdeltaeta = (pow(hMass,2)+pow(jetslorentz->Pt(),2)-pow(tmplorentz->M(),2)-pow(tmplorentz->Pt(),2)-pow(nu_tmp_pt,2))/(2*tmp2lorentz->Pt()*nu_tmp_pt);

    delete tmplorentz;
    delete tmp2lorentz;

    // place the cuts we may need 
    //
    //at present if (|chdeltaeta|>1) return true; 

    return (fabs(chdeltaeta)>1);
}


// THE END
