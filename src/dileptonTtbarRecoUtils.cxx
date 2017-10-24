/*
Created:        --
Last Updated:   15 October  2017

Dan Marley
daniel.p4.Edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

 Methods for dileptonTtbarReco
  Class for dilepton ttbar reconstruction.
  Imported from
   https://gitlab.cern.ch/cms-desy-top/TopAnalysis/blob/master/
           Configuration/analysis/common/src/KinematicReconstruction_LSroutines.cc
  on 15 October 2017

  This is similar to the existing src/MassSolver.cxx class

*/
#include "cms-ttbarAC/CyMiniAna/interface/dileptonTtbarRecoUtils.h"


dileptonTtbarRecoUtils::dileptonTtbarRecoUtils(configuration &cmaConfig,
                                               const double& mass_top, const double& mass_topbar, 
                                               const double& mass_b,   const double& mass_bbar, 
                                               const double& mass_Wp,  const double& mass_Wm, 
                                               const double& mass_al,  const double& mass_l) :
  m_config(&cmaConfig){
    m_mt     = mass_top;
    m_mtbar  = mass_topbar;
    m_mb     = mass_b;
    m_mbbar  = mass_bbar;
    m_m_w    = mass_Wp;
    m_m_wbar = mass_Wm;
    m_mal    = mass_al;
    m_ml     = mass_l;
    m_mv     = 0;
    m_mav    = 0;

    m_beamEnergy = m_config->beamEnergy();  // 13000.;
  }

dileptonTtbarRecoUtils::~dileptonTtbarRecoUtils() {}

void dileptonTtbarRecoUtils::fDelete() const {
    delete this;
}


void dileptonTtbarRecoUtils::setTruthInfo(const Top& top, const Top& antiTop, 
                                          const Neutrino& neutrino, const Neutrino& antiNeutrino) {
    m_true_top         = top;
    m_true_topbar      = antiTop;
    m_true_neutrino    = neutrino;
    m_true_neutrinobar = antiNeutrino);

    filldR();
    filldN();

    return;
}


void dileptonTtbarRecoUtils::sortBy(std::string ch) {
    /* Re-sort the top solutions depending on parameters */
    if (m_ttSol.size()<1) return;


    if(ch.compare("dR")==0) {
        for(uint i=0;i<m_ttSol.size()-1;++i) {
            if(m_ttSol[i].dR > m_ttSol[i+1].dR){ 
                swapTopSol(m_ttSol[i],m_ttSol[i+1]);
                i=-1; // force the loop to restart
            }
        }
    }

    else if(ch.compare("dN")==0) {
        for(uint i=0;i<m_ttSol.size()-1;++i){
            if(m_ttSol[i].dN > m_ttSol[i+1].dN){ 
                swapTopSol(m_ttSol[i],m_ttSol[i+1]);
                i=-1; 
            }
        }
    }

    else if(ch.compare("dRN")==0) {
        for(uint i=0;i<m_ttSol.size()-1;++i){
            if(m_ttSol[i].dN*m_ttSol[i].dR > m_ttSol[i+1].dN*m_ttSol[i+1].dR){ 
                swapTopSol(m_ttSol[i],m_ttSol[i+1]);
                i=-1;
            }
        }
    }

    return;
}


void dileptonTtbarRecoUtils::sortTopSol(std::vector<ttbarDilepton>& v) const {
    //std::vector< dileptonTtbarRecoUtils::ttbarDilepton > result;
    for(uint i=0;i<v.size()-1;++i){
        if(v[i].weight < v[i+1].weight){
            this->swapTopSol(v[i],v[i+1]);
            i=-1; // force the loop to restart (to re-check everything)
        }
    }

    //v.swap(result);
    return;
}


void dileptonTtbarRecoUtils::setdR() {
    /* Set the DeltaR values for top and anti-top wrt truth partons */
    for(int i=0;i<m_NSol;++i){
        m_ttSol[i].dR = sqrt( pow(m_ttSol[i].top.DeltaR(m_true_top),2) +
                              pow(m_ttSol[i].topbar.DeltaR(m_true_topbar),2) );
    }

    return;
}


void dileptonTtbarRecoUtils::setdN() {
    /* Set the difference between reco and truth neutrino 3-vector components */
    for(int i=0;i<m_NSol;++i){
        m_ttSol[i].dN = sqrt(pow((m_ttSol[i].neutrino.p4.Px()    - m_true_neutrino.p4.Px()),2) + 
                             pow((m_ttSol[i].neutrino.p4.Py()    - m_true_neutrino.p4.Py()),2) + 
                             pow((m_ttSol[i].neutrino.p4.Pz()    - m_true_neutrino.p4.Pz()),2) + 
                             pow((m_ttSol[i].neutrinobar.p4.Px() - m_true_neutrinobar.p4.Px()),2) + 
                             pow((m_ttSol[i].neutrinobar.p4.Py() - m_true_neutrinobar.p4.Py()),2) + 
                             pow((m_ttSol[i].neutrinobar.p4.Pz() - m_true_neutrinobar.p4.Pz()),2));
    }

    return;
}


int dileptonTtbarRecoUtils::getNsol() const {
    return m_NSol;
}


const std::vector<ttbarDilepton> dileptonTtbarRecoUtils::getTtSol() const {
    return m_ttSol;
}


void dileptonTtbarRecoUtils::execute(const DileptonReco& ttSystem) {
    /* Set the objects and then solve! */
    m_l    = ttSystem.lepton_neg;
    m_al   = ttSystem.lepton_pos;
    m_b    = ttSystem.bJet;
    m_bbar = ttSystem.bbarJet;
    m_px_miss = ttSystem.met.Px();
    m_py_miss = ttSystem.met.Py();

    execute();

    return;
}


void dileptonTtbarRecoUtils::execute() {
    /* Run the algorithm */
    findCoeff(coeffs_);
    quartic_equation(coeffs_[0],coeffs_[1],coeffs_[2],coeffs_[3],coeffs_[4],vect_pxv_);

    // loop over solutions
    for(int i=1;i<=vect_pxv_[0];++i) {

        topRec(vect_pxv_[i]);     // assign values for physics objects

        ttbarDilepton TS_temp;

        TS_temp.top    = m_top;
        TS_temp.topbar = m_topbar;
        TS_temp.Wplus  = m_w;
        TS_temp.Wminus = m_wbar;
        TS_temp.neutrino    = m_neutrino;
        TS_temp.neutrinobar = m_neutrinobar;

        TS_temp.x1 = (m_top.p4.E()+m_topbar.p4.E()+m_top.p4.Pz()+m_topbar.p4.Pz())/(m_beamEnergy);
        TS_temp.x2 = (m_top.p4.E()+m_topbar.p4.E()-m_top.p4.Pz()-m_topbar.p4.Pz())/(m_beamEnergy);
        TS_temp.mtt = m_tt.M();

        TS_temp.weight = 1.0/m_tt.M();  // landau2D(m_neutrino.p4.E(),m_neutrinobar.p4.E());

        m_ttSol.p4.push_back(TS_temp);
    }

    m_NSol = m_ttSol.size();  // should be the same as vect_pxv_[0]?

    if(m_NSol>0)
        sortTopSol(m_ttSol);

    return;
}


void dileptonTtbarRecoUtils::topRec(const double& px_neutrino) {
    /* Build the ttbar system */
    double pxp, pyp, pzp, pup, pvp, pwp;

    d0_ = d00_;
    d1_ = d11_+d10_*px_neutrino;
    d2_ = d22_+d21_*px_neutrino+d20_*px_neutrino*px_neutrino;

    c0_ = c00_;
    c1_ = c11_+c10_*px_neutrino;
    c2_ = c22_+c21_*px_neutrino+c20_*px_neutrino*px_neutrino;

    pup = px_neutrino;
    pvp = (c0_*d2_-c2_*d0_)/(c1_*d0_-c0_*d1_);
    pwp = (-1)*(a1_+a2_*pup+a3_*pvp)/a4_;

    pxp = m_px_miss-pup;   
    pyp = m_py_miss-pvp;
    pzp = (-1)*(b1_+b2_*pxp+b3_*pyp)/b4_;

    // Define the neutrinos
    m_neutrinobar.SetXYZM(pxp, pyp, pzp, m_mav);
    m_neutrino.SetXYZM(pup, pvp, pwp, m_mv);

    // Define the W bosons
    m_w      = m_al.p4 + m_neutrino.p4;
    m_wbar   = m_l.p4  + m_neutrinobar.p4;

    // Define the ttbar quark system
    m_top.lepton   = m_al;
    m_top.neutrino = m_neutrino;
    m_top.jet      = m_b;
    m_top.set_p4_lep();

    m_topbar.lepton   = m_l;
    m_topbar.neutrino = m_neutrinobar;
    m_topbar.jet      = m_bbar;
    m_topbar.set_p4_lep();

    m_tt = m_top.p4 + m_topbar.p4;

    return;
}


void dileptonTtbarRecoUtils::findCoeff(double* const koeficienty) {
    /* Find the coefficients */
    a1_ = ((m_b.p4.P4.E()+m_al.p4.P4.E())*(mm_w*mm_w-m_mal*m_mal-m_mv*m_mv)-m_al.p4.P4.E()*(m_mt*m_mt-m_mb*m_mb-m_mal*m_mal-m_mv*m_mv)+2*m_b.p4.P4.E()*m_al.p4.P4.E()*m_al.p4.P4.E()-2*m_al.p4.P4.E()*(m_al.Vect().Dot(m_b.Vect())))/(2*m_al.p4.P4.E()*(m_b.p4.P4.E()+m_al.p4.P4.E()));
    a2_ = 2*(m_b.p4.P4.E()*m_al.p4.Px()-m_al.p4.P4.E()*m_b.p4.Px())/(2*m_al.p4.P4.E()*(m_b.p4.P4.E()+m_al.p4.P4.E()));
    a3_ = 2*(m_b.p4.P4.E()*m_al.p4.Py()-m_al.p4.P4.E()*m_b.p4.Py())/(2*m_al.p4.P4.E()*(m_b.p4.P4.E()+m_al.p4.P4.E()));
    a4_ = 2*(m_b.p4.P4.E()*m_al.p4.Pz()-m_al.p4.P4.E()*m_b.p4.Pz())/(2*m_al.p4.P4.E()*(m_b.p4.P4.E()+m_al.p4.P4.E()));

    b1_ = ((m_bbar.p4.P4.E()+m_l.p4.P4.E())*(mm_wbar*mm_wbar-m_ml*m_ml-m_mav*m_mav)-m_l.p4.P4.E()*(mtbar_*mtbar_-m_mbbar*m_mbbar-m_ml*m_ml-m_mav*m_mav)+2*m_bbar.p4.P4.E()*m_l.p4.P4.E()*m_l.p4.P4.E()-2*m_l.p4.P4.E()*(m_l.Vect().Dot(m_bbar.Vect())))/(2*m_l.p4.P4.E()*(m_bbar.p4.E()+m_l.p4.p4.E()));
    b2_ = 2*(m_bbar.p4.P4.E()*m_l.p4.Px()-m_l.p4.P4.E()*m_bbar.p4.Px())/(2*m_l.p4.P4.E()*(m_bbar.p4.P4.E()+m_l.p4.P4.E()));
    b3_ = 2*(m_bbar.p4.P4.E()*m_l.p4.Py()-m_l.p4.P4.E()*m_bbar.p4.Py())/(2*m_l.p4.P4.E()*(m_bbar.p4.P4.E()+m_l.p4.P4.E()));
    b4_ = 2*(m_bbar.p4.P4.E()*m_l.p4.Pz()-m_l.p4.P4.E()*m_bbar.p4.Pz())/(2*m_l.p4.P4.E()*(m_bbar.p4.P4.E()+m_l.p4.P4.E()));


    c22_ = (sqr((mm_w*mm_w-m_mal*m_mal-m_mv*m_mv))-4*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Pz()))*sqr(a1_/a4_)-4*(mm_w*mm_w-m_mal*m_mal-m_mv*m_mv)*m_al.p4.Pz()*(a1_/a4_))/sqr(2*(m_b.p4.P4.E()+m_al.p4.P4.E())); 
    c21_ = (4*(mm_w*mm_w-m_mal*m_mal-m_mv*m_mv)*(m_al.p4.Px()-m_al.p4.Pz()*(a2_/a4_))-8*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Pz()))*(a1_*a2_/sqr(a4_))-8*m_al.p4.Px()*m_al.p4.Pz()*(a1_/a4_))/sqr(2*(m_b.p4.P4.E()+m_al.p4.P4.E())); 
    c20_ = (-4*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Px()))-4*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Pz()))*sqr(a2_/a4_)-8*m_al.p4.Px()*m_al.p4.Pz()*(a2_/a4_))/sqr(2*(m_b.p4.P4.E()+m_al.p4.P4.E())); 
    c11_ = (4*(mm_w*mm_w-m_mal*m_mal-m_mv*m_mv)*(m_al.p4.Py()-m_al.p4.Pz()*(a3_/a4_))-8*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Pz()))*(a1_*a3_/sqr(a4_))-8*m_al.p4.Py()*m_al.p4.Pz()*(a1_/a4_))/sqr(2*(m_b.p4.P4.E()+m_al.p4.P4.E())); 
    c10_ = (-8*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Pz()))*(a2_*a3_/sqr(a4_)) + 8*m_al.p4.Px()*m_al.p4.Py() - 8*m_al.p4.Px()*m_al.p4.Pz()*(a3_/a4_) - 8*m_al.p4.Py()*m_al.p4.Pz()*(a2_/a4_))/sqr(2*(m_b.p4.P4.E()+m_al.p4.P4.E()));
    c00_ = (-4*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Py())) -4*(sqr(m_al.p4.P4.E())-sqr(m_al.p4.Pz()))*sqr(a3_/a4_)-8*m_al.p4.Py()*m_al.p4.Pz()*(a3_/a4_))/sqr(2*(m_b.p4.P4.E()+m_al.p4.P4.E()));


    double D22,D21,D20,D11,D10,D00;
    D22 = (sqr((mm_wbar*mm_wbar-m_ml*m_ml-m_mav*m_mav))-4*(sqr(m_l.p4.P4.E())-sqr(m_l.p4.Pz()))*sqr(b1_/b4_)-4*(mm_wbar*mm_wbar-m_ml*m_ml-m_mav*m_mav)*m_l.p4.Pz()*(b1_/b4_))/sqr(2*(m_bbar.p4.P4.E()+m_l.p4.P4.E())); 
    D21 = (4*(mm_wbar*mm_wbar-m_ml*m_ml-m_mav*m_mav)*(m_l.p4.Px()-m_l.p4.Pz()*(b2_/b4_))-8*(sqr(m_l.p4.P4.E())-sqr(m_l.p4.Pz()))*(b1_*b2_/sqr(b4_))-8*m_l.p4.Px()*m_l.p4.Pz()*(b1_/b4_))/sqr(2*(m_bbar.p4.P4.E()+m_l.p4.P4.E())); 
    D20 = (-4*(sqr(m_l.p4.P4.E())-sqr(m_l.p4.Px()))-4*(sqr(m_l.p4.E())-sqr(m_l.p4.Pz()))*sqr(b2_/b4_)-8*m_l.p4.Px()*m_l.p4.Pz()*(b2_/b4_))/sqr(2*(m_bbar.p4.E()+m_l.p4.E())); 
    D11 = (4*(mm_wbar*mm_wbar-m_ml*m_ml-m_mav*m_mav)*(m_l.p4.Py()-m_l.p4.Pz()*(b3_/b4_))-8*(sqr(m_l.p4.E())-sqr(m_l.p4.Pz()))*(b1_*b3_/sqr(b4_))-8*m_l.p4.Py()*m_l.p4.Pz()*(b1_/b4_))/sqr(2*(m_bbar.p4.E()+m_l.p4.E())); 
    D10 = (-8*(sqr(m_l.p4.E())-sqr(m_l.p4.Pz()))*(b2_*b3_/sqr(b4_)) + 8*m_l.p4.Px()*m_l.p4.Py() - 8*m_l.p4.Px()*m_l.p4.Pz()*(b3_/b4_) - 8*m_l.p4.Py()*m_l.p4.Pz()*(b2_/b4_))/sqr(2*(m_bbar.p4.E()+m_l.p4.E()));
    D00  = (-4*(sqr(m_l.p4.E())-sqr(m_l.p4.Py())) -4*(sqr(m_l.p4.E())-sqr(m_l.p4.Pz()))*sqr(b3_/b4_)-8*m_l.p4.Py()*m_l.p4.Pz()*(b3_/b4_))/sqr(2*(m_bbar.p4.E()+m_l.p4.E()));


    d22_ = D22+sqr(m_px_miss)*D20+sqr(m_py_miss)*D00+m_px_miss*m_py_miss*D10+m_px_miss*D21+m_py_miss*D11;
    d21_ = -D21-2*m_px_miss*D20-m_py_miss*D10;
    d20_ = D20;
    d11_ = -D11-2*m_py_miss*D00-m_px_miss*D10;
    d10_ = D10;
    d00_  = D00;

    // Set coefficient values
    koeficienty[4] = sqr(c00_)*sqr(d22_)+c11_*d22_*(c11_*d00_-c00_*d11_)+c00_*c22_*(sqr(d11_)-2*d00_*d22_)+c22_*d00_*(c22_*d00_-c11_*d11_);
    koeficienty[3] = c00_*d21_*(2*c00_*d22_-c11_*d11_)+c00_*d11_*(2*c22_*d10_+c21_*d11_)+c22_*d00_*(2*c21_*d00_-c11_*d10_) - 
                     c00_*d22_*(c11_*d10_+c10_*d11_)-2*c00_*d00_*(c22_*d21_+c21_*d22_)-d00_*d11_*(c11_*c21_+c10_*c22_)+c11_*d00_*(c11_*d21_+2*c10_*d22_);
    koeficienty[2] = sqr(c00_)*(2*d22_*d20_+sqr(d21_))-c00_*d21_*(c11_*d10_+c10_*d11_)+c11_*d20_*(c11_*d00_-c00_*d11_) + 
                     c00_*d10_*(c22_*d10_-c10_*d22_)+c00_*d11_*(2*c21_*d10_+c20_*d11_)+(2*c22_*c20_+sqr(c21_))*sqr(d00_) - 
                     2*c00_*d00_*(c22_*d20_+c21_*d21_+c20_*d22_)+c10_*d00_*(2*c11_*d21_+c10_*d22_)-d00_*d10_*(c11_*c21_+c10_*c22_)-d00_*d11_*(c11_*c20_+c10_*c21_);
    koeficienty[1] = c00_*d21_*(2*c00_*d20_-c10_*d10_)-c00_*d20_*(c11_*d10_+c10_*d11_)+c00_*d10_*(c21_*d10_+2*c20_*d11_) - 
                     2*c00_*d00_*(c21_*d20_+c20_*d21_)+c10_*d00_*(2*c11_*d20_+c10_*d21_)+c20_*d00_*(2*c21_*d00_-c10_*d11_)-d00_*d10_*(c11_*c20_+c10_*c21_);
    koeficienty[0] = sqr(c00_)*sqr(d20_)+c10_*d20_*(c10_*d00_-c00_*d10_)+c20_*d10_*(c00_*d10_-c10_*d00_)+c20_*d00_*(c20_*d00_-2*c00_*d20_);

    return;
}


void dileptonTtbarRecoUtils::quartic_equation(const double& h0, const double& h1, 
                                              const double& h2, const double& h3, 
                                              const double& h4, std::vector<double>& v) const {
    /* Quartic equation for ttbar solution */
    std::vector<double> result;

    if(sign(a4_)==0 || sign(b4_)==0) {
        result.p4.Push_back(0);
        v.swap(result);
    }
    else{
        if(sign(h0)==0) {
            this->cubic_equation(h1,h2,h3,h4,result);
            v.swap(result);
        }
        else{
            if(sign(h4)==0){
                this->cubic_equation(h0,h1,h2,h3,result);
                result[0]=result[0]+1;
                result.p4.Push_back(0);
                v.swap(result);
            }
            else{
                double H1=h1/h0;
                double H2=h2/h0;
                double H3=h3/h0;
                double H4=h4/h0;
                double K1 = H2 -3*sqr(H1)/8;
                double K2 = H3 + H1*sqr(H1)/8-H1*H2/2;
                double K3 = H4-3*sqr(sqr(H1))/256+sqr(H1)*H2/16-H1*H3/4;

                if(sign(K3)==0){
                    this->cubic_equation(1,0,K1,K2,result);
                    for(int i=1;i<=result[0];++i){
                        result[i]=result[i]-H1/4;
                    }
                    result[0]=result[0]+1;
                    result.p4.Push_back(-H1/4);
                    v.swap(result);
                }
                else{
                    std::vector<double> result_t12;
                    std::vector<double> result_t1;
                    std::vector<double> result_t2;

                    result_t1.p4.Push_back(0);
                    result_t2.p4.Push_back(0);

                    this->cubic_equation(1,2*K1,(K1*K1-4*K3),(-1)*K2*K2,result_t12); 

                    for(int i=1;i<=result_t12[0];++i){
                        if(result_t12[i]>=0){
                            result_t1[0]=result_t1[0]+2;
                            result_t1.p4.Push_back(sqrt(result_t12[i]));
                            result_t1.p4.Push_back((-1)*sqrt(result_t12[i]));
                            result_t2[0]=result_t2[0]+2;
                            result_t2.p4.Push_back((K1+result_t12[i]-K2/sqrt(result_t12[i]))/2);
                            result_t2.p4.Push_back((K1+result_t12[i]+K2/sqrt(result_t12[i]))/2);
                        }
                    }  

                    std::vector<double> pre_result1;

                    result.p4.Push_back(0);
                    for(int i=1;i<=result_t1[0];++i){
                        this->quadratic_equation(1,result_t1[i],result_t2[i],pre_result1);
 
                        for(int j=1;j<=pre_result1[0];++j){
                            int flag=1;
                            for(int r=1;r<=result.at(0);++r) {
                                if(fabs(result.at(r)-pre_result1[j])<0.02)flag=0;
                            }
                            if(flag){
                                result.at(0)=result.at(0)+1;
                                result.p4.Push_back(pre_result1[j]);
                            }
                        }
                        pre_result1.clear();  
                    }  
                    for(int k=1;k<=result.at(0);++k){
                       result.at(k)=result.at(k)-H1/4;
                    }
                    v.swap(result);   
                } // end else sign(K3)!=0
            } // end else sign(h4)!=0
        } // end else sign(h0)!=0
    } // end !(sign(a4_)==0 || sign(b4_)==0)

    return;
}



void dileptonTtbarRecoUtils::cubic_equation(const double& a, const double& b, 
                                            const double& c, const double& d, std::vector<double>& v) const {
    /* Cubic equation for ttbar solution */
    std::vector<double> result;
    if(a==0) {
        this->quadratic_equation(b,c,d,result);
        v.swap(result);
    }
    else{
        double s1 = b/a;
        double s2 = c/a;
        double s3 = d/a;

        double q = (s1*s1-3*s2)/9;
        double q3 = q*q*q;
        double r = (2*s1*s1*s1-9*s1*s2+27*s3)/54;
        double r2 = r*r;
        double S = r2-q3;
        if(sign(S)<0){
            result.p4.Push_back(3);
            double F = acos(r/sqrt(q3));
            result.p4.Push_back(-2*sqrt(fabs(q))*cos(F/3)-s1/3);
            result.p4.Push_back(-2*sqrt(fabs(q))*cos((F+2*TMath::Pi())/3)-s1/3);
            result.p4.Push_back(-2*sqrt(fabs(q))*cos((F-2*TMath::Pi())/3)-s1/3);  
            v.swap(result);
        }
        else{
            if(sign(S)==0){
                long double A = r+sqrt(fabs(r2-q3));
                A = A<0 ? pow(fabs(A),(long double)1.0/3) : -pow(fabs(A),(long double)1.0/3);
                long double B = sign(A) == 0 ? 0 : q/A; 
                result.p4.Push_back(2);
                result.p4.Push_back(A+B-s1/3);
                result.p4.Push_back(-0.5*(A+B)-s1/3);  //!!!
                v.swap(result);
            }
            else{
                long double A = r+sqrt(fabs(r2-q3));
                A = A<0 ? pow(fabs(A),(long double)1.0/3) : -pow(fabs(A),(long double)1.0/3);
                long double B = sign(A) == 0 ? 0 : q/A; 
                result.p4.Push_back(1);
                result.p4.Push_back(A+B-s1/3);
                v.swap(result);
            } // end else sign(S)!=0
        } // end else sign(S)<=0
    } // end

    return;
}



void dileptonTtbarRecoUtils::quadratic_equation(const double& a, const double& b, 
                                                const double& c, std::vector<double>& v) const {
    /* Quadratic equation for reconstruction */
    std::vector<double> result;

    if(a==0){
        this->linear_equation(b,c,result);
        v.swap(result);
    }
    else{
        double D = b*b-4*a*c;

        if(this->sign(D)<0){
            result.p4.Push_back(0);
            v.swap(result);
        }
        else {
            if(sign(D)==0){
                result.p4.Push_back(1);
                result.p4.Push_back((-1)*b/(2*a));
                v.swap(result);
            }
            else{
                result.p4.Push_back(2);
                result.p4.Push_back((-b-sqrt(D))/(2*a));
                result.p4.Push_back((-b+sqrt(D))/(2*a));
                v.swap(result);
            } // end else sign(D)!=0
        } // end else sign(D)>=0
    } // end else a!=0

    return;
}



void dileptonTtbarRecoUtils::linear_equation(const double& a, const double& b, std::vector<double>& v) const {
    /* Linear equation for reconstruction */
    std::vector<double> result;
    if(a==0){
        result.p4.Push_back(0);
        v.swap(result);
    }
    else{
        result.p4.Push_back(1);
        result.p4.Push_back((-1)*(b/a));
        v.swap(result);
    }

    return;
}


void dileptonTtbarRecoUtils::angle_rot(const double& alpha, const double& e, 
                                       const Jet& inJet,    Jet& jet_sm) const {
    /* Angle rotation */
    double px_1, py_1, pz_1; // Coordinate system where momentum is along Z-axis

    // Transition matrix detector -> syst1 ...
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;

    Jet jet = inJet;
    if( std::abs(jet.p4.Px())<=e ) jet.p4.SetPx(0);
    if( std::abs(jet.p4.Py())<=e ) jet.p4.SetPy(0);
    if( std::abs(jet.p4.Pz())<=e ) jet.p4.SetPz(0);

    // Rotation in syst 1 ...
    double phi = 2*TMath::Pi()*m_r3->Rndm();
    pz_1 = jet.p4.Vect().Mag()*cos(alpha);
    px_1 = - jet.p4.Vect().Mag()*sin(alpha)*sin(phi);
    py_1 = jet.p4.Vect().Mag()*sin(alpha)*cos(phi);  

    // Transition detector <- syst1 ...
    if ( std::abs(jet.p4.Py())>0 || std::abs(jet.p4.Pz())>0) {
        double d = sqrt(jet.p4.Pz()*jet.p4.Pz() + jet.p4.Py()*jet.p4.Py());
        double p = jet.p4.Vect().Mag();

        x1 = d/p;
        y1 = 0;
        z1 = jet.p4.Px()/p;

        x2 = - jet.p4.Px()*jet.p4.Py()/d/p;
        y2 = jet.p4.Pz()/d;
        z2 = jet.p4.Py()/p;

        x3 = - jet.p4.Px()*jet.p4.Pz()/d/p;
        y3 = - jet.p4.Py()/d;
        z3 = jet.p4.Pz()/p;

        jet_sm.p4.SetPx(x1*px_1+y1*py_1+z1*pz_1);
        jet_sm.p4.SetPy(x2*px_1+y2*py_1+z2*pz_1);
        jet_sm.p4.SetPz(x3*px_1+y3*py_1+z3*pz_1);
        jet_sm.p4.SetE(jet.p4.E());
    }

    // Smearing for jets with no y/z-components
    if (std::abs(jet.p4.Py())<1e-6 && std::abs(jet.p4.Pz())<1e-6) {
        if ( std::abs(jet.p4.Px())<1e-6 )
            jet_sm.p4.SetPxPyPzE(jet.p4.Px(),jet.p4.Py(),jet.p4.Pz(),jet.p4.E());
        else
            jet_sm.SetPxPyPzE(pz_1,px_1,py_1,jet.p4.E()); // not sure why the momenta are switched...
    }

    return;
}


double dileptonTtbarRecoUtils::sqr(const double& x) const {
    /* Square a value */
    return (x*x);
}


void dileptonTtbarRecoUtils::swap(double& realone, double& realtwo) const {
    /* Swap two variables */
    if(realtwo < realone){
        double aux = realtwo;
        realtwo = realone;
        realone = aux;
    }

    return;
}


int dileptonTtbarRecoUtils::sign(const long double& ld) const {
    /* Get the sign of a value */
    if(fabs(ld)<1e-8) return 0;
    return (ld>0) ? 1 : -1;
}


double dileptonTtbarRecoUtils::landau2D(const double& xv, const double& xvbar) const {
    /* Landau Likelihood function for dilepton reconstruction */
    return 30.641*TMath::Landau(xv,57.941,22.344,0)*TMath::Landau(xvbar,57.533,22.232,0);
    // where do these numbers come from? 
    // This isn't used, so I will ignore it
}

// THE END