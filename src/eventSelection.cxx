/*
Created:        26 August 2017
Last Updated:    2 March  2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Event Selection script

-----
m_level        String passed to this class
m_selection    Defined in configuration file

Recommended Use:
  - If you have a 'general' selection, e.g., "QCD", 
    consisting of many regions (a la ABCD method),
    set 'm_selection' to 'qcd' and make different instances
    of 'eventSelection' for each region.
    Then, level will represent each region ("A", "B", "C", "D").
    You can make more descriptive names as well.
  - If you just have one region to define, do not pass an
    argument for 'level' and use 'm_selection' to define the cuts
*/
#include "Analysis/CyMiniAna/interface/eventSelection.h"


eventSelection::eventSelection(configuration &cmaConfig, const std::string &level) :
  m_config(&cmaConfig),
  m_level(level),
  m_selection("SetMe"),
  m_cutsfile("SetMe"),
  m_numberOfCuts(0),
  m_dummySelection(false),
  m_isZeroLeptonAnalysis(false),
  m_isOneLeptonAnalysis(false),
  m_isTwoLeptonAnalysis(false),
  m_isCWoLaAnalysis(false),
  m_allHadDNNSelection(false){
    m_cuts.resize(0);
    m_cutflowNames.clear();
  }

eventSelection::~eventSelection() {}

void eventSelection::initialize(const std::string& selection, const std::string& cutsfile) {
    /* Build the cuts using the cut file from configuration */
    m_selection = selection;
    m_cutsfile  = cutsfile;

    m_zeroLeptonTriggers = m_config->zeroLeptonTriggers();
    m_ejetsTriggers      = m_config->ejetsTriggers();
    m_mujetsTriggers     = m_config->mujetsTriggers();
    m_dileptonTriggers   = m_config->dileptonTriggers();

    initialize( m_cutsfile );

    return;
}

void eventSelection::initialize(const std::string &cutsfile) {
    /* Load cut values using specific name for cutsfile */
    std::ifstream file = cma::open_file(cutsfile);

    // Read one line at a time into the vector of Cut structs:
    // this only stores information, but can be expanded
    m_cuts.clear();
    std::string line;
    if (file.is_open()){
        while(std::getline(file, line)){
            std::stringstream  lineStream(line);
            Cut tmp_cut;
            // read line
            lineStream >> tmp_cut.name >> tmp_cut.comparison >> tmp_cut.value;
            m_cuts.push_back(tmp_cut);
        }
        file.close();
    } // end reading cuts file

    // Get the number of cuts (for cutflow histogram binning)
    m_numberOfCuts = m_cuts.size();

    // Get the names of cuts (for cutflow histogram bin labeling)
    m_cutflowNames.clear();
    getCutNames();

    // Identify the selection this instance will apply
    identifySelection();

    return;
}


void eventSelection::identifySelection(){
    /* Set the booleans for applying the selection below */
    m_dummySelection       = m_selection.compare("none")==0;          // no selection

    // DNN cuts -- for preparing DNN ntuples/histograms/etc.
    m_allHadDNNSelection   = m_selection.compare("allHadDNN")==0;     // setup for all-had DNN

    // Analysis cuts :: all-had, l+jets, dilepton
    m_isZeroLeptonAnalysis = m_selection.compare("allhad")==0;
    m_isOneLeptonAnalysis  = (m_selection.compare("ejets")==0 || m_selection.compare("mujets")==0 || m_selection.compare("ljets")==0);
    m_isTwoLeptonAnalysis  = m_selection.compare("dilepton")==0;

    m_isCWoLaAnalysis = m_selection.compare("cwola")==0 || m_selection.compare("cwolaejets")==0 || m_selection.compare("cwolamujets")==0;
    m_isAFBAnalysis   = m_selection.compare("afb")==0;

    return;
}


void eventSelection::setCutflowHistograms(TFile& outputFile){
    /* Set the cutflow histograms to use in the framework -- 
       can modify this function to generate histograms with different names
       e.g., based on the name of the TTree 

       Two cutflows:  
         "cutflow"            event weights
         "cutflow_unweighted" no event weights -> raw event numbers
    */
    outputFile.cd();

    m_cutflow     = new TH1D( (m_selection+"_cutflow").c_str(),(m_selection+"_cutflow").c_str(),m_numberOfCuts+1,0,m_numberOfCuts+1);
    m_cutflow_unw = new TH1D( (m_selection+"_cutflow_unweighted").c_str(),(m_selection+"_cutflow_unweighted").c_str(),m_numberOfCuts+1,0,m_numberOfCuts+1);

    m_cutflow->GetXaxis()->SetBinLabel(1,"INITIAL");
    m_cutflow_unw->GetXaxis()->SetBinLabel(1,"INITIAL");

    for (unsigned int c=1;c<=m_numberOfCuts;++c){
        m_cutflow->GetXaxis()->SetBinLabel(c+1,m_cutflowNames.at(c-1).c_str());
        m_cutflow_unw->GetXaxis()->SetBinLabel(c+1,m_cutflowNames.at(c-1).c_str());
    }

    return;
}


void eventSelection::finalize() {
    /* Clean-up */
    return;
}


void eventSelection::setObjects(const Event& event) {
    /* Set member variables of the physics objects for all selections */
    m_nominal_weight = event.nominal_weight();

    // FIRST CHECK IF VALID EVENT FROM TREE
    m_valid = event.isValidRecoEntry();

    // set physics objects
    m_jets  = event.jets();
    m_ljets = event.ljets();
    m_leptons = event.leptons();
    //m_muons = event.muons();
    //m_electrons = event.electrons();
    m_neutrinos = event.neutrinos();
    m_met = event.met();
    m_ht  = event.HT();
    m_st  = event.ST();

    m_triggers = event.triggers();
    m_filters  = event.filters();
    // add more objects as needed

    m_Nbtags     = event.btag_jets().size();
    m_NLjets     = m_ljets.size();
    m_NJets      = m_jets.size();
    m_NLeptons   = m_leptons.size();

    m_NMuons     = 0;   //m_muons.size();
    m_NElectrons = 0;   //m_electrons.size();
    for (const auto x : m_leptons){
        if (x.isMuon) m_NMuons++;
        else m_NElectrons++;
    }

    // ttbar system(s)
    m_ttbar0L = event.ttbar0L();
    m_ttbar1L = event.ttbar1L();
    m_ttbar2L = event.ttbar2L();

    return;
}


bool eventSelection::applySelection() {
    /* Apply cuts 
       Example Cut::
          if (n_jets==3 && n_ljets<1)  FAIL
          else :                       PASS & fill cutflows
    */
    bool passSelection(false);

    float cf_bin(0.5);            // bin value in cutflow histogram ("INITIAL")

    // FIRST CHECK IF VALID EVENT FROM TREE
    if(!m_valid)
        return false;             // skip event


    // fill cutflow histograms with initial value (before any cuts)
    fillCutflows(cf_bin);      // fillCutflows() iterates 'cf_bin'


    // no selection applied
    if (m_dummySelection)
        return true;              // event 'passed'  


    // Perform selections
    // -- use if/else if statements to maintain orthogonality

    // -- Filter (only necessary for data)
    bool passFilter(true);
    if (!m_config->isMC()){
        for (const auto& x : m_filters){
            if (!x.second){
                passFilter = false;
                break;
            }
        }
    }
    if (!passFilter) return false;
    fillCutflows(cf_bin);


    // -- All-hadronic DNN 
    if (m_allHadDNNSelection)
        passSelection = allHadDNNSelection(cf_bin);

    // -- All-hadronic analysis
    else if (m_isZeroLeptonAnalysis){
        passSelection = zeroLeptonSelection(cf_bin);
    }

    // -- Single lepton analyses
    else if (m_isOneLeptonAnalysis){
        // Standard 1-lepton selection
        passSelection = (m_selection.compare("ejets")==0) ? ejetsSelection(cf_bin) : mujetsSelection(cf_bin);
    }
    else if (m_isCWoLaAnalysis){
        // inputs to CWoLa
        passSelection = cwolaSelection(cf_bin);
    }
    else if (m_isAFBAnalysis){
        // mimc the AFB analysis as a cross-check
        passSelection = afbSelection(cf_bin);
    }

    // -- Dilepton analysis
    else if (m_isTwoLeptonAnalysis){
        passSelection = twoLeptonSelection(cf_bin);
    }

    return passSelection;
}


// ******************************************************* //

// Put selections in functions (allow other selections to call them!)
bool eventSelection::allHadDNNSelection(float& cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: >=2 jets pT>400 GeV -- trim down the files
    if ( m_NLjets<2 )     // m_ljets.size()<2
        pass = false;
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // truth information?

    return pass;
}


// ******************************************************* //
bool eventSelection::zeroLeptonSelection(float& cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: triggers
    unsigned int passTrig(0);
    for (const auto& trig : m_zeroLeptonTriggers){
        if (m_triggers.at(trig)) passTrig++;
    }

    if (passTrig<1)
        return false;
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }


    // cut1 :: No leptons
    if ( m_NLeptons>0 )    // m_NElectrons+m_NMuons
        return false;      // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // cut2 :: >=2 ljets  (same as allHadDNN for now)
    pass = allHadDNNSelection(cutflow_bin);        // increment the cutflow bin

    // b-tagging cuts? others?

    return pass;
}



// ******************************************************* //
bool eventSelection::afbSelection(float& cutflow_bin){
    /* Forward-backward asymmetry selection (cross-check of framework) */
    /**** PRE SELECTION ****/
    // cut0 :: >=1 One lepton
    if (m_NLeptons==0)
        return false;          // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut1 :: >=1 top tags
    unsigned int n_toptags(0);
    for (const auto& ljet : m_ljets) {
        if (ljet.tau32 < 0.80 && 105<ljet.softDropMass && ljet.softDropMass<220.)
            n_toptags++;
    }

    if ( n_toptags < 1 )
        return false;  // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut2 :: Minimal AK4 kinematics
    if ( m_jets.size() < 1 )
        return false;  // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    /**** REGULAR SELECTION ****/
    // cut3 :: One lepton (kinematics only)
    if (m_NMuons!=1 || m_NElectrons>0)
        return false;          // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut4 :: 2D isolation
    Lepton lep = m_leptons.at(0);
    if ( !lep.iso ) 
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut5 :: triggers -- ejets is lepton==electron else mujets
    unsigned int passTrig(0);
    for (const auto& trig : m_mujetsTriggers){
        if (m_triggers.at(trig)) passTrig++;
    }

    if (passTrig<1)
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut6 :: >=2 jets (should have 1 AK4 near lepton & 1 AK4 inside the AK8)
    if ( m_NJets < 2 )
        return false;  // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut7 :: leading AK4 pT>150; sub-leading AK4 pT>50
    if ( m_jets.at(0).p4.Pt() < 150. || m_jets.at(1).p4.Pt() < 50.)
        return false;  // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut8 :: >= 1 b-tags
    if ( m_Nbtags<1 )
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut9 :: MET > 50 GeV
    if ( m_met.p4.Pt() < 50 )
        return false;
    else
        fillCutflows(cutflow_bin);

    return true;
}


bool eventSelection::cwolaSelection(float& cutflow_bin){
    /* Check if event passes selection */
    // Standard 1-lepton selection
    bool pass_basic(false);
    if (m_selection.find("ejets")!=std::string::npos)
        pass_basic = ejetsSelection(cutflow_bin,"ejets");
    else
        pass_basic = mujetsSelection(cutflow_bin,"mujets");

    if (!pass_basic) 
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut0 :: 1 b-tag
    if ( m_Nbtags<1 )
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut1 :: BEST(top)>0.2
    Ljet hadtop_ak8 = m_ttbar1L.ljet;
    if ( hadtop_ak8.BEST_t<0.1 )
        return false;
    else
        fillCutflows(cutflow_bin);

    return true;
}


bool eventSelection::ejetsSelection(float& cutflow_bin, std::string selection){
    /* Check if event passes selection; called from 1-lepton selection
       -- Following CMS AN-2016/174
    */
    if (selection.size()<1) selection = m_selection;

    // Standard 1-lepton selection
    bool pass = oneLeptonSelection(cutflow_bin,selection);
    if (!pass)
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut5 :: MET > 50 GeV
    if ( m_met.p4.Pt() < 50 )
        return false;
    else
        fillCutflows(cutflow_bin);

    Lepton lep = m_leptons.at(0);
    // cut6 :: DeltaPhi(e,MET)
    float met_triangle = 1.5*m_met.p4.Pt() / 110.;
    if ( std::abs(lep.p4.DeltaPhi(m_met.p4)-1.5) > met_triangle )
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut7 :: DeltaPhi(leading AK4,MET)
    if ( std::abs(m_jets.at(0).p4.DeltaPhi(m_met.p4)-1.5) > met_triangle )
        return false;
    else
        fillCutflows(cutflow_bin);

    return true;
}


bool eventSelection::mujetsSelection(float& cutflow_bin, std::string selection){
    /* Check if event passes selection; called from 1-lepton selection
       -- Following CMS AN-2016/174
    */
    if (selection.size()<1) selection = m_selection;

    // Standard 1-lepton selection
    bool pass = oneLeptonSelection(cutflow_bin,selection);
    if (!pass)
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut5 :: MET > 35 GeV
    if ( m_met.p4.Pt() < 35 )
        return false;
    else
        fillCutflows(cutflow_bin);

    return true;
}


bool eventSelection::oneLeptonSelection(float& cutflow_bin, std::string selection){
    /* Single lepton selection following CMS AN-2016/174 (el+jets or mu+jets)
       -- add b-tagging, others?
       -- separate boosted (AK8) & resolved (4 AK4)?
       Only perform relevant selection
       e.g., if user selected "ejets", don't do mu+jets selection!
    */
    if (selection.size()<1) selection = m_selection;

    bool ejets  = selection.compare("ejets")==0;
    bool mujets = selection.compare("mujets")==0;
    bool ljets  = ((!ejets && !mujets) || selection.compare("ljets")==0);  // general "lepton+jets" selection

    // cut0 :: One lepton
    bool nLeptons(false);
    if (ejets)       nLeptons = (m_NElectrons==1 && m_NMuons==0);
    else if (mujets) nLeptons = (m_NMuons==1 && m_NElectrons==0);
    else if (ljets)  nLeptons = (m_NLeptons==1);

    if ( !nLeptons )
        return false;               // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut1 :: triggers -- ejets is lepton==electron else mujets
    unsigned int passTrig(0);
    std::vector<std::string> oneLeptonTriggers = (m_leptons.at(0).isElectron) ? m_ejetsTriggers : m_mujetsTriggers;
    for (const auto& trig : oneLeptonTriggers){
        if (m_triggers.at(trig)) passTrig++;
    }

    if (passTrig<1)
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut2 :: >=1 ljets -- Assuming boosted final state
    if ( m_NLjets < 1 )
        return false;  // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut3 :: >=2 jets (should have 1 AK4 near lepton & 1 AK4 inside the AK8)
    if ( m_NJets < 2 )
        return false;  // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);


    // cut4 :: DeltaR(AK4,lepton)
    //         >=1 AK4 jet in the same hemisphere as the electron, 0.3 < R(l,jet) < pi/2
    Jet leptop_ak4 = m_ttbar1L.jet;
    if (!leptop_ak4.isGood)
        return false;
    else
        fillCutflows(cutflow_bin);


    // cut5 :: DeltaR(AK8,lepton)
    //         >=1 AK8 jet in the opposite hemisphere from the electron, R(l,jet) > pi/2
    //         mark any AK8 jets that don't meet this requirement as "isGood=false"
    Ljet hadtop_ak8 = m_ttbar1L.ljet;
    if (!hadtop_ak8.isGood)
        return false;
    else
        fillCutflows(cutflow_bin);

    return true;
}



// ******************************************************* //
bool eventSelection::twoLeptonSelection(float& cutflow_bin){
    /* Check if event passes selection */
    // cut0 :: triggers
    bool passTrig(true);
    for (const auto& trig : m_dileptonTriggers)
        if (!m_triggers.at(trig)) passTrig = false;

    if (!passTrig)
        return false;
    else
        fillCutflows(cutflow_bin);

    // cut1 :: Two leptons
    if ( m_NLeptons != 2 )  // m_NElectrons+m_NMuons
        return false;       // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);

    // cut2 :: >=2 jets (2 b-jets)
    if ( m_NJets < 2 )
        return false;  // exit the function now; no need to test other cuts!
    else
        fillCutflows(cutflow_bin);

    // dilepton/b-tagging/MET/DeltaPhi(MET,jets)/AK4 cuts? others?

    return true;
}



// -- Helper functions

void eventSelection::fillCutflows(float& cutflow_bin){
    /* Fill cutflow histograms with weight at specific bin */
    m_cutflow->Fill(cutflow_bin,m_nominal_weight);  // fill cutflow
    m_cutflow_unw->Fill(cutflow_bin);

    cutflow_bin++;                                  // iterate the bin here (don't have to keep track elsewhere)
    return;
}

void eventSelection::getCutNames(){
    /* Get the cut names (for labeling bins in cutflow histograms) and store in vector */
    m_cutflowNames.clear();
    for (const auto& cut : m_cuts)
        m_cutflowNames.push_back( cut.name );

    return;
}

// the end
