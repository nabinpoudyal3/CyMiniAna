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
    m_isCWoLaAnalysis      = m_selection.compare("cwola")==0;

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


bool eventSelection::applySelection(const Event &event) {
    /* Apply cuts 
       Example Cut::
          if (n_jets==3 && n_ljets<1)  FAIL
          else :                       PASS & fill cutflows
    */
    bool passSelection(false);

    m_nominal_weight = event.nominal_weight();
    double first_bin(0.5);            // first bin value in cutflow histogram ("INITIAL")

    // FIRST CHECK IF VALID EVENT FROM TREE
    if(!event.isValidRecoEntry())
        return false;             // skip event


    // fill cutflow histograms with initial value (before any cuts)
    fillCutflows(first_bin);


    // no selection applied
    if (m_dummySelection)
        return true;              // event 'passed'  


    // set physics objects
    m_jets  = event.jets();
    m_ljets = event.ljets();
    m_leptons = event.leptons();
//    m_muons = event.muons();
//    m_electrons = event.electrons();
    m_neutrinos = event.neutrinos();
    m_met = event.met();
    m_ht  = event.HT();
    m_st  = event.ST();

    m_triggers = event.triggers();
    m_filters  = event.filters();
    // add more objects as needed


    m_NLjets     = m_ljets.size();
    m_NJets      = m_jets.size();
    m_NLeptons   = m_leptons.size();

    m_NMuons     = 0;   //m_muons.size();
    m_NElectrons = 0;   //m_electrons.size();
    for (const auto x : m_leptons){
        if (x.isMuon) m_NMuons++;
        else m_NElectrons++;
    }

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
    fillCutflows(first_bin+1);



    // -- All-hadronic DNN 
    if (m_allHadDNNSelection)
        passSelection = allHadDNNSelection(first_bin+2);

    // -- All-hadronic analysis
    else if (m_isZeroLeptonAnalysis){
        m_ttbar0L = event.ttbar0L();
        passSelection = zeroLeptonSelection(first_bin+2);
    }

    // -- Single lepton analysis
    else if (m_isOneLeptonAnalysis){
        m_ttbar1L = event.ttbar1L();
        passSelection = oneLeptonSelection(first_bin+2);
    }

    else if (m_isCWoLaAnalysis){
        m_ttbar1L = event.ttbar1L();
        passSelection = (first_bin+2);
    }

    // -- Dilepton analysis
    else if (m_isTwoLeptonAnalysis){
        m_ttbar2L = event.ttbar2L();
        passSelection = twoLeptonSelection(first_bin+2);
    }

    return passSelection;
}


// ******************************************************* //

// Put selections in functions (allow other selections to call them!)
bool eventSelection::allHadDNNSelection(double cutflow_bin){
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


bool eventSelection::cwoalaSelection(double cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: 1 b-tag
    Jet leptop_ak4  = m_ttbar1L.jet;
    Ljet hadtop_ak8 = m_ttbar1L.ljet;

    bool btagged(false);
    if (leptop_ak4.bdisc>0.5803) btagged=true;
    if (hadtop_ak8.subjet0_bdisc>0.5803) btagged=true;
    if (hadtop_ak8.subjet1_bdisc>0.5803) btagged=true;

    if ( !btagged )     // m_ljets.size()<2
        pass = false;
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    return pass;
}



// ******************************************************* //

bool eventSelection::zeroLeptonSelection(double cutflow_bin){
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
        fillCutflows(cutflow_bin+1);
        pass = true;
    }

    // cut2 :: >=2 ljets  (same as allHadDNN for now)
    pass = allHadDNNSelection(cutflow_bin+2);        // increment the cutflow bin

    // b-tagging cuts? others?

    return pass;
}



// ******************************************************* //

bool eventSelection::ejetsSelection(double cutflow_bin, const Lepton& lep){
    /* Check if event passes selection; called from 1-lepton selection
       -- Following CMS AN-2016/174
    */
    bool pass(false);

    // cut5 :: MET > 50 GeV
    if ( m_met.p4.Pt() < 50 )
        return false;
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // cut6 :: DeltaPhi(e,MET)
    float met_triangle = 1.5*m_met.p4.Pt() / 110.;
    if ( std::abs(lep.p4.DeltaPhi(m_met.p4)-1.5) > met_triangle )
        return false;
    else{
        fillCutflows(cutflow_bin+1);
        pass = true;
    }

    // cut7 :: DeltaPhi(leading AK4,MET)
    if ( std::abs(m_jets.at(0).p4.DeltaPhi(m_met.p4)-1.5) > met_triangle )
        return false;
    else{
        fillCutflows(cutflow_bin+2);
        pass = true;
    }

    return pass;
}


bool eventSelection::mujetsSelection(double cutflow_bin){
    /* Check if event passes selection; called from 1-lepton selection
       -- Following CMS AN-2016/174
    */
    bool pass(false);

    // cut5 :: MET > 35 GeV
    if ( m_met.p4.Pt() < 35 )
        return false;
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    return pass;
}


bool eventSelection::oneLeptonSelection(double cutflow_bin){
    /* Check if event passes selection 
       -- Following CMS AN-2016/174
       -- add b-tagging, others?
       -- separate boosted (AK8) & resolved (4 AK4)?
    */
    bool pass(false);


    // selection based on lepton -- e+jets or mu+jets
    // only do selection if the user requested a specific
    // lepton flavor
    // e.g., if user selected "ejets", don't do mu+jets selection!
    bool ljets  = m_selection.compare("ljets")==0;    // general "lepton+jets" selection
    bool ejets  = m_selection.compare("ejets")==0;
    bool mujets = m_selection.compare("mujets")==0;


    // cut0 :: One lepton
    bool nLeptons(false);
    if (ejets)       nLeptons = (m_NElectrons==1 && m_NMuons==0);
    else if (mujets) nLeptons = (m_NMuons==1 && m_NElectrons==0);
    else if (ljets)  nLeptons = (m_NLeptons==1);

    if ( !nLeptons )
        return false;          // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    Lepton lep = m_leptons.at(0);

    // cut1 :: triggers -- ejets is lepton==electron else mujets
    unsigned int passTrig(0);
    std::vector<std::string> oneLeptonTriggers = (lep.isElectron) ? m_ejetsTriggers : m_mujetsTriggers;
    for (const auto& trig : oneLeptonTriggers){
        if (m_triggers.at(trig)) passTrig++;
    }

    if (passTrig<1)
        return false;
    else{
        fillCutflows(cutflow_bin+1);
        pass = true;
    }


    // cut2 :: >=1 ljets -- Assuming boosted final state
    if ( m_NLjets < 1 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin+2);
        pass = true;
    }

    // cut3 :: >=2 jets (should have 1 AK4 near lepton & 1 AK4 inside the AK8)
    if ( m_NJets < 2 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin+3);
        pass = true;
    }


    // cut4 :: DeltaR(AK4,lepton)
    //         >=1 AK4 jet in the same hemisphere as the electron, 0.3 < R(l,jet) < pi/2
    Jet leptop_ak4 = m_ttbar1L.jet;

    if (!leptop_ak4.isGood)
        return false;
    else{
        fillCutflows(cutflow_bin+4);
        pass = true;
    }

    // cut5 :: DeltaR(AK8,lepton)
    //         >=1 AK8 jet in the opposite hemisphere from the electron, R(l,jet) > pi/2
    //         mark any AK8 jets that don't meet this requirement as "isGood=false"
    Ljet hadtop_ak8 = m_ttbar1L.ljet;
    if (!hadtop_ak8.isGood)
        return false;
    else{
        fillCutflows(cutflow_bin+5);
        pass = true;
    }

    if (ljets){
        // Do the selection based on which lepton flavor we have in the event
        if (m_NElectrons==1) pass = ejetsSelection(cutflow_bin+6,lep);
        else pass = mujetsSelection(cutflow_bin+6);
    }
    else{
        // Only do the el (mu) selection if the user request e+jets (mu+jets)
        pass = (ejets) ? ejetsSelection(cutflow_bin+6,lep) : mujetsSelection(cutflow_bin+6);
    }

    return pass;
}



// ******************************************************* //

bool eventSelection::twoLeptonSelection(double cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: triggers
    bool passTrig(true);
    for (const auto& trig : m_dileptonTriggers)
        if (!m_triggers.at(trig)) passTrig = false;

    if (!passTrig)
        return false;
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // cut1 :: Two leptons
    if ( m_NLeptons != 2 )  // m_NElectrons+m_NMuons
        return false;       // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin+1);
        pass = true;
    }

    // cut2 :: >=2 jets (2 b-jets)
    if ( m_NJets < 2 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin+2);
        pass = true;
    }

    // dilepton/b-tagging/MET/DeltaPhi(MET,jets)/AK4 cuts? others?

    return pass;
}



// -- Helper functions

void eventSelection::fillCutflows(double cutflow_bin){
    /* Fill cutflow histograms with weight at specific bin */
    m_cutflow->Fill(cutflow_bin,m_nominal_weight);  // fill cutflow
    m_cutflow_unw->Fill(cutflow_bin);
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
