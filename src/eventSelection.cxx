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
  m_allHadDNNSelection(false){
    m_cuts.resize(0);
    m_cutflowNames.clear();

    m_selection = m_config->selection();
    m_cutsfile  = m_config->cutsfile();
  }

eventSelection::~eventSelection() {}


void eventSelection::initialize() {
    /* Build the cuts using the cut file from configuration */
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

    // Analysis cuts
    m_isZeroLeptonAnalysis = m_selection.compare("allhad")==0;        // all-had
    m_isOneLeptonAnalysis  = m_selection.compare("singlelepton")==0;  // single lepton
    m_isTwoLeptonAnalysis  = m_selection.compare("dilepton")==0;      // dilepton

    return;
}


void eventSelection::setCutflowHistograms(TH1D &cutflow, TH1D &cutflow_unweighted) {
    /* Set the cutflow histograms to use in the framework -- 
       don't need to pass them for every event! 

       Two cutflows:  
         "cutflow"            event weights
         "cutflow_unweighted" no event weights -> raw event numbers
    */
    m_cutflow     = cutflow;
    m_cutflow_unw = cutflow_unweighted;
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
    double cutflow_bin(1.5);          // bin value for subsequent cuts (set below)

    // fill cutflow histograms with initial value (before any cuts)
    m_cutflow.Fill(first_bin,m_nominal_weight); // event weights
    m_cutflow_unw.Fill( first_bin );   // raw event numbers

    // FIRST CHECK IF VALID EVENT FROM TREE
    if(!event.isValidRecoEntry())
        return false;             // skip event


    // no selection applied
    if (m_dummySelection)
        return true;              // event 'passed'  


    // set physics objects
    m_jets  = event.jets();
    m_ljets = event.ljets();
    m_muons = event.muons();
    m_electrons = event.electrons();
    m_neutrinos = event.neutrinos();
    m_met = event.met();
    m_ht  = event.HT();
    m_st  = event.ST();
    // add more objects as needed


    // selection for all-hadronic DNN 
    if (m_allHadDNNSelection){
        passSelection = allHadDNNSelection(cutflow_bin);
    }

    // selection for all-hadronic analysis
    else if (m_isZeroLeptonAnalysis){
        passSelection = zeroLeptonSelection(cutflow_bin);
    }

    // selection for single lepton analysis
    else if (m_isOneLeptonAnalysis){
        passSelection = oneLeptonSelection(cutflow_bin);
    }

    // selection for dilepton analysis
    else if (m_isTwoLeptonAnalysis){
        passSelection = twoLeptonSelection(cutflow_bin);
    }

    return passSelection;
}


// Put selections in functions (allow other selections to call them!)
bool eventSelection::allHadDNNSelection(double cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: >=2 jets 
    if ( m_ljets.size()<2 )
        pass = false;
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // truth information?

    return pass;
}

bool eventSelection::zeroLeptonSelection(double cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: No leptons
    if ( m_electrons.size()>0 || m_muons.size()>0 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // cut1 :: >=2 ljets  (same as allHadDNN for now)
    pass = allHadDNNSelection(cutflow_bin+1);        // increment the cutflow bin

    // b-tagging cuts? others?

    return pass;
}

bool eventSelection::oneLeptonSelection(double cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: One lepton
    unsigned int nElectrons( m_electrons.size() );
    unsigned int nMuons( m_muons.size() );
    if ( nElectrons+nMuons != 1 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // Assuming boosted final state
    // cut1 :: >=1 ljets
    if ( m_ljets.size() < 1 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin+1);
        pass = true;
    }

    // cut2 :: >=2 jets (should have 1 AK4 near lepton & 1 AK4 inside the AK8)
    if ( m_jets.size() < 2 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin+2);
        pass = true;
    }

    // b-tagging/MET/DeltaPhi(MET,jets)/AK4 cuts? others?
    // boosted (AK8) vs resolved (4 AK4)?

    return pass;
}


bool eventSelection::twoLeptonSelection(double cutflow_bin){
    /* Check if event passes selection */
    bool pass(false);

    // cut0 :: One lepton
    unsigned int nElectrons( m_electrons.size() );
    unsigned int nMuons( m_muons.size() );
    if ( nElectrons+nMuons != 2 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin);
        pass = true;
    }

    // cut1 :: >=2 jets (should have 1 AK4 near lepton & 1 AK4 inside the AK8)
    if ( m_jets.size() < 2 )
        return false;  // exit the function now; no need to test other cuts!
    else{
        fillCutflows(cutflow_bin+1);
        pass = true;
    }

    // dilepton/b-tagging/MET/DeltaPhi(MET,jets)/AK4 cuts? others?

    return pass;
}



// -- Helper functions

void eventSelection::fillCutflows(double cutflow_bin){
    /* Fill cutflow histograms with weight at specific bin */
    m_cutflow.Fill(cutflow_bin,m_nominal_weight);  // fill cutflow
    m_cutflow_unw.Fill(cutflow_bin);
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
