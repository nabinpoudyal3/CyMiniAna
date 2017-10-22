/*
Created:        --
Last Updated:   26 August    2017

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

  Contact Dan with any questions about using this
  or suggestions for your specific use case
*/
#include "cms-ttbarAC/CyMiniAna/interface/eventSelection.h"



eventSelection::eventSelection(configuration &cmaConfig, const std::string &level) :
  m_config(&cmaConfig),
  m_level(level),
  m_selection("SetMe"),
  m_cutsfile("SetMe"),
  m_numberOfCuts(0),
  m_exampleSelection(false),
  m_example2Selection(false){
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
    m_exampleSelection  = false;
    m_example2Selection = false;

    // Define selections to be exclusive (if/else if), but this could be modified
    // Use either m_selection or m_level to make these decisions
    if (m_selection.compare("example")==0)       m_exampleSelection  = true;
    else if (m_selection.compare("example2")==0) m_example2Selection = true;

    return;
}


bool eventSelection::applySelection(Event &event, TH1D &cutflow, TH1D &cutflow_unweighted) {
    /* Apply cuts 

       Two cutflows:  
         "cutflow"            event weights
         "cutflow_unweighted" no event weights -> raw event numbers

       Example Cut::
          if (n_jets==3 && n_ljets<1)  FAIL
          else :                       PASS & fill cutflows
    */
    float nominal_weight  = event.nominal_weight();
    double first_bin(0.5);            // first bin value in cutflow histogram ("INITIAL")
                                      // easily increment by 1 for each cut (don't need to remeber values)

    // fill cutflow histograms with initial value (before any cuts)
    cutflow.Fill(first_bin,nominal_weight); // event weights
    cutflow_unweighted.Fill( first_bin );   // raw event numbers

    // FIRST CHECK IF VALID EVENT FROM TREE
    // this is particularly important when looping over truth tree,
    // where a match in reco tree is not guaranteed
    if(!event.isValidRecoEntry()) return false;

    // Example #1
    if (m_exampleSelection){
        std::vector<Jet> jets = event.jets();  // access some event information
        // cut0 :: >=1 jets 
        if ( jets.size()<1 )                   // check if the event passes the cut!
            return false;
        else{
            cutflow.Fill(first_bin+1,nominal_weight);  // fill cutflow
            cutflow_unweighted.Fill(first_bin+1);
        }
    }

    return true;
}


void eventSelection::getCutNames(){
    /* Get the cut names (for labeling bins in cutflow histograms) and store in vector */
    m_cutflowNames.clear();
    for (const auto& cut : m_cuts)
        m_cutflowNames.push_back( cut.name );

    return;
}

std::vector<std::string> eventSelection::cutNames(){
    /* Return a vector of the cut names (for labeling bins in cutflow histograms) */
    return m_cutflowNames;
}

unsigned int eventSelection::numberOfCuts(){
    /* Return the number of cuts (number of bins in cutflow histograms) */
    return m_numberOfCuts;
}

// the end
