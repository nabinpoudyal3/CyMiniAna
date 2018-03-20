/*
Created:        20 August    2016
Last Updated:   16 October   2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Configuration class
  -- Read config file and use functions
     to return configurations later

*/
#include "Analysis/CyMiniAna/interface/configuration.h"


configuration::configuration(const std::string &configFile) : 
  m_configFile(configFile),
  m_isMC(false),
  m_isGridFile(false),
  m_isZeroLeptonAnalysis(false),
  m_isOneLeptonAnalysis(false),
  m_isTwoLeptonAnalysis(false),
  m_useTruth(false),
  m_useJets(false),
  m_useLeptons(false),
  m_useLargeRJets(false),
  m_useNeutrinos(false),
  m_input_selection("SetMe"),
  m_treename("SetMe"),
  m_filename("SetMe"),
  m_verboseLevel("SetMe"),
  m_nEventsToProcess(0),
  m_firstEvent(0),
  m_outputFilePath("SetMe"),
  m_customDirectory("SetMe"),
  m_makeTTree(false),
  m_makeHistograms(false),
  m_sumWeightsFiles("SetMe"),
  m_cma_absPath("SetMe"),
  m_metadataFile("SetMe"),
  m_useDNN(false),
  m_DNNinference(false),
  m_DNNtraining(false),
  m_dnnFile("SetMe"),
  m_dnnKey("SetMe"),
  m_doRecoEventLoop(false),
  m_matchTruthToReco(true),
  m_kinematicReco(true),
  m_jet_btag_wkpt("SetMe"),
  m_calcWeightSystematics(false),
  m_listOfWeightSystematicsFile("SetMe"),
  m_listOfWeightVectorSystematicsFile("SetMe"){
    m_selections.clear();
    m_cutsfiles.clear();

    m_XSection.clear();
    m_KFactor.clear();
    m_AMI.clear();
    m_map_config.clear();
  }

configuration::~configuration() {}

configuration &configuration::operator=(const configuration &rhs) { return *this; }

void configuration::initialize() {
    /* Initialize the configurations */
    std::vector<std::string> configurations; 
    cma::read_file( m_configFile, configurations ); // read config file into vector

    // fill map with values from configuration file
    for (const auto& config : configurations){
        // split config items by space
        std::istringstream cfg(config);
        std::istream_iterator<std::string> start(cfg), stop;
        std::vector<std::string> tokens(start, stop);

        m_map_config.insert( std::pair<std::string,std::string>(tokens.at(0),tokens.at(1)) );
    }

    // Protection against default settings missing in custom configuration
    // -- map of defaultConfigs defined in header (can't use 'verbose' tools, not defined yet!)
    for (const auto& defaultConfig : m_defaultConfigs){
        if ( m_map_config.find(defaultConfig.first) == m_map_config.end() ){ // item isn't in config file
            std::cout << " WARNING :: CONFIG : Configuration " << defaultConfig.first << " not defined" << std::endl;
            std::cout << " WARNING :: CONFIG : Setting value to default " << defaultConfig.second << std::endl;
            m_map_config[defaultConfig.first] = defaultConfig.second;
        }
    }


    // Set the verbosity level (the amount of output to the console)
    std::map<std::string,unsigned int> verboseMap = cma::verboseMap(); // load mapping of string to integer
    m_verboseLevel = getConfigOption("verboseLevel");
    if (verboseMap.find(m_verboseLevel)==verboseMap.end()){
        m_verboseLevel = "INFO";
        cma::setVerboseLevel(m_verboseLevel);

        cma::WARNING( "CONFIG : Verbose level selected, "+m_verboseLevel+", is not supported " );
        cma::WARNING( "CONFIG : Please select one of the following: " );
        for (const auto& dm : verboseMap)
            cma::WARNING( "CONFIG :          "+dm.first);
        cma::WARNING( "CONFIG : Continuing; setting verbose level to "+m_verboseLevel);
    }
    else{
        cma::setVerboseLevel(m_verboseLevel);
    }


    // Get the absolute path to CyMiniAna for loading
    char* cma_path = getenv("CYMINIANADIR");
    if (cma_path==NULL){
        cma::WARNING("CONFIG : environment variable " );
        cma::WARNING("CONFIG :    'CYMINIANADIR' " );
        cma::WARNING("CONFIG : is not set.  Using PWD to set path." );
        cma::WARNING("CONFIG : This may cause problems submitting batch jobs." );
        cma_path = getenv("PWD");
    }
    m_cma_absPath = cma_path;
    cma::DEBUG("CONFIG : path set to: "+m_cma_absPath );

    // Assign values
    m_nEventsToProcess = std::stoi(getConfigOption("NEvents"));
    m_firstEvent       = std::stoi(getConfigOption("firstEvent"));
    m_input_selection  = getConfigOption("input_selection"); // "grid", "pre", etc.
    cma::split( m_map_config.at("selection"), ',', m_selections );  // different event selections
    cma::split( m_map_config.at("cutsfile"), ',', m_cutsfiles );  // different event selections

    m_isZeroLeptonAnalysis = cma::str2bool( getConfigOption("isZeroLeptonAnalysis") );
    m_isOneLeptonAnalysis  = cma::str2bool( getConfigOption("isOneLeptonAnalysis") );
    m_isTwoLeptonAnalysis  = cma::str2bool( getConfigOption("isTwoLeptonAnalysis") );

    if ( (m_isZeroLeptonAnalysis + m_isOneLeptonAnalysis + m_isTwoLeptonAnalysis) != 1 ){
        cma::ERROR("CONFIG : Must choose only one of 'isZeroLeptonAnalysis', 'isOneLeptonAnalysis', 'isTwoLeptonAnalysis'");
        exit(1);
    }


    // check that b-tag and top-tag WPs are recognized as one of supported values
    check_btag_WP(getConfigOption("jet_btag_wkpt"));

    m_jet_btag_wkpt    = getConfigOption("jet_btag_wkpt");
    m_outputFilePath   = getConfigOption("output_path");
    m_customDirectory  = getConfigOption("customDirectory");
    m_sumWeightsFiles  = getConfigOption("sumWeightsFiles");
    m_useTruth         = cma::str2bool( getConfigOption("useTruth") );
    m_useJets          = cma::str2bool( getConfigOption("useJets") );
    m_useLeptons       = cma::str2bool( getConfigOption("useLeptons") );
    m_useLargeRJets    = cma::str2bool( getConfigOption("useLargeRJets") );
    m_useNeutrinos     = cma::str2bool( getConfigOption("useNeutrinos") );
    m_makeTTree        = cma::str2bool( getConfigOption("makeTTree") );
    m_makeHistograms   = cma::str2bool( getConfigOption("makeHistograms") );
    m_makeEfficiencies = cma::str2bool( getConfigOption("makeEfficiencies") );
    m_dnnFile          = getConfigOption("dnnFile");
    m_dnnKey           = getConfigOption("dnnKey");
    m_useDNN           = cma::str2bool( getConfigOption("useDNN") );
    m_DNNinference     = cma::str2bool( getConfigOption("DNNinference") );
    m_DNNtraining      = cma::str2bool( getConfigOption("DNNtraining") );
    m_doRecoEventLoop  = cma::str2bool( getConfigOption("doRecoEventLoop") );
    m_matchTruthToReco = true;  // not needed in this analysis (so it's not a config option) but here in case we do later
    m_kinematicReco    = cma::str2bool( getConfigOption("kinematicReco") );
    m_metadataFile     = getConfigOption("metadataFile");
    m_calcWeightSystematics             = cma::str2bool( getConfigOption("calcWeightSystematics") );
    m_listOfWeightSystematicsFile       = getConfigOption("weightSystematicsFile");
    m_listOfWeightVectorSystematicsFile = getConfigOption("weightVectorSystematicsFile");

    cma::read_file( getConfigOption("inputfile"), m_filesToProcess );
    cma::read_file( getConfigOption("treenames"), m_treeNames );

    m_isGridFile = (m_input_selection.compare("grid")==0) ? true : false;

    m_mapOfSamples.clear();
    cma::getSampleWeights( m_metadataFile,m_mapOfSamples );

    // systematics that are weights in the nominal tree
    m_listOfWeightSystematics.resize(0);
    cma::read_file(m_listOfWeightSystematicsFile,m_listOfWeightSystematics);

    std::vector<std::string> weightVectorSystematics; // vector btagging SF
    cma::read_file(m_listOfWeightVectorSystematicsFile,weightVectorSystematics);

    for (const auto& weightVectorSystematic : weightVectorSystematics){
        // split config items by space
        std::istringstream cfg(weightVectorSystematic);
        std::istream_iterator<std::string> start(cfg), stop;
        std::vector<std::string> tokens(start, stop);

        m_mapOfWeightVectorSystematics.insert( std::pair<std::string,unsigned int>( tokens.at(0),std::stoi(tokens.at(1)) ) );
    }


    // dilepton ttbar reco (not used right now)
    m_NJetSmear   = std::stoi( getConfigOption("NJetSmear") );
    m_NMassPoints = std::stoi( getConfigOption("NMassPoints") );
    m_massMin     = std::stoi( getConfigOption("massMin") );
    m_massMax     = std::stoi( getConfigOption("massMax") );

    return;
}


void configuration::print(){
    // -- Print the configuration
    std::cout << " ** CyMiniAna ** " << std::endl;
    std::cout << " --------------- " << std::endl;
    std::cout << " CONFIGURATION :: Printing configuration " << std::endl;
    std::cout << " " << std::endl;
    for (const auto& config : m_map_config){
        std::cout << " " << config.first << "\t\t\t" << config.second << std::endl;
    }
    std::cout << " --------------- " << std::endl;

    return;
}


std::string configuration::getConfigOption( std::string item ){
    /* Check that the item exists in the map & return it; otherwise throw exception  */
    std::string value("");

    try{
        value = m_map_config.at(item);
    }
    catch(const std::exception&){
        cma::ERROR("CONFIG : Option "+item+" does not exist in configuration.");
        cma::ERROR("CONFIG : This does not exist in the default configuration either.");
        cma::ERROR("CONFIG : Returing an empty string.");
    }

    return value;
}


void configuration::inspectFile( TFile& file ){
    // -- Check the sum of weights tree DSIDs (to determine Data || MC)
    m_isMC = true; // only MC for now -- need to know how CMS does this!
/*
    TTreeReader sumWeights("sumWeights", &file);
    TTreeReaderValue<int> dsid(sumWeights, "dsid");

    std::vector<int> dsids;   // keep track of dsids (just in case)
    unsigned int mc_dsid(0);  // count number of events that have DSID > 0

    while (sumWeights.Next()){
        dsids.push_back(*dsid);
        if (*dsid>0){
            ++mc_dsid;
        } // MC sample (dsid>0)
    }
    m_isMC = (mc_dsid > 0) ? true : false;
*/
    return;
}


void configuration::setTreename(std::string treeName){
    m_treename = treeName;
    return;
}

void configuration::setFilename(std::string fileName){
    m_filename = fileName;
    return;
}

bool configuration::isNominalTree(){
    return isNominalTree( m_treename );
}

bool configuration::isNominalTree( const std::string &tree_name ){
    /* Check if tree is a nominal one */
    bool isNominal(false);
    if (tree_name.compare("nominal")==0 || tree_name.compare("Nominal")==0)
        isNominal = true;
    else
        isNominal = false;

    return isNominal;
}

bool configuration::isMC( TFile& file ){
    /* Check the sum of weights tree DSIDs (to determine Data || MC) */
    inspectFile( file );
    return m_isMC;
}

void configuration::check_btag_WP(const std::string &wkpt){
    /* Check the b-tagging working point */
    if(! std::any_of(m_btag_WPs.begin(), m_btag_WPs.end(), [&](const std::string& s){return (s.compare(wkpt) == 0);} ) ) {
        cma::ERROR("CONFIG : Unknown b-tagging WP: "+wkpt+". Aborting!");
        cma::ERROR("CONFIG : Available b-tagging WPs: "+cma::vectorToStr(m_btag_WPs));
        exit(EXIT_FAILURE);
    }

    return;
}


configuration::Era configuration::convert(const std::string& era) {
    /* Convert string to era enum */
    if(era == "run2_13tev_25ns") return configuration::run2_13tev_25ns;
    else if(era == "run2_13tev_2015_25ns") return configuration::run2_13tev_2015_25ns;
    else if(era == "run2_13tev_2016_25ns") return configuration::run2_13tev_2016_25ns;
    else if(era == "run2_13tev_25ns_74X") return configuration::run2_13tev_25ns_74X;
    else {
        cma::ERROR("CONFIGURATION : ERA convert conversion is not implemented: "+era);
        exit(97);
    }
}

std::string configuration::convert(const Era& era) {
    /* Convert era to string */
    if(era == run2_13tev_25ns) return "run2_13tev_25ns";
    else if(era == run2_13tev_2015_25ns) return "run2_13tev_2015_25ns";
    else if(era == run2_13tev_2016_25ns) return "run2_13tev_2016_25ns";
    else if(era == run2_13tev_25ns_74X) return "run2_13tev_25ns_74X";
    else{
        cma::ERROR("CONFIGURATION : ERA convert conversion is not implemented: "+std::to_string(era));
        exit(97);
    }
}

void configuration::setMatchTruthToReco(bool truthToReco){
    m_matchTruthToReco = truthToReco;
    return;
}

// THE END
