#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "TROOT.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <sstream>

#include "Analysis/CyMiniAna/interface/tools.h"


class configuration {
  public:
    // Default - so root can load based on a name;
    configuration( const std::string &configFile );
    //configuration( const configuration& );
    configuration& operator=( const configuration& rhs );

    // Default - so we can clean up;
    virtual ~configuration();

    // Run once at the start of the job;
    virtual void initialize();
    std::string getConfigOption( std::string item );

    // Print configuration
    virtual void print();

    // Type of File(s) being processed
    virtual bool isMC();              // must call "checkFileType(file)" or "isMC(file)" first!
    virtual bool isMC( TFile& file );
    bool isGridFile();

    // Type of analysis (all-hadronic, semi-leptonic, or di-leptonic)
    virtual bool isZeroLeptonAnalysis();
    virtual bool isOneLeptonAnalysis();
    virtual bool isTwoleptonAnalysis();

    // object declarations
    virtual bool useJets();
    virtual bool useNeutrinos();
    virtual bool useLeptons();
    virtual bool useLargeRJets();
    virtual bool useRCJets();
    virtual bool useTruth();
    virtual bool useFlags();
    virtual bool useTtbar();

    std::string jet_btagWkpt();
    std::vector<std::string> btagWkpts();
    float cMVAv2L() {return m_cMVAv2L;}
    float cMVAv2M() {return m_cMVAv2M;}
    float cMVAv2T() {return m_cMVAv2T;}

    // functions about the TTree
    virtual bool isNominalTree();
    virtual bool isNominalTree( const std::string &tree_name );
    std::vector<std::string> treeNames();
    void setTreename(std::string treeName);
    std::string treename();

    // functions about the file
    virtual void checkFileType( TFile& file );
    std::vector<std::string> filesToProcess();
    void setFilename(std::string fileName);
    std::string filename();

    // return some values from config file
    std::string verboseLevel();
    std::string selection();
    std::vector<std::string> qcdSelections();
    std::string cutsfile();
    std::string outputFilePath();
    std::string customFileEnding();
    std::string configFileName();
    std::string getAbsolutePath();
    int nEventsToProcess();
    unsigned long long firstEvent();
    bool makeNewFile();
    bool makeHistograms();
    bool makeEfficiencies();

    // information for event weights
    std::string metadataFile();
    double XSectionMap ( std::string mcChannelNumber);
    double KFactorMap ( std::string mcChannelNumber );
    double sumWeightsMap ( std::string mcChannelNumber );
    virtual double LUMI();

    // weight systematics
    bool calcWeightSystematics();
    std::map<std::string,unsigned int> mapOfWeightVectorSystematics();
    std::vector<std::string> listOfWeightSystematics();
    std::string listOfWeightSystematicsFile();
    std::string listOfWeightVectorSystematicsFile();

    // DNN & HME
    std::string dnnFile();
    bool getDNN();
    double minDNN();
    double maxDNN();
    std::string dnnKey();   // key for lwtnn
    bool getHME();

    // Reco/Truth event loops
    bool doRecoEventLoop();
    bool doTruthEventLoop();
    bool matchTruthToReco();
    void setMatchTruthToReco(bool truthToReco);

    // misc. for dilepton ttbar
    bool kinematicReco();
    unsigned int NJetSmear();
    unsigned int NMassPoints();
    unsigned int massMin();
    unsigned int massMax();

    float beamEnergy() {return m_beamEnergy;}           // 13000.;
    double topQuarkMass() {return m_topQuarkMass;}      // 172.5
    double bQuarkMass() {return m_bQuarkMass;}          // 4.18
    double WMass() {return m_WMass;}                    // 80.2

    /// All analysis eras as needed
    enum Era{run2_13tev_25ns,     run2_13tev_2015_25ns, run2_13tev_2016_25ns, 
             run2_13tev_25ns_74X, undefined};
    Era convert(const std::string& era);       /// Convert an era from string to enum
    std::string convert(const Era& era);   /// Convert an era from enum to string
    double energyInTev(const Era era) {return 13.;}     /// Return energy for given era in TeV

  protected:

    void check_btag_WP(const std::string &wkpt);

    std::map<std::string,std::string> m_map_config;
    const std::string m_configFile;

    // type of file(s)
    bool m_isMC;
    bool m_isGridFile;

    // type of analysis
    bool m_isZeroLeptonAnalysis;
    bool m_isOneLeptonAnalysis;
    bool m_isTwoleptonAnalysis;

    // object declarations
    bool m_useTruth;
    bool m_useJets;
    bool m_useLeptons;
    bool m_useLargeRJets;
    bool m_useRCJets;
    bool m_useNeutrinos;
    bool m_useFlags;
    bool m_useTtbar;

    // luminosity
    double m_LUMI      = 36074.56; // 2015+2016 luminosity
    double m_LUMI_2015 = 3212.96;
    double m_LUMI_2016 = 32861.6; // OflLumi-13TeV-008

    // return some values from config file
    std::string m_input_selection;
    std::string m_selection;
    std::string m_cutsfile;
    std::string m_treename;
    std::string m_filename;
    std::string m_verboseLevel;
    int m_nEventsToProcess;
    unsigned long long m_firstEvent;
    std::string m_outputFilePath;
    std::string m_customFileEnding;
    bool m_makeNewFile;
    bool m_makeHistograms;
    bool m_makeEfficiencies;
    std::string m_sumWeightsFiles;
    std::string m_cma_absPath;
    std::string m_metadataFile;
    bool m_getDNN;
    std::string m_dnnFile;
    std::string m_dnnKey;
    bool m_getHME;
    bool m_doRecoEventLoop;
    bool m_doTruthEventLoop;
    bool m_matchTruthToReco;

    std::string m_jet_btag_wkpt;   // "L","M","T"
    std::string m_tjet_btag_wkpt;
    std::string m_toptag_wkpt;

    // b-tagging (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco)
    // isBTagged = (jet.cMVAv2 > wkpt)
    std::vector<std::string> m_btag_WPs = {"L","M","T"};
    float m_cMVAv2L=-0.5884;
    float m_cMVAv2M=0.4432;
    float m_cMVAv2T=0.9432;

    std::vector<std::string> m_filesToProcess;
    std::vector<std::string> m_treeNames;

    bool m_calcWeightSystematics;
    std::map<std::string,unsigned int> m_mapOfWeightVectorSystematics;
    std::vector<std::string> m_listOfWeightSystematics;
    std::string m_listOfWeightSystematicsFile;
    std::string m_listOfWeightVectorSystematicsFile;

    std::map<std::string, float> m_XSection; // map file to XSection
    std::map<std::string, float> m_KFactor;  // map file to KFactor
    std::map<std::string, float> m_AMI;      // map file to sum of weights
    std::map<std::string, unsigned int> m_NEvents;   // map file to total number of events

    std::vector<std::string> m_qcdSelections = {"0b0t","0b1t","0b2t","1b0t",
                                                "1b1t","1b2t","2b0t","2b1t", "2b2t"};

    double m_minDNN  = 0.0;   // min. value in the DNN discriminant
    double m_maxDNN  = 1.0;   // max. value in the DNN discriminant

    // -- Top Mass Variables -- //
    const double m_electronMass = 0.000511;
    const double m_muonMass     = 0.105658;
    const double m_bQuarkMass   = 4.8;
    const double m_WMass        = 80.4;
    const double m_topQuarkMass = 172.5;
    const float m_beamEnergy    = 13000.;
    const int SENTINEL    = -1000;
    const int NCHAN       = 4;
    const double m_sqrt_s = 13000;      // center-of-mass energy

    bool m_kinematicReco;
    unsigned int m_NJetSmear;    // 500
    unsigned int m_NMassPoints;  // 500
    unsigned int m_massMin;      // 100
    unsigned int m_massMax;      // 300

    std::map<std::string,std::string> m_defaultConfigs = {
             {"isZeroLeptonAnalysis",  "false"},
             {"isOneLeptonAnalysis",   "false"},
             {"isTwoLeptonAnalysis",   "false"},
             {"useJets",               "false"},
             {"useLeptons",            "false"},
             {"useLargeRJets",         "false"},
             {"useRCJets",             "false"},
             {"useNeutrinos",          "false"},
             {"useTruth",              "false"},
             {"useFlags",              "false"},
             {"jet_btag_wkpt",         "M"},
             {"makeNewFile",           "false"},
             {"makeHistograms",        "false"},
             {"makeEfficiencies",      "false"},
             {"NEvents",               "-1"},
             {"firstEvent",            "0"},
             {"input_selection",       "grid"},
             {"selection",             "example"},
             {"output_path",           "./"},
             {"customFileEnding",      ""},
             {"calcWeightSystematics", "false"},
             {"weightSystematicsFile",       "config/weightSystematics.txt"},
             {"weightVectorSystematicsFile", "config/weightVectorSystematics.txt"},
             {"cutsfile",              "examples/config/cuts_example.txt"},
             {"inputfile",             "examples/config/miniSL_ALLfiles.txt"},
             {"treenames",             "examples/config/treenames_nominal"},
             {"sumWeightsFiles",       "examples/config/miniSL_ALLMCFiles.txt"},
             {"verboseLevel",          "INFO"},
             {"dnnFile",               "config/keras_ttbar_DNN.json"},
             {"dnnKey",                "dnn"},
             {"getDNN",                "false"},
             {"getHME",                "false"},
             {"doRecoEventLoop",       "true"},
             {"doTruthEventLoop",      "false"},
             {"kinematicReco",        "true"},
             {"NJetSmear",             "500"},
             {"NMassPoints",           "1"},
             {"massMin",               "172"},
             {"massMax",               "301"} };
};

#endif
