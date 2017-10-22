#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "TROOT.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <sstream>

#include "cms-ttbarAC/CyMiniAna/interface/tools.h"


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

    virtual bool isMC();              // must call "checkFileType(file)" or "isMC(file)" first!
    virtual bool isMC( TFile& file );
    bool isGridFile();

    // object declarations
    virtual bool useJets();
    virtual bool useNeutrinos();
    virtual bool useLeptons();
    virtual bool useLargeRJets();
    virtual bool useRCJets();
    virtual bool useTruth();

    std::string jet_btagWkpt();
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
    bool makeNewFile();
    bool makeHistograms();
    bool makeEfficiencies();

    // information for event weights
    double XSectionMap ( unsigned int );
    double KFactorMap ( unsigned int );
    double sumWeightsMap ( unsigned int );
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
    bool buildNeutrinos();
    unsigned int NJetSmear();
    unsigned int NMassPoints();
    unsigned int massMin();
    unsigned int massMax();

    float beamEnergy() {return m_beamEnergy;}           // 13000.;
    double topQuarkMass() {return m_topQuarkMass;}      // 172.5
    double bQuarkMass() {return m_bQuarkMass;}          // 4.18
    double WMass() {return m_WMass;}                    // 80.2

  protected:

    void check_btag_WP(const std::string &wkpt);

    std::map<std::string,std::string> m_map_config;
    const std::string m_configFile;

    bool m_isMC;
    bool m_isGridFile;
    bool m_useTruth;

    // object declarations
    bool m_useJets;
    bool m_useLeptons;
    bool m_useLargeRJets;
    bool m_useRCJets;
    bool m_useNeutrinos;

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
    std::string m_outputFilePath;
    std::string m_customFileEnding;
    bool m_makeNewFile;
    bool m_makeHistograms;
    bool m_makeEfficiencies;
    std::string m_sumWeightsFiles;
    std::string m_cma_absPath;
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

    std::map<unsigned int, float> m_XSection; // map DSID to XSection
    std::map<unsigned int, float> m_KFactor;  // map DSID to KFactor
    std::map<unsigned int, float> m_AMI;      // map DSID to sum of weights

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

    bool m_buildNeutrinos;
    unsigned int m_NJetSmear;    // 500
    unsigned int m_NMassPoints;  // 500
    unsigned int m_massMin;      // 100
    unsigned int m_massMax;      // 300

    std::map<std::string,std::string> m_defaultConfigs = {
             {"useJets",               "true"},
             {"useLeptons",            "true"},
             {"useLargeRJets",         "true"},
             {"useRCJets",             "false"},
             {"useNeutrinos",          "true"},
             {"useTruth",              "false"},
             {"jet_btag_wkpt",         "70"},
             {"makeNewFile",           "true"},
             {"makeHistograms",        "true"},
             {"makeEfficiencies",      "true"},
             {"NEvents",               "-1"},
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
             {"buildNeutrinos",        "true"},
             {"NJetSmear",             "500"},
             {"NMassPoints",           "1"},
             {"massMin",               "172"},
             {"massMax",               "301"} };
};

#endif
