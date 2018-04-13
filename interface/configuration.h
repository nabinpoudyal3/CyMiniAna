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
    virtual bool isMC() {return m_isMC;}              // must call "checkFileType(file)" or "isMC(file)" first!
    virtual bool isMC( TFile& file );
    bool isGridFile(){ return m_isGridFile;}
    bool isTtbar(){ return m_isTtbar;}
    bool isQCD(){ return m_isQCD;}


    // Type of analysis (all-hadronic, semi-leptonic, or di-leptonic)
    virtual bool isZeroLeptonAnalysis() {return m_isZeroLeptonAnalysis;}
    virtual bool isOneLeptonAnalysis() {return m_isOneLeptonAnalysis;}
    virtual bool isTwoLeptonAnalysis() {return m_isTwoLeptonAnalysis;}

    // object declarations
    virtual bool useJets() {return m_useJets;}
    virtual bool useNeutrinos() {return m_useNeutrinos;}
    virtual bool useLeptons() {return m_useLeptons;}
    virtual bool useLargeRJets() {return m_useLargeRJets;}
    virtual bool useTruth() {return m_useTruth;}
    bool kinematicReco() {return m_kinematicReco;}
    bool neutrinoReco() {return m_neutrinoReco;}

    std::string jet_btagWkpt() {return m_jet_btag_wkpt;}
    std::vector<std::string> btagWkpts() {return m_btag_WPs;}
    float cMVAv2L() {return m_cMVAv2L;}
    float cMVAv2M() {return m_cMVAv2M;}
    float cMVAv2T() {return m_cMVAv2T;}
    float CSVv2L() {return m_CSVv2L;}
    float CSVv2M() {return m_CSVv2M;}
    float CSVv2T() {return m_CSVv2T;}


    std::vector<std::string> zeroLeptonTriggers() {return m_zeroLeptonTriggers;}
    std::vector<std::string> ejetsTriggers() {return m_ejetsTriggers;}
    std::vector<std::string> mujetsTriggers() {return m_mujetsTriggers;}
    std::vector<std::string> dileptonTriggers() {return m_dileptonTriggers;}

    // functions about the TTree
    virtual bool isNominalTree();
    virtual bool isNominalTree( const std::string &tree_name );
    std::vector<std::string> treeNames() {return m_treeNames;}
    void setTreename(std::string treeName);
    std::string treename() {return m_treename;}

    // functions about the file
    bool checkPrimaryDataset(const std::vector<std::string>& files);
    void readMetadata(TFile& file, const std::string& metadataTreeName);
    virtual void inspectFile( TFile& file, const std::string& metadataTreeName="" );
    std::vector<std::string> filesToProcess() {return m_filesToProcess;}
    bool recalculateMetadata() {return m_recalculateMetadata;}
    void setFilename(std::string fileName);
    std::string filename(){ return m_filename;}
    std::string primaryDataset(){ return m_primaryDataset;}
    unsigned int NTotalEvents(){ return m_NTotalEvents;}

    // return some values from config file
    std::string verboseLevel() {return m_verboseLevel;}
    std::vector<std::string> selections() {return m_selections;}
    std::vector<std::string> qcdSelections() {return m_qcdSelections;}
    std::vector<std::string> cutsfiles() {return m_cutsfiles;}
    std::string outputFilePath() {return m_outputFilePath;}
    std::string customDirectory() {return m_customDirectory;}
    std::string configFileName() {return m_configFile;}
    std::string getAbsolutePath() {return m_cma_absPath;}
    int nEventsToProcess() {return m_nEventsToProcess;}
    unsigned long long firstEvent() {return m_firstEvent;}
    bool makeTTree() {return m_makeTTree;}
    bool makeHistograms() {return m_makeHistograms;}
    bool makeEfficiencies() {return m_makeEfficiencies;}

    // information for event weights
    std::string metadataFile() {return m_metadataFile;}
    std::map<std::string,Sample> mapOfSamples(){return m_mapOfSamples;}
    Sample sample(){return m_sample;}
    virtual double LUMI() {return m_LUMI;}

    // weight systematics
    bool calcWeightSystematics() {return m_calcWeightSystematics;}
    std::map<std::string,unsigned int> mapOfWeightVectorSystematics() {return m_mapOfWeightVectorSystematics;}
    std::vector<std::string> listOfWeightSystematics() {return m_listOfWeightSystematics;}
    std::string listOfWeightSystematicsFile() {return m_listOfWeightSystematicsFile;}
    std::string listOfWeightVectorSystematicsFile() {return m_listOfWeightVectorSystematicsFile;}

    // DNN
    std::string dnnFile() {return m_dnnFile;}
    bool useDNN() {return m_useDNN;}
    bool DNNinference() {return m_DNNinference;}
    bool DNNtraining() {return m_DNNtraining;}
    double minDNN() {return m_minDNN;}
    double maxDNN() {return m_maxDNN;}
    std::string dnnKey() {return m_dnnKey;}   // key for lwtnn

    // Reco/Truth event loops
    bool doRecoEventLoop() {return m_doRecoEventLoop;}
    bool matchTruthToReco() {return m_matchTruthToReco;}
    void setMatchTruthToReco(bool truthToReco);

    // truth-reco matching
    std::map<std::string,int> mapOfPartonContainment() {return m_containmentMap;}
    std::map<int,std::string> mapOfPartonContainmentRev() {return m_containmentMapRev;}
    std::map<std::string,int> mapOfTargetValues() {return m_targetMap;}

    // dilepton ttbar reco
    unsigned int NJetSmear() {return m_NJetSmear;}
    unsigned int NMassPoints() {return m_NMassPoints;}
    unsigned int massMin() {return m_massMin;}
    unsigned int massMax() {return m_massMax;}

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
    bool m_isQCD;
    bool m_isTtbar;
    bool m_isWjets;
    bool m_isSingleTop;

    // type of analysis
    bool m_isZeroLeptonAnalysis;
    bool m_isOneLeptonAnalysis;
    bool m_isTwoLeptonAnalysis;

    // object declarations
    bool m_useTruth;
    bool m_useJets;
    bool m_useLeptons;
    bool m_useLargeRJets;
    bool m_useNeutrinos;
    bool m_neutrinoReco;

    // luminosity
    double m_LUMI      = 36074.56; // 2015+2016 luminosity
    double m_LUMI_2015 = 3212.96;
    double m_LUMI_2016 = 32861.6; // OflLumi-13TeV-008

    // return some values from config file
    std::string m_input_selection;
    std::vector<std::string> m_selections;
    std::vector<std::string> m_cutsfiles;
    std::string m_treename;
    std::string m_filename;
    std::string m_primaryDataset;
    unsigned int m_NTotalEvents;

    std::string m_verboseLevel;
    int m_nEventsToProcess;
    unsigned long long m_firstEvent;
    std::string m_outputFilePath;
    std::string m_customDirectory;
    bool m_makeTTree;
    bool m_makeHistograms;
    bool m_makeEfficiencies;
    std::string m_cma_absPath;
    std::string m_metadataFile;
    bool m_useDNN;
    bool m_DNNinference;
    bool m_DNNtraining;
    std::string m_dnnFile;
    std::string m_dnnKey;
    bool m_doRecoEventLoop;
    bool m_matchTruthToReco;
    bool m_kinematicReco;

    std::string m_jet_btag_wkpt;   // "L","M","T"
    std::string m_tjet_btag_wkpt;
    std::string m_toptag_wkpt;

    // b-tagging (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco)
    // isBTagged = (jet.cMVAv2 > wkpt)
    std::vector<std::string> m_btag_WPs = {"L","M","T"};
    float m_cMVAv2L=-0.5884;
    float m_cMVAv2M=0.4432;
    float m_cMVAv2T=0.9432;
    float m_CSVv2L=0.5426;
    float m_CSVv2M=0.8484;
    float m_CSVv2T=0.9535;

    std::vector<std::string> m_filters = {"goodVertices",
        "eeBadScFilter",
        "HBHENoiseFilter",
        "HBHENoiseIsoFilter",
        "globalTightHalo2016Filter",
        "EcalDeadCellTriggerPrimitiveFilter"};

    std::vector<std::string> m_zeroLeptonTriggers = {"HLT_PFHT800","HLT_PFHT900","HLT_AK8PFJet450","HLT_PFHT700TrimMass50","HLT_PFJet360TrimMass30"};
    std::vector<std::string> m_ejetsTriggers  = {"HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50","HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
    std::vector<std::string> m_mujetsTriggers = {"HLT_Mu40_Eta2P1_PFJet200_PFJet50","HLT_Mu50","HLT_TkMu50"};
    std::vector<std::string> m_dileptonTriggers = {};

    bool m_recalculateMetadata;
    std::vector<std::string> m_filesToProcess;
    std::vector<std::string> m_treeNames;

    bool m_calcWeightSystematics;
    std::map<std::string,unsigned int> m_mapOfWeightVectorSystematics;
    std::vector<std::string> m_listOfWeightSystematics;
    std::string m_listOfWeightSystematicsFile;
    std::string m_listOfWeightVectorSystematicsFile;

    Sample m_sample;                               // struct of information for current sample
    std::map<std::string,Sample> m_mapOfSamples;   // map of Sample structs
    std::map<std::string, float> m_XSection;       // map file to XSection
    std::map<std::string, float> m_KFactor;        // map file to KFactor
    std::map<std::string, float> m_AMI;            // map file to sum of weights
    std::map<std::string, unsigned int> m_NEvents; // map file to total number of events

    std::vector<std::string> m_qcdSelections = {"0b0t","0b1t","0b2t","1b0t",
                                                "1b1t","1b2t","2b0t","2b1t", "2b2t"};

    double m_minDNN  = 0.0;   // min. value in the DNN discriminant
    double m_maxDNN  = 1.0;   // max. value in the DNN discriminant

    // -- Top Mass Variables (dilepton ttbar reco) -- //
    const double m_electronMass = 0.000511;
    const double m_muonMass     = 0.105658;
    const double m_bQuarkMass   = 4.8;
    const double m_WMass        = 80.4;
    const double m_topQuarkMass = 172.5;
    const float m_beamEnergy    = 13000.;
    const int SENTINEL    = -1000;
    const int NCHAN       = 4;
    const double m_sqrt_s = 13000;      // center-of-mass energy
    unsigned int m_NJetSmear;    // 500
    unsigned int m_NMassPoints;  // 500
    unsigned int m_massMin;      // 100
    unsigned int m_massMax;      // 300

    // Primary dataset names for different samples in analysis
    std::map<std::string,std::string> m_mapOfPrimaryDatasets = {
        {"ttbarGOOD","TT_TuneCUETP8M1_13TeV-powheg-pythia8"},
        {"ttbarGEN","TT_TuneCUETP8M1_13TeV-powheg-pythia8"},
        {"singletop_schan","ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"},
        {"singletop_tchan_top","ST_t-channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV-powhegV2-madspin"},
        {"singletop_tchan_antitop","ST_t-channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV-powhegV2-madspin"},
        {"singletop_tWchan_antitop","ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"},
        {"singletop_tWchan_top","ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"},
        {"wjets1","WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"},
        {"wjets2","WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"},
        {"wjets3","WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"},
        {"wjets4","WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"},
        {"qcd","QCD_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt0080","QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt0120","QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt0170","QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt0300","QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt0470","QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt0600","QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt0800","QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_mu_pt1000","QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_el_pt0080","QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_el_pt0120","QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_el_pt0170","QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_el_pt0300","QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_b_pt080","QCD_Pt_80to170_bcToE_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_b_pt170","QCD_Pt_170to250_bcToE_TuneCUETP8M1_13TeV_pythia8"},
        {"qcd_b_pt250","QCD_Pt_250toInf_bcToE_TuneCUETP8M1_13TeV_pythia8"},
    };

    std::vector<std::string> m_qcdFiles   = {"qcd","qcd_mu_pt0080","qcd_mu_pt0120","qcd_mu_pt0170","qcd_mu_pt0300",
                                             "qcd_mu_pt0470","qcd_mu_pt0600","qcd_mu_pt0800","qcd_mu_pt1000",
                                             "qcd_el_pt0080","qcd_el_pt0120","qcd_el_pt0170","qcd_el_pt0300",
                                             "qcd_b_pt080","qcd_b_pt170","qcd_b_pt250"};
    std::vector<std::string> m_ttbarFiles = {"ttbarGOOD","ttbarGEN"};
    std::vector<std::string> m_wjetsFiles = {"wjets1","wjets2","wjets3","wjets4"};
    std::vector<std::string> m_singleTopFiles = {"singletop_schan","singletop_tWchan_antitop",
                                                 "singletop_tWchan_top","singletop_tchan_antitop",
                                                 "singletop_tchan_top","singletop_tchan_top"};

    // Degrees of 'containment' for parton matching to jets
    std::map<std::string,int> m_containmentMap = {
                    {"NONE",  0},
                    {"BONLY", 1},
                    {"QONLY", 2},
                    {"BQ",    3},   // B+Q
                    {"W",     4},   // Q+Q
                    {"FULL",  5}};  // (BQ)+Q || (W+B) || Q+Q+B
    std::map<int,std::string> m_containmentMapRev = {
                    {0,"NONE"},
                    {1,"BONLY"},
                    {2,"QONLY"},
                    {3,"BQ"},     // B+Q
                    {4,"W"},      // Q+Q
                    {5,"FULL"}};  // (BQ)+Q || (W+B) || Q+Q+B
    // map of target values used in training ML of large-R jets
    std::map<std::string,int> m_targetMap = {
      {"none",0},    // :: NONE = QCD (background)
      {"QB",1},      // :: QB-Q = Signal AK8(QB) + AK4(Q)
      {"W",2},       // :: QQ-B = Signal AK8(W)  + AK4(B)
      {"full",3},    // :: FULL = Signal AK8
      {"other",4} }; // :: Other = placeholder (resolved/Q-only/B-only/etc.)


    std::map<std::string,std::string> m_defaultConfigs = {
             {"isZeroLeptonAnalysis",  "false"},
             {"isOneLeptonAnalysis",   "false"},
             {"isTwoLeptonAnalysis",   "false"},
             {"useJets",               "false"},
             {"useLeptons",            "false"},
             {"useLargeRJets",         "false"},
             {"useNeutrinos",          "false"},
             {"neutrinoReco",          "false"},
             {"useTruth",              "false"},
             {"jet_btag_wkpt",         "M"},
             {"makeTTree",             "false"},
             {"makeHistograms",        "false"},
             {"makeEfficiencies",      "false"},
             {"NEvents",               "-1"},
             {"firstEvent",            "0"},
             {"input_selection",       "grid"},
             {"selection",             "example"},
             {"output_path",           "./"},
             {"customDirectory",      ""},
             {"calcWeightSystematics", "false"},
             {"weightSystematicsFile",       "config/weightSystematics.txt"},
             {"weightVectorSystematicsFile", "config/weightVectorSystematics.txt"},
             {"cutsfile",              "examples/config/cuts_example.txt"},
             {"inputfile",             "examples/config/miniSL_ALLfiles.txt"},
             {"treenames",             "examples/config/treenames_nominal"},
             {"treename",              "tree/eventVars"},
             {"metadataFile",          "config/sampleMetaData.txt"},
             {"verboseLevel",          "INFO"},
             {"dnnFile",               "config/keras_ttbar_DNN.json"},
             {"dnnKey",                "dnn"},
             {"useDNN",                "false"},
             {"DNNinference",          "false"},
             {"DNNtraining",           "false"},
             {"doRecoEventLoop",       "true"},
             {"kinematicReco",         "false"},
             {"NJetSmear",             "500"},
             {"NMassPoints",           "1"},
             {"massMin",               "172"},
             {"massMax",               "301"} };
};

#endif
