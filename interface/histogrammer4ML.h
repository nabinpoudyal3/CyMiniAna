#ifndef HISTOGRAMMER4ML_H_
#define HISTOGRAMMER4ML_H_

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <string>
#include <map>
#include <vector>

#include "Analysis/CyMiniAna/interface/histogrammer.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/Event.h"

class histogrammer4ML : public histogrammer {
  public:

    // Default - so root can load based on a name;
    histogrammer4ML( configuration& cmaConfig, std::string name="" );

    // Default - so we can clean up;
    virtual ~histogrammer4ML();

    /* fill histograms */
    void fill( const std::map<std::string,double> top, double weight=1.0 );

    /* Book histograms */
    void initialize( TFile& outputFile );
    void bookHists();

  protected:

    configuration *m_config;
    std::string m_name;

    // Target values for system
    std::vector<std::string> m_targets = {"0","1","2"};
};

#endif

