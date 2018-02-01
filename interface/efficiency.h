#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <string>
#include <map>
#include <vector>

#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/configuration.h"

class efficiency {
  public:

    // Default
    efficiency(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~efficiency();

    /* initialize efficiencies (1D & 2D) */
    virtual void init_eff(const TH1 &passed, const TH1 &total); // from existing histograms
    virtual void init_eff( const std::string &name, 
                                  const unsigned int nBins, const double x_min, const double x_max );
    virtual void init_eff( const std::string &name, 
                                  const unsigned int nBins, const double *xbins );

    virtual void init_eff( const std::string &name, 
                                  const unsigned int nBinsX, const double x_min, const double x_max,
                                  const unsigned int nBinsY, const double y_min, const double y_max );
    virtual void init_eff( const std::string &name, 
                                  const unsigned int nBinsX, const double *xbins,
                                  const unsigned int nBinsY, const double *ybins );

    /* fill efficiencies */
    virtual void fill( Event &event );
    virtual void fill( const std::string &name, const double &value, const bool &decision, const double &weight );
    virtual void fill( const std::string &name, const double &xvalue, const double &yvalue, const bool &decision, const double &weight );

    /* Book efficiencies */
    virtual void bookEffs( TFile& outputFile );

  protected:

    configuration *m_config;

    std::map<std::string, TEfficiency*> m_map_efficiencies;

    std::vector<std::string> m_names;
};

#endif
