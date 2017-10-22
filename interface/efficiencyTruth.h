#ifndef EFFICIENCY_TRUTH_H
#define EFFICIENCY_TRUTH_H

#include "cms-ttbarAC/CyMiniAna/interface/efficiency.h"

class efficiencyTruth : efficiency {
  public:
    
    efficiencyTruth(configuration &cmaConfig);

    /* fill efficiencies */
    using efficiency::fill;
    void fill( Event &event, bool decision );

    /* Book efficiencies */
    void bookEffs( TFile& outputFile );
};

#endif
