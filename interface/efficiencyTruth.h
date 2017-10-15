#ifndef EFFICIENCY_TRUTH_H_
#define EFFICIENCY_TRUTH_H_

#include "diHiggs/CyMiniAna/interface/efficiency.h"

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
