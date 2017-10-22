#ifndef HISTOGRAMMERTRUTH_H_
#define HISTOGRAMMERTRUTH_H_

#include "cms-ttbarAC/CyMiniAna/interface/histogrammer.h"

class histogrammerTruth : histogrammer {

  public:

    histogrammerTruth( configuration &cmaConfig );

    void initialize( TFile& outputFile );
    void bookHists();
    using histogrammer::fill;
    using histogrammer::overUnderFlow;
    void fill( Event &event );
};

#endif
