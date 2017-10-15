/*
Created:        --
Last Updated:   29 August    2017

Oliver Majersky
oliver.majersky@cernSPAMNOT.ch

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

For creating histograms of ONLY truth information.
No systematics, other aspects necessary.

Originally developed by Oliver for ATLAS tasks.
Kept here for completion & no need to duplicate work in the future.

*/
#include "diHiggs/CyMiniAna/interface/histogrammerTruth.h"


histogrammerTruth::histogrammerTruth( configuration &cmaConfig ) : 
  histogrammer::histogrammer(cmaConfig, "_truth"){}

void histogrammerTruth::initialize( TFile& outputFile ) {
    /* Initialize truth histograms -- no complicated setup needed */
    bookHists();
}

void histogrammerTruth::bookHists() {
    /* Book histograms */
    // init_hist();
}

void histogrammerTruth::fill( Event &event ) {
    /* fill truth histograms -- don't need multiple functions to handle systematics/weights */
    // calculate value from event object
    // fill();
}
