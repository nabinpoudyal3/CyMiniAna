// Example for making simple histograms
// :: g++ simpleHistogrammer.cxx -o simpleHistogrammer `root-config --cflags --glibs` -lTreePlayer
// :: ./simpleHistogrammer <filename>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <sstream>
#include <stdio.h>



//----------------------------------------------------------------------
int main(int argc, char* argv[]){

    std::string sample_name = std::string( argv[1] );          // e.g., "TTS_M1100.root"
    std::string tmp_sample_name(sample_name);

    // protection against having the file in a directory somewhere
    if (tmp_sample_name.find("/")!=std::string::npos){
        std::size_t found = tmp_sample_name.find_last_of("/");
        tmp_sample_name   = tmp_sample_name.substr(found+1);
    }
    std::string outfile_name = "simpleHists_"+tmp_sample_name; // e.g., "simpleHists_TTS_M1100.root"

    // outfile of histograms
    TFile f(outfile_name.c_str(),"recreate");

    // create a simple histogram
    auto h_pt0 = TH1F("signal_pt0", ";Leading Small-R Jet p_T;Events", 40,0,2000);


    // -- load root file with data we are analyzing
    auto file = TFile::Open(sample_name.c_str());
    if (!file || file->IsZombie()) {
        std::cout << "   -- File: " << sample_name << std::endl;
        std::cout << "      does not exist or it is a Zombie. " << std::endl;
        std::cout << "      Exiting. " << std::endl;
        return -1;
    }

    // -- Create a TTreeReader for the tree
    TTreeReader myReader("nominal", file); // "nominal" tree is the one we want

    // -- Access branches in the tree
    TTreeReaderValue<std::vector<float>> jet_pt(myReader,   "jet_pt");


    // --------------------------------------------------------------
    // Loop over all entries in the input root file
    // --------------------------------------------------------------
    unsigned int jentry(0.);
    while (myReader.Next()) {

        if (jentry%100000==0)
            std::cout << "Entry: " << jentry << std::endl;

        // check for at least 1 large-R jet in the event; only if you need large-R jets!
        if (jet_pt->size() < 3)
            continue;

        // histogram of pT
        h_pt0.Fill( jet_pt->at(0)*1e-3 ); // leading large-R jet pT

        ++jentry;
    }
    f.Write();

    return 0;
}

