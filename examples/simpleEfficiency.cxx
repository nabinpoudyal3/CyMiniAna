// Simple TEfficiency maker
// :: g++ simpleEfficiency.cxx `root-config --cflags --glibs` -o simpleEfficiency -lTreePlayer
// :: ./simpleEfficiency <filename>
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
    std::string outfile_name = "simpleEff_"+tmp_sample_name; // e.g., "trigger_eff_TTS_M1100.root"

    // outfile of TEfficiency histograms
    TFile f(outfile_name.c_str(),"recreate");

    // -- create a simple pT histogram to plot against efficiency
    auto h_pt0 = TH1F("signal_pt0", ";Leading Small-R Jet p_T;Events", 38,100,2000);

    // -- leading pT efficiency
    auto eff_pt0 = TEfficiency("eff_pt0",";Leading Small-R Jet p_T;#epsilon",38,100,2000);


    // -- load root file with data we are analyzing
    auto file = TFile::Open(sample_name.c_str());
    if (!file || file->IsZombie()) {
        return -1;
    }
    // -- Create a TTreeReader for the tree
    TTreeReader myReader("nominal", file); // "nominal" tree is the one we want


    // -- Access branches in the tree
    //    The triggers are type <char> (not <int> or <bool> as you might expect...)
    TTreeReaderValue<std::vector<float>> jet_pt(myReader,"jet_pt");
    TTreeReaderValue<char> hlt_j380(myReader,"HLT_j380");

    // -- Values for tracking whether or not the object passed/failed
    bool passed_j380(false); // boolean to check small-R jet pT trigger
    char success = 0x1;      // success value for 'char'

    // --------------------------------------------------------------
    // Loop over all entries in the input root file
    // --------------------------------------------------------------
    unsigned int jentry(0.);
    while (myReader.Next()) {

        if (jentry%100000==0)
            std::cout << "Entry: " << jentry << std::endl;

        // histogram of pT
        h_pt0.Fill( jet_pt->at(0)*1e-3 ); // leading small-R jet pT

        // booleans for triggers
        passed_j380 = (*hlt_j380)==success;

        //
        // Fill Efficiency histograms
        //

        // small-R jets
        eff_pt0.Fill(passed_j380, jet_pt->at(0)*1e-3 );


        ++jentry;
    }
    f.Write();

    return 0;
}

// THE END
