/*
Created:        --
Last Updated:   20 August    2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

Basic steering macro for running CyMiniAna
 - Make new ntuples
 - Make histograms/efficiencies

*/
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <map>
#include <fstream>
#include <string>
#include <math.h>

#include "diHiggs/CyMiniAna/interface/configuration.h"
#include "diHiggs/CyMiniAna/interface/Event.h"
#include "diHiggs/CyMiniAna/interface/eventSelection.h"
#include "diHiggs/CyMiniAna/interface/miniTree.h"
#include "diHiggs/CyMiniAna/interface/tools.h"
#include "diHiggs/CyMiniAna/interface/histogrammer.h"
#include "diHiggs/CyMiniAna/interface/efficiency.h"


int main(int argc, char** argv) {
    /* Steering macro for CyMiniAna */
    if (argc < 2) {
        cma::HELP();
        return -1;
    }

    unsigned int maxEntriesToRun(0);     // maximum number of entries in TTree
    unsigned int numberOfEventsToRun(0); // number of events to run
    bool passEvent(false);   // event passed selection
    bool isMC(false);      // Running over MC file

    // configuration
    configuration config(argv[1]);                         // configuration file
    config.initialize();

    int p_nEvents       = config.nEventsToProcess(); // requested number of events to run
    std::string outpath = config.outputFilePath();   // directory for output files
    std::vector<std::string> filenames = config.filesToProcess();
    std::vector<std::string> treenames = config.treeNames();
    std::string selection(config.selection());

    bool makeNewFile      = config.makeNewFile();
    bool makeHistograms   = config.makeHistograms();
    bool makeEfficiencies = config.makeEfficiencies();
    bool doSystWeights    = config.calcWeightSystematics(); // systemaics associated with scale factors

    std::string customFileEnding( config.customFileEnding() );
    if (customFileEnding.length()>0  && customFileEnding.substr(0,1).compare("_")!=0){
        customFileEnding = "_"+customFileEnding; // add '_' to beginning of string, if needed
    }

    // event selection
    eventSelection evtSel( config );
    evtSel.initialize();
    unsigned int ncuts = evtSel.numberOfCuts();            // number of cuts in selection
    std::vector<std::string> cutNames = evtSel.cutNames(); // names of cuts


    // --------------- //
    // -- File loop -- //
    // --------------- //
    unsigned int numberOfFiles(filenames.size());
    unsigned int currentFileNumber(0);
    cma::INFO("RUN : *** Starting file loop *** ");
    for (const auto& filename : filenames) {

        ++currentFileNumber;
        cma::INFO("RUN :   Opening "+filename+"   ("+std::to_string(currentFileNumber)+"/"+std::to_string(numberOfFiles)+")");

        auto file = TFile::Open(filename.c_str());
        if (!file || file->IsZombie()){
            cma::WARNING("RUN :  -- File: "+filename);
            cma::WARNING("RUN :     does not exist or it is a Zombie. ");
            cma::WARNING("RUN :     Continuing to next file. ");
            continue;
        }

        // -- Output file -- //
        std::size_t pos   = filename.find_last_of(".");     // the last ".", i.e., ".root"
        std::size_t found = filename.find_last_of("/");     // the last "/"
        std::string outputFilename = filename.substr(found+1,pos-1-found); // betwee "/" and "."
        // hopefully this returns: "diboson_WW_361082" given something like:
        // "/some/path/to/file/diboson_WW_361082.root"

        std::string fullOutputFilename = outpath+"/"+outputFilename+"_"+selection+customFileEnding+".root";
        std::unique_ptr<TFile> outputFile(TFile::Open( fullOutputFilename.c_str(), "RECREATE"));
        cma::INFO("RUN :   >> Saving to "+fullOutputFilename);

        // check the file type
        isMC = config.isMC( *file );

        std::vector<std::string> fileKeys;
        cma::getListOfKeys(file,fileKeys); // keep track of ttrees in file

/*
        // -- Extra trees to just clone into new file
        TChain * truth_chain;
        TChain * sumWeights_chain;

        if (makeNewFile){
            if (isMC){
                truth_chain = new TChain("truth");
                truth_chain->Add( filename.c_str() );
                truth_chain->CloneTree(-1,"fast");      // clone the truth tree
            }

            sumWeights_chain = new TChain("sumWeights");
            sumWeights_chain->Add( filename.c_str() );
            sumWeights_chain->CloneTree(-1,"fast");    // clone the sumWeights tree
        }
        // -- End tree cloning
*/

        histogrammer histMaker(config);      // initialize histogrammer
        efficiency effMaker(config);         // initialize efficiency class
        if (makeHistograms)
            histMaker.initialize( *outputFile,doSystWeights );
        if (makeEfficiencies)
            effMaker.bookEffs( *outputFile );

        // -- Cutflow histograms
        std::map<std::string, TH1D*>  h_cutflows;            // map of cutflow histograms (weights applied)
        std::map<std::string, TH1D*>  h_cutflows_unweighted; // map of cutflow histograms (raw # of events)

        // -- Loop over treenames
        for (const auto& treename : treenames) {

            // check that the ttree exists in this file before proceeding
            if (std::find(fileKeys.begin(), fileKeys.end(), treename) == fileKeys.end()){
                cma::INFO("SKIM : TTree "+treename+" is not present in this file, continuing to next TTree");
                continue;
            }

            // -- Cutflow histogram [initialize and label bins]
            h_cutflows[treename] = new TH1D( (treename+"_cutflow").c_str(),(treename+"_cutflow").c_str(),ncuts+1,0,ncuts+1);
            h_cutflows_unweighted[treename] = new TH1D( (treename+"_cutflow_unweighted").c_str(),(treename+"_cutflow_unweighted").c_str(),ncuts+1,0,ncuts+1);

            h_cutflows[treename]->GetXaxis()->SetBinLabel(1,"INITIAL");
            h_cutflows_unweighted[treename]->GetXaxis()->SetBinLabel(1,"INITIAL");

            for (unsigned int c=1;c<=ncuts;++c){
                h_cutflows[treename]->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
                h_cutflows_unweighted[treename]->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
            }

            // -- Load TTree to loop over
            cma::INFO("RUN :      TTree "+treename);
            TTreeReader myReader(treename.c_str(), file);

            // -- Make new Tree in Root file

            miniTree miniTTree(config);          // initialize TTree for new file
            if (makeNewFile)
                miniTTree.initialize( myReader.GetTree(), *outputFile );

            // -- Number of Entries -- //
            maxEntriesToRun = myReader.GetEntries(true);
            if (maxEntriesToRun<1) // skip files with no entries
                continue;

            if (p_nEvents < 0 || (unsigned int)p_nEvents > maxEntriesToRun)
                numberOfEventsToRun = maxEntriesToRun;
            else
                numberOfEventsToRun = p_nEvents;

            // -- Event Loop -- //
            Long64_t imod = 1;                     // print to the terminal
            Event event = Event(myReader, config);

            Long64_t entry = 0;
            while (myReader.Next()) {

                if (entry+1 > numberOfEventsToRun)
                    break;
                if (entry%imod==0){
                    cma::INFO("RUN :       Processing event "+std::to_string(entry) );
                    if(imod<2e4) imod *=10;
                }

                // -- Build Event -- //
                cma::DEBUG("RUN : Execute event");
                event.execute(entry);
                // now we have event object that has the event-level objects in it
                // pass this to the selection tools

                // -- Event Selection -- //
                // can do separate cutflows for each lepton flavor (just add if statement)
                // can do separate event selections as well, e.g.,
                //   evtSel.applySelectionBoosted(...) || evtSel.applySelectionResolved(...)
                cma::DEBUG("RUN : Apply event selection");
                passEvent = evtSel.applySelection(event,*h_cutflows.at( treename.c_str() ),*h_cutflows_unweighted.at( treename.c_str() ));

                if (passEvent){
                    cma::DEBUG("RUN : Passed selection, now save information");
                    if (makeNewFile)       miniTTree.saveEvent(event);
                    if (makeHistograms)    histMaker.fill( event );
                    if (makeEfficiencies)  effMaker.fill( event );
                }

                ++entry;
            } // end event loop

            event.finalize();
            miniTTree.finalize();
        } // end tree loop

        // put overflow/underflow content into the first and last bins
        histMaker.overUnderFlow();

        cma::INFO("RUN :   END Running  "+filename);
        cma::INFO("RUN :   >> Output at "+fullOutputFilename);

        outputFile->Write();
        outputFile->Close();

        // -- Clean-up stuff
        delete file;          // free up some memory 
        file = ((TFile *)0);  // (no errors for too many root files open)
    } // end file loop

    cma::INFO("RUN : *** End of file loop *** ");
}

// THE END
