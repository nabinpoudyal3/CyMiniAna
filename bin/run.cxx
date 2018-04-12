/*
Created:        --
Last Updated:    2 March 2018

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
#include <cstdlib>
#include <stdio.h>
#include <map>
#include <fstream>
#include <string>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <boost/algorithm/string/join.hpp>

#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/eventSelection.h"
#include "Analysis/CyMiniAna/interface/miniTree.h"
#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/histogrammer.h"
#include "Analysis/CyMiniAna/interface/efficiency.h"


int main(int argc, char** argv) {
    /* Steering macro for CyMiniAna */
    if (argc < 2) {
        cma::HELP();
        return -1;
    }

    unsigned long long maxEntriesToRun(0);                           // maximum number of entries in TTree
    unsigned int numberOfEventsToRun(0);                             // number of events to run
    bool passEvent(false);                                           // event passed selection

    // configuration
    configuration config(argv[1]);                                   // configuration file
    config.initialize();

    int nEvents         = config.nEventsToProcess();                 // requested number of events to run
    std::string outpathBase = config.outputFilePath();               // directory for output files
    std::string outpath = config.outputFilePath();                   // directory for output files
    unsigned long long firstEvent       = config.firstEvent();       // first event to begin running over
    std::vector<std::string> filenames  = config.filesToProcess();
    std::vector<std::string> treenames  = config.treeNames();
    std::vector<std::string> selections = config.selections();
    std::vector<std::string> cutfiles   = config.cutsfiles();
    std::string selection = boost::algorithm::join(selections, "-");

    bool generateCutsFiles = (cutfiles.size()!=selections.size());   // user did not provide different cuts files

    bool makeTTree        = config.makeTTree();
    bool makeHistograms   = config.makeHistograms();
    bool makeEfficiencies = config.makeEfficiencies();
    bool doSystWeights    = config.calcWeightSystematics();          // systemaics associated with scale factors

    std::string customDirectory( config.customDirectory() );
    if (customDirectory.length()>0  && customDirectory.substr(0,1).compare("_")!=0){
        customDirectory = "_"+customDirectory; // add '_' to beginning of string, if needed
    }

    // event selection(s) -- support for multiple event selections simulataneously
    std::vector<eventSelection> evtSels;
    std::vector<unsigned int> ncuts;                         // number of cuts in selection
    std::vector< std::vector<std::string> > namesOfCuts;     // names of cuts in selection
    for (unsigned int ss=0, size=selections.size(); ss<size; ss++) {
        std::string sel      = selections.at(ss);
        std::string cutsfile = (generateCutsFiles) ? "config/cuts_"+sel+".txt" : cutfiles.at(ss);

        eventSelection evtSel_tmp( config );
        evtSel_tmp.initialize( sel, cutsfile );

        evtSels.push_back(evtSel_tmp);
        ncuts.push_back(evtSel_tmp.numberOfCuts());
        namesOfCuts.push_back( evtSel_tmp.cutNames() );
    }


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
        struct stat dirBuffer;
        std::string outpath = outpathBase+"/"+selection+customDirectory;
        if ( !(stat((outpath).c_str(),&dirBuffer)==0 && S_ISDIR(dirBuffer.st_mode)) ){
            cma::DEBUG("RUN : Creating directory for storing output: "+outpath);
            system( ("mkdir "+outpath).c_str() );  // make the directory so the files are grouped together
        }

        std::size_t pos   = filename.find_last_of(".");     // the last ".", i.e., ".root"
        std::size_t found = filename.find_last_of("/");     // the last "/"
        std::string outputFilename = filename.substr(found+1,pos-1-found); // between "/" and "."
        // hopefully this returns: "diboson_WW" given something like:  "/some/path/to/file/diboson_WW.root"

        std::string fullOutputFilename = outpath+"/"+outputFilename+".root";
        std::unique_ptr<TFile> outputFile(TFile::Open( fullOutputFilename.c_str(), "RECREATE"));
        cma::INFO("RUN :   >> Saving to "+fullOutputFilename);

        // check the file type
        config.setFilename( filename );      // set the filename for the configuration
        config.inspectFile( *file );         // check the type of file this is

        std::vector<std::string> fileKeys;
        cma::getListOfKeys(file,fileKeys);   // keep track of ttrees in file

        histogrammer histMaker(config);      // initialize histogrammer
        efficiency effMaker(config);         // initialize efficiency class
        if (makeHistograms)
            histMaker.initialize( *outputFile,doSystWeights );
        if (makeEfficiencies)
            effMaker.bookEffs( *outputFile );

        for (auto& x : evtSels) x.setCutflowHistograms( *outputFile );  // setup cutflow histograms

        // -- Loop over treenames -> usually only one tree
        for (const auto& treename : treenames) {

            // check that the ttree exists in this file before proceeding
            if (std::find(fileKeys.begin(), fileKeys.end(), treename) == fileKeys.end()){
                cma::INFO("RUN : TTree "+treename+" is not present in this file, continuing to next TTree");
                continue;
            }


            // -- Load TTree to loop over
            cma::INFO("RUN :      TTree "+treename);
            TTreeReader myReader(treename.c_str(), file);

            // -- Make new Tree in Root file
            miniTree miniTTree(config);          // initialize TTree for new file
            if (makeTTree){
                // Setup subdirectories, if they exist
                std::string subdir("");
                std::size_t found = treename.find("/");
                if (found!=std::string::npos){
                    outputFile->cd();
                    subdir = treename.substr(0,found);
                    gDirectory->mkdir(subdir.c_str());
                }
                miniTTree.initialize( myReader.GetTree(), *outputFile, subdir );
            }

            // -- Number of Entries to Process -- //
            maxEntriesToRun = myReader.GetEntries(true);
            if (maxEntriesToRun<1) // skip files with no entries
                continue;

            numberOfEventsToRun = (nEvents<0 || ((unsigned int)nEvents+firstEvent)>maxEntriesToRun) ? maxEntriesToRun - firstEvent : nEvents;
            cma::INFO("RUN :      Processing "+std::to_string(numberOfEventsToRun)+" events ");

            // -- Event Loop -- //
            Long64_t imod = 1;                     // print to the terminal
            Event event = Event(myReader, config);

            Long64_t eventCounter = 0;    // counting the events processed
            Long64_t entry = firstEvent;  // start at a different event!
            while (myReader.Next()) {

                // Check number of events processed against number of events to run
                if (eventCounter+1 > numberOfEventsToRun){
                    cma::INFO("RUN :      Processed the desired number of events: "+std::to_string(eventCounter)+"/"+std::to_string(numberOfEventsToRun));
                    break;
                }

                // Update status on the console
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
                // can do separate cutflows by creating multiple instances of eventSelection()
                cma::DEBUG("RUN : Apply event selection");
                std::vector<unsigned int> passEvents;
                unsigned int passedEvents(0);
                for (unsigned int ss=0,size=selections.size();ss<size;ss++){
                    passEvent = evtSels.at(ss).applySelection(event);
                    passEvents.push_back( passEvent );
                    passedEvents += passEvent;
                }

                if (passedEvents>0){
                    // at least 1 selection passed
                    // share information on which selection passed in case these classes
                    // want to use that information, e.g., special branch or histogram name
                    cma::DEBUG("RUN : Passed selection, now reconstruct ttbar & save information");
                    event.ttbarReconstruction();

                    if (makeTTree)        miniTTree.saveEvent(event,passEvents);
                    if (makeHistograms)   histMaker.fill(event,passEvents);
                    if (makeEfficiencies) effMaker.fill(event,passEvents);
                }

                // iterate the entry and number of events processed
                ++entry;
                ++eventCounter;
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

    for (auto& evtSel : evtSels)
        evtSel.finalize();
    evtSels.clear();

    cma::INFO("RUN : *** End of file loop *** ");
}

// THE END
