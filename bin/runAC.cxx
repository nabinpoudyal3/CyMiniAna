/*
Created:        16 February  2016
Last Updated:   20 August    2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Steering macro for running CyMiniAna
 - Make new ntuples
 - Make histograms

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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/eventSelection.h"
#include "Analysis/CyMiniAna/interface/miniTree.h"
#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/histogrammer.h"
#include "Analysis/CyMiniAna/interface/histogrammerTruth.h"
#include "Analysis/CyMiniAna/interface/efficiency.h"
#include "Analysis/CyMiniAna/interface/efficiencyTruth.h"


int main(int argc, char** argv) {
    /* Steering macro for CyMiniAna */
    if (argc < 2) {
        cma::HELP();
        return -1;
    }

    unsigned int maxEntriesToRun(0);     // maximum number of entries in TTree
    unsigned int numberOfEventsToRun(0); // number of events to run
    bool passEvent(false);    // event passed selection
    bool isMC(false);         // Running over MC file
    bool qcdSelection(false); // qcd regions

    configuration config(argv[1]);                         // configuration file
    config.initialize();

    cma::INFO("RUN : Configuration initialized");

    // -- Config arguments
    int p_nEvents = config.nEventsToProcess(); // requested number of events to run
    std::string outpathBase = config.outputFilePath();   // directory for output files
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

    // Event selection(s)
    std::vector<eventSelection> evtSels;
    unsigned int ncuts(0);
    std::vector<std::string> cutNames;
    std::vector<std::string> regions;

    std::map<std::string,bool> passQCD;
    if (selection.compare("qcd")==0){
        // this setup assumes applying selections across many regions
        // -- all regions have the same number and names of cuts
        //    only the values are different
        qcdSelection = true;
        regions = config.qcdSelections();
        for (const auto& region : regions){
            eventSelection evtSel_tmp( config, region );
            evtSel_tmp.initialize( "share/cuts_qcd.txt" ); // just need the cut names, not actual values
            evtSels.push_back(evtSel_tmp);
            passQCD[region] = false;                       // for saving to TTree later
        }
    }
    else{
        regions.push_back(config.selection());
        eventSelection evtSel_tmp( config );               // initialize event selection
        evtSel_tmp.initialize();
        evtSels.push_back(evtSel_tmp);
    }

    ncuts    = evtSels.at(0).numberOfCuts();
    cutNames = evtSels.at(0).cutNames();


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
        // CMS doesn't use 'mcChannelNumber', need to keep the same file names
        // therefore, make new directories for the different selections
        struct stat dirBuffer;
        std::string outpath = outpathBase+"/"+selection+customFileEnding;
        if ( !(stat((outpath).c_str(),&dirBuffer)==0 && S_ISDIR(dirBuffer.st_mode)) ){
            cma::DEBUG("RUN : Creating directory for storing output: "+outpath);
            system( ("mkdir "+outpath).c_str() );  // make the directory so the files are grouped together
        }

        std::size_t pos   = filename.find_last_of(".");     // the last ".", i.e., ".root"
        std::size_t found = filename.find_last_of("/");     // the last "/"
        std::string outputFilename = filename.substr(found+1,pos-1-found); // betwee "/" and "."
        // hopefully this returns: "diboson_WW_361082" given something like:
        // "/some/path/to/file/diboson_WW_361082.root"

        std::string fullOutputFilename = outpath+"/"+outputFilename+".root";
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
            cma::DEBUG("RUN : Clone truth and sumWeights trees");
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

        if(config.doRecoEventLoop()) {
            // -- Cutflow histograms
            std::map<std::string, std::vector<TH1D*> > h_cutflows;            // map of cutflow histograms (weights applied)
            std::map<std::string, std::vector<TH1D*> > h_cutflows_unweighted; // map of cutflow histograms (raw # of events)

            std::vector<histogrammer> histMakers;
            if (makeHistograms){
                cma::DEBUG("RUN : Book histograms");
                for (const auto& region : regions){
                    histogrammer histMaker_tmp(config,region);  // initialize histogrammer for this region
                    histMaker_tmp.initialize( *outputFile,doSystWeights );
                    histMakers.push_back(histMaker_tmp);
                }
            }

            efficiency effMaker(config);                         // initialize efficiency class (only signal region considered for now)
            if (makeEfficiencies) effMaker.bookEffs( *outputFile );

            // -- Loop over treenames
            cma::DEBUG("RUN : Begin loop over TTrees");
            for (const auto& treename : treenames) {

                // check that the ttree exists in this file before proceeding
                if (std::find(fileKeys.begin(), fileKeys.end(), treename) == fileKeys.end()){
                    cma::INFO("SKIM : TTree "+treename+" is not present in this file, continuing to next TTree");
                    continue;
                }

                // -- Cutflow histogram [initialize and label bins]
                h_cutflows[treename].clear();
                h_cutflows_unweighted[treename].clear();

                for (auto& reg : regions) {
                    TH1D* h_cutflows_tmp     = new TH1D( (treename+"_"+reg+"_cutflow").c_str(),(treename+"_"+reg+"_cutflow").c_str(),ncuts+1,0,ncuts+1);
                    TH1D* h_cutflows_tmp_unw = new TH1D( (treename+"_"+reg+"_cutflow_unweighted").c_str(),(treename+"_"+reg+"_cutflow_unweighted").c_str(),ncuts+1,0,ncuts+1);

                    // label bins with cut names (from cut files)
                    h_cutflows_tmp->GetXaxis()->SetBinLabel(1,"INITIAL");
                    h_cutflows_tmp_unw->GetXaxis()->SetBinLabel(1,"INITIAL");

                    for (unsigned int c=1;c<=ncuts;++c){
                        h_cutflows_tmp->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
                        h_cutflows_tmp_unw->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
                    }

                    h_cutflows[treename].push_back(h_cutflows_tmp);
                    h_cutflows_unweighted[treename].push_back(h_cutflows_tmp_unw);
                }


                // -- Load TTree to loop over
                cma::INFO("RUN :      TTree "+treename);
                TTreeReader myReader(treename.c_str(), file);

                // -- Make new Tree in Root file
                miniTree miniTTree(config);                  // initialize TTree for new file
                if (makeNewFile)      miniTTree.initialize(myReader.GetTree(),*outputFile);

                // -- Number of Entries -- //
                maxEntriesToRun = myReader.GetEntries(true);
                // skip files with no entries
                if (maxEntriesToRun<1) continue;

                if (p_nEvents < 0 || (unsigned int)p_nEvents > maxEntriesToRun){
                    numberOfEventsToRun = maxEntriesToRun;
                }
                else{
                    numberOfEventsToRun = p_nEvents;
                }
                cma::INFO("RUN :      #events to loop over: "+std::to_string(numberOfEventsToRun));

                // -- Event Loop -- //
                Long64_t imod = 1;                     // print to the terminal
                Event event = Event(myReader, config);

                Long64_t entry = 0;
                while (myReader.Next()) {

                    if (entry+1 > numberOfEventsToRun) break;
                    if (entry%imod==0){
                        cma::INFO("RUN :      Processing event "+std::to_string(entry) );
                        if(imod<2e4) imod *=10;
                    }

                    // -- Build Event -- //
                    cma::DEBUG("RUN : Execute event");
                    event.execute(entry);
                    // now we have the 'event' that has the event-level objects in it
                    // pass this to the selection tools

                    // -- Event Selection -- //
                    // can do separate cutflows & event selections as well, e.g.,
                    //   evtSelBoosted.applySelection(...) || evtSelResolved.applySelection(...)
                    // for more than a few different selections (or 'regions') recommend using vectors to hold these
                    // rather than coding them individually
                    unsigned int idx(0);
                    passEvent = false;
                    bool passSingleEvent(false);

                    cma::DEBUG("RUN : Apply event selection");
                    for (auto& sel : evtSels){
                        passSingleEvent = sel.applySelection(event, *h_cutflows[treename].at(idx),*h_cutflows_unweighted[treename].at(idx));
                        if( passSingleEvent ){
                            passEvent = true;
                            if (makeHistograms) histMakers.at(idx).fill( event ); // fill relevant histograms if this selection successful
                        }
                        if (qcdSelection) passQCD[regions.at(idx)] = passSingleEvent;  // save booleans of qcd regions
                        idx++;
                    }
                    if (!qcdSelection){
                        for (const auto& x : config.qcdSelections()) passQCD[x] = false;
                    }

                    if (passEvent){
                        cma::DEBUG("RUN : Pass event, now save information");
                        if (makeNewFile) miniTTree.saveEvent(event);
                        if (makeEfficiencies)  effMaker.fill( event );
                    }
                    ++entry;
                } // end event loop

                event.finalize();
                miniTTree.finalize();
            } // end tree loop

            // put overflow/underflow bins in the last/first bins
            if (makeHistograms){
                cma::DEBUG("RUN : Fill histograms with overflow & underflow content");
                for (unsigned int reg=0,size=regions.size();reg<size;reg++){
                    histMakers.at(reg).overUnderFlow();
                }
            }
        } // end if running reco event loop

        /**** LOOP OVER TRUTH TREE FOR TTBAR SAMPLES WITH TRUTH TREE (MUST HAVE TOP PARTONS HISTORY!) ****/
        // FOR NOW ONLY DO ACCEPTANCE STUDIES ON NOMINAL TTREE
        if(config.useTruth() && config.doTruthEventLoop()) {

            TTreeReader truthReader("truth", file);
            TTreeReader recoReader("nominal", file);

            config.setMatchTruthToReco(false); // switch config option so that we run over truth events

            // -- Number of Entries -- //
            numberOfEventsToRun = truthReader.GetEntries(true);

            // skip files with no entries
            if (numberOfEventsToRun<1) continue;
            cma::INFO("RUN : TRUTH :     TTree truth");
            cma::INFO("RUN : TRUTH :     #events to loop over: "+std::to_string(numberOfEventsToRun));

            // event selector -- only signal region for now
            eventSelection evtSel(config, "signal");
            evtSel.initialize("share/cuts_signal.txt");

            // Objects for making truth efficiencies / histograms
            efficiencyTruth effMaker(config);
            if (makeEfficiencies) effMaker.bookEffs( *outputFile );

            histogrammerTruth histMaker(config);
            if (makeHistograms) histMaker.initialize( *outputFile );

            // cutflow histograms for truth
            ncuts    = evtSel.numberOfCuts();
            cutNames = evtSel.cutNames();
            TH1D* h_cutflow_truth = new TH1D("truth_signal_cutflow", "truth_signal_cutflow", ncuts+1, 0, ncuts+1);
            TH1D* h_cutflow_truth_unw = new TH1D("truth_signal_cutflow_unweighted", "truth_signal_cutflow_unweighted", ncuts+1, 0, ncuts+1);

            // label bins with cut names (from cut files)
            h_cutflow_truth->GetXaxis()->SetBinLabel(1,"INITIAL");
            h_cutflow_truth_unw->GetXaxis()->SetBinLabel(1,"INITIAL");

            for (unsigned int c=1;c<=ncuts;++c){
                h_cutflow_truth->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
                h_cutflow_truth_unw->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
            }

            // -- Event Loop -- //
            Long64_t imod = 1;                     // print to the terminal
            Event event = Event(recoReader, config);

            Long64_t entry = 0;
            while (truthReader.Next()) {

                if (entry+1 > numberOfEventsToRun) break;
                if (entry%imod==0){
                    cma::INFO("RUN :      Processing event "+std::to_string(entry) );
                    if(imod<2e4) imod *=10;
                }

                // -- Build Event -- //
                cma::DEBUG("RUN : TRUTH : Execute event");
                event.execute(entry);
                // this event object has initialized truth, and initialized reco level
                // *if* and only *if* there is a matching reco event to truth event
                // pass this to the selection tools

                // -- Event Selection (only signal region for acceptance studies now) -- //
                passEvent = false;
                cma::DEBUG("RUN : Apply event selection");

                // TODO we also need to create some fiducial selection (like truth mttbar cut)
                // fill underflow bin of cutflow with number of events in truth
                h_cutflow_truth_unw->Fill(-0.5);
                h_cutflow_truth->Fill(-0.5, event.truth_weight_mc());

                passEvent = evtSel.applySelection(event, *h_cutflow_truth, *h_cutflow_truth_unw);

                if (passEvent)
                    cma::DEBUG("RUN : TRUTH : Event "+std::to_string(entry)+" passed selection.");
                else 
                    cma::DEBUG("RUN : TRUTH : Event "+std::to_string(entry)+" failed selection.");

                if(makeEfficiencies) effMaker.fill( event, passEvent );
                if(makeHistograms)   histMaker.fill( event );

                ++entry;
            } // end event loop

            histMaker.overUnderFlow();
        }

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

// The End
