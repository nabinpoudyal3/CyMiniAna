/*
Created:        --
Last Updated:   22 May 2017

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
#include "TSystem.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"

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

#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/Event.h"
#include "Analysis/CyMiniAna/interface/eventSelection.h"



int main(int argc, char** argv){
    /* Steering macro for skimming in CyMiniAna */
    if (argc < 2) {
        cma::HELP("skim");
        return -1;
    }

    configuration config(argv[1]);                   // configuration file
    config.initialize();

    // config arguments
    // -- for skimming, set the 'selection' to the kind of skimming you want
    std::vector<std::string> filenames = config.filesToProcess();
    std::vector<std::string> treenames = config.treeNames();
    std::string outpath  = config.outputFilePath();   // directory for output files
    int p_nEvents        = config.nEventsToProcess();       // requested number of events to run
    std::string skimMode = config.selection();        // kind of skimming to perform
    std::vector<std::string> skimModes;
    cma::split(skimMode, ',', skimModes);

    bool skimTtbarXS(false);
    bool skimDNN(false);
    bool skimPre(false);
    for (const auto& oneSkimMode : skimModes){
        if (oneSkimMode.compare("ttbarXS")==0) skimTtbarXS = true;
        else if (oneSkimMode.compare("DNN")==0) skimDNN = true;
        else if (oneSkimMode.compare("pre")==0) skimPre = true;    // AllHadAC skimming is the 'pre' selection
    }
    if (!skimTtbarXS && !skimDNN && !skimPre){
        cma::HELP("skim");
        cma::ERROR(" SKIM : The skimming choice "+skimMode+" is not supported. ");
        cma::ERROR(" SKIM : For skimming, must select from the following: ");
        cma::ERROR(" SKIM :    'ttbarXS', 'dnn', 'pre' ");
        cma::ERROR(" SKIM : as the 'selection' comma separated in the configuration file ");

        return -1;
    }
    cma::INFO("SKIM :   Skimming selected = "+skimMode );


    // Event selection
    cma::INFO("RUN : Setup event selection");
    eventSelection evtSel( config );                           // initialize event selection
    evtSel.initialize();

    unsigned int ncuts = evtSel.numberOfCuts();                // number of cuts in selection
    std::vector<std::string> cutNames = evtSel.cutNames();     // names of cuts


    // -- Loop over Filenames
    for (const auto& filename : filenames){

        cma::INFO("SKIM :   Opening "+filename);

        auto oldfile = TFile::Open(filename.c_str());
        if (!oldfile || oldfile->IsZombie()){
            cma::WARNING("SKIM :  -- File: "+filename);
            cma::WARNING("SKIM :     does not exist or it is a Zombie. ");
            cma::WARNING("SKIM :     Continuing to next file. ");
            continue;
        }
        config.setFilename(filename);
        config.checkFileType( *oldfile );
        bool isMC = config.isMC();

        std::vector<std::string> fileKeys;
        cma::getListOfKeys(oldfile,fileKeys); // keep track of ttrees in file

        // clone the truth tree (if MC sample) into each outputfile
        TTree *ttbarXSTruthTree;
        TTree *dnnTruthTree;
        TTree *preTruthTree;
        TTree *oldTruthTree;
        if (isMC) oldTruthTree = (TTree*)oldfile->Get("truth");

        // clone the sum of weights ttree into each outputfile
        TTree *ttbarXSSumWeights;
        TTree *dnnSumWeights;
        TTree *preSumWeights;
        TTree *oldSumWeights = (TTree*)oldfile->Get("sumWeights");


        // -- Output file -- //
        // directory to put all files under one samples in same directory (don't want all files in a single directory together)
        std::vector<std::string> directories;
        cma::split(filename,'/',directories);          // split path by '/' delimeter
        std::string originalDir("/");                  // keep the '/' in case this isn't needed
        if (directories.size()>1){
            originalDir = "/"+directories.rbegin()[1]; // next-to-last element in vector
        }
        // check if the directory already exists
        // http://stackoverflow.com/questions/3828192/checking-if-a-directory-exists-in-unix-system-call
        struct stat dirBuffer;
        if ( !(stat((outpath+originalDir).c_str(),&dirBuffer)==0 && S_ISDIR(dirBuffer.st_mode)) ){
            system( ("mkdir "+outpath+originalDir).c_str() );  // make the directory so the files are grouped together
        }
        // make the new outputfile
        std::size_t pos   = filename.find_last_of(".");     // the last ".", i.e., ".root" (sometimes '.2' or '.1')
        std::size_t found = filename.find_last_of("/");     // the last "/"
        std::string outfilename = filename.substr(found+1,pos-1-found); // betwee "/" and "."

        std::size_t rootEnding = outfilename.find(".root"); // remove '.root' if it is still there
        if (rootEnding!=std::string::npos){
            outfilename = outfilename.substr(0,rootEnding);
        }

        std::string outputFileName(outpath+originalDir+"/"+outfilename+"_");
        cma::INFO("SKIM :   >> Saving to "+outputFileName);
        for (const auto& skimMode : skimModes){
            cma::INFO("SKIM :             -- "+skimMode);
        }

        // -- Cutflow histograms
        std::map<std::string,TH1D*> h_cutflow_pre;
        std::map<std::string,TH1D*> h_cutflow_pre_unweighted;

        // make new files -- multiple skims at once!
        std::unique_ptr<TFile> outputFileTtbarXS;
        std::unique_ptr<TFile> outputFileDNN;
        std::unique_ptr<TFile> outputFilePre;

        if (skimPre){
            TFile* preFile = new TFile( (outputFileName+"pre.root").c_str(), "RECREATE");
            outputFilePre  = std::unique_ptr<TFile>(preFile);
            preSumWeights  = oldSumWeights->CloneTree(-1,"fast");        // sum of weights
            if (isMC) preTruthTree = oldTruthTree->CloneTree(-1,"fast"); // truth
        }
        if (skimTtbarXS){
            TFile* ttbarXSFile = new TFile( (outputFileName+"ttbarXS.root").c_str(), "RECREATE");
            outputFileTtbarXS  = std::unique_ptr<TFile>(ttbarXSFile);
            ttbarXSSumWeights  = oldSumWeights->CloneTree(-1,"fast");        // sum of weights
            if (isMC) ttbarXSTruthTree = oldTruthTree->CloneTree(-1,"fast"); // truth
        }
        if (skimDNN){
            TFile* dnnFile = new TFile( (outputFileName+"dnn.root").c_str(), "RECREATE");
            outputFileDNN  = std::unique_ptr<TFile>(dnnFile);
            dnnSumWeights  = oldSumWeights->CloneTree(-1,"fast");        // sum of weights
            if (isMC) dnnTruthTree = oldTruthTree->CloneTree(-1,"fast"); // truth
        }


        // -- Loop over TTrees
        for (const auto& treename : treenames) {

            // check that the ttree exists in this file before proceeding
            if (std::find(fileKeys.begin(), fileKeys.end(), treename) == fileKeys.end()){
                cma::INFO("SKIM : TTree "+treename+" is not present in this file, continuing to next TTree");
                continue;
            }

            // -- Make new Tree in Root file
            cma::INFO("SKIM :      TTree "+treename);

            TTree *originalTree = (TTree*)oldfile->Get(treename.c_str());

            Long64_t numberOfEventsToRun(0);
            Long64_t maxEntriesToRun = originalTree->GetEntries();
            if (p_nEvents < 0 || p_nEvents > maxEntriesToRun) numberOfEventsToRun = maxEntriesToRun;
            else numberOfEventsToRun = p_nEvents;

            if (maxEntriesToRun<1 || originalTree->IsZombie()){
                cma::WARNING("SKIM : TTree "+treename+" has no entries or is a Zombie, continuing to next TTree");
                continue;
            }
            config.setTreename(treename);

            // -- Cutflows (only supporting 'pre' right now)
            if (skimPre){
                outputFilePre->cd();

                h_cutflow_pre[treename] = new TH1D( (treename+"_cutflow").c_str(),(treename+"_cutflow").c_str(),ncuts+1,0,ncuts+1);
                h_cutflow_pre_unweighted[treename] = new TH1D( (treename+"_cutflow_unweighted").c_str(),(treename+"_cutflow_unweighted").c_str(),ncuts+1,0,ncuts+1);

                // label bins with cut names (from cut files)
                h_cutflow_pre[treename]->GetXaxis()->SetBinLabel(1,"INITIAL");
                h_cutflow_pre_unweighted[treename]->GetXaxis()->SetBinLabel(1,"INITIAL");
                for (unsigned int c=1;c<=ncuts;++c){
                    h_cutflow_pre[treename]->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
                    h_cutflow_pre_unweighted[treename]->GetXaxis()->SetBinLabel(c+1,cutNames.at(c-1).c_str());
                }
            }

            // -- Load variables for skimming -- //

            // TRIGGERS -- A
            char HLT_j360_a10r_L1J100 = 0;
            char HLT_j420_a10r_L1J100 = 0;
            char HLT_ht1000_L1J100    = 0;

            // AnalysisTop SELECTIONS
            int jet_2015 = 0;
            int jet_2016 = 0;

            // LARGE-R JETS
            std::vector<float> *ljet_pt  = 0;
            std::vector<float> *ljet_eta = 0;

            // -- SetBranchAddresses
            originalTree->SetBranchAddress("HLT_j360_a10r_L1J100",&HLT_j360_a10r_L1J100);
            originalTree->SetBranchAddress("HLT_j420_a10r_L1J100",&HLT_j420_a10r_L1J100);
            originalTree->SetBranchAddress("HLT_ht1000_L1J100",&HLT_ht1000_L1J100);

            originalTree->SetBranchAddress("jets_2015",&jet_2015);
            originalTree->SetBranchAddress("jets_2016",&jet_2016);

            originalTree->SetBranchAddress("ljet_pt",&ljet_pt);
            originalTree->SetBranchAddress("ljet_eta",&ljet_eta);

            // -- New TTree which represents skimmed results
            TTree *ttbarXSTree = 0;
            if (skimTtbarXS){
                outputFileTtbarXS->cd();
                ttbarXSTree = originalTree->CloneTree(0);
            }

            TTree *dnnTree = 0;
            if (skimDNN){
                outputFileDNN->cd();
                dnnTree = originalTree->CloneTree(0);
            }

            TTree *preTree = 0;
            if (skimPre){
                outputFilePre->cd();
                preTree = originalTree->CloneTree(0);
            }


            // -- Event loop -- //
            TTreeReader myReader( (treename).c_str(), oldfile );
            Event event = Event(myReader, config);  // to load Event objects

            cma::INFO("SKIM :      Number of events in original tree "+std::to_string(maxEntriesToRun) );
            Long64_t imod = 1;
            for (Long64_t i=0;i<numberOfEventsToRun; i++) {
                if (i%imod==0){
                    cma::INFO("SKIM :      Processing event "+std::to_string(i) );
                    if(imod<2e4) imod *=10;
                }
                originalTree->GetEntry(i);
                event.execute(i);

                // -- Setup skimming
                bool passSkim(false);
                unsigned int passedCuts(0);

                // https://cds.cern.ch/record/2198341/files/ATL-COM-PHYS-2016-902.pdf
                if (skimTtbarXS){
                    if (HLT_j360_a10r_L1J100==0x1 ||  // 2015: HLT_j360_a10r_L1J100
                        HLT_j420_a10r_L1J100==0x1 ||  // 2016: HLT_j400_a10r_L1J100
                        HLT_ht1000_L1J100==0x1 ){     // both: HLT_ht1000 
                        ++passedCuts;
                    }
                    if (ljet_pt->size()>1){
                        if (ljet_pt->at(0)>500000. && ljet_pt->at(1)>350000.)
                            ++passedCuts;
                    }

                    // Fill the tree
                    if (passedCuts>1){        // 1 cuts for skimming
                        ttbarXSTree->Fill();
                    }
                }

                // Just events for making DNN
                if (skimDNN){
                    if (ljet_pt->size()>1){
                        if (ljet_pt->at(0)>300000. && ljet_pt->at(1)>300000.) dnnTree->Fill();
                    }
                }

                // Use eventSelection class for AllHadAC "pre" skimming
                if (skimPre){
                    cma::DEBUG("RUN : Apply pre event selection ");
                    passSkim = evtSel.applySelection(event,*h_cutflow_pre.at( treename ),*h_cutflow_pre_unweighted.at( treename ) );
                    if (passSkim) preTree->Fill();
                }
            } // end event loop

            event.finalize();
            if (skimTtbarXS){
                cma::INFO("SKIM :      Number of events in skimTtbarXS tree "+std::to_string(ttbarXSTree->GetEntries()) );
            }
            if (skimDNN){
                cma::INFO("SKIM :      Number of events in skimDNNS tree "+std::to_string(dnnTree->GetEntries()) );
            }
            if (skimPre){
                cma::INFO("SKIM :      Number of events in skimPreS tree "+std::to_string(preTree->GetEntries()) );
            }
        } // end tree loop

        cma::INFO("SKIM : Write output file");
        if (skimTtbarXS){
            outputFileTtbarXS->Write();
            outputFileTtbarXS->Close();
        }
        if (skimDNN){
            outputFileDNN->Write();
            outputFileDNN->Close();
        }
        if (skimPre){
            outputFilePre->Write();
            outputFilePre->Close();
        }

        cma::DEBUG("SKIM : Delete old file");
        delete oldfile;   // delete other files (limits "ROOT: too many open files" error)
        oldfile = ((TFile *)0);
    } // end file loop
    cma::INFO("SKIM : Finished.");

    return 0;
}

// the end
