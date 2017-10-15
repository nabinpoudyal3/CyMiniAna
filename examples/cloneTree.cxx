// Simple script to clone a tree (or more) from a file
// Example: File that contains ~100 systematics, but you 
//          just want to study the nominal tree on your laptop
//
// To compile:
//   g++ cloneTree.cxx -o cloneTree `root-config --cflags --glibs`
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"

#include <fstream>
#include <iostream>

#include <stdlib.h>


int main(){
  std::string path("/lustre/umt3/user/demarley/VLQ0L_AT2.4.13_dev/grid_trigger_downloads/");
  std::string dir("user.dmarley.302472.ProtosLHEFPythia8EvtGen.DAOD_EXOT7.e4112_s2608_r7772_r7676_p2669.AT2.4.13_dev_AllHad_trig_v0_output.root/");

  TChain * chain = new TChain("nominal");
  chain->Add( (path+dir+"user.dmarley*.root").c_str() );

  TChain * chain2 = new TChain("sumWeights");
  chain2->Add( (path+dir+"user.dmarley*.root").c_str() );

  TFile *file = TFile::Open("output.root","RECREATE");

  chain->CloneTree(-1,"fast");
  chain2->CloneTree(-1,"fast");

  file->Write();

  delete file;

  return 0;
}