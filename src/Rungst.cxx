#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <TROOT.h>
#include <TApplication.h>
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "gstana.h"

Int_t flag = -1;

int main(int argc,char** argv) {

  TApplication theApp("App", &argc, argv);
  gApplication->Init();

  argc=theApp.Argc();
  argv=theApp.Argv();

  char *s = getenv( "ANA_LIST_DIR" );
  std::string filelist = s + std::string("/diffractive_events_file_list.txt");

  std::cout << "List file = " << filelist << std::endl;
  std::ifstream filestream( filelist.c_str() );
  TChain* chain = new TChain("gst");
  // Char_t fname[200];
  Char_t filename[200];
  while ( filestream >> filename ){
    // sprintf(fname, filename);
    std::cout << "  Adding " << filename << " to list..." << std::endl;
    chain->Add(filename);
  }

  std::cout << "Make new gstana pointer..." << std::endl;
  gstana * jjj = new gstana(chain);
  std::cout << "Starting loop..." << std::endl;
  jjj->Loop();
  std::cout << "Exiting loop..." << std::endl;
  std::cout << "Exiting Run Code..." << std::endl;
  return 0;
}
