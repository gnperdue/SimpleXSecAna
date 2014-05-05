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
#include "grooana.h"

Int_t flag = -1;

int main(int argc,char** argv) {

  TApplication theApp("App", &argc, argv);
  gApplication->Init();

  argc=theApp.Argc();
  argv=theApp.Argv();

  char *s = getenv( "ANA_LIST_DIR" );
  std::string filelist = s + std::string("/gtrac_list.txt");

  std::cout << "List file = " << filelist << std::endl;
  std::ifstream filestream( filelist.c_str() );
  TChain* chain = new TChain("gRooTracker");
  Char_t fname[200];
  Char_t filename[200];
  while ( filestream >> filename ){
    sprintf(fname, "$ANA_DATA_DIR/%s",filename);
    std::cout << "  Adding " << fname << " to list..." << std::endl;
    chain->Add(fname);
  }

  std::cout << "Make new grooana pointer..." << std::endl;
  grooana * jjj = new grooana(chain);
  std::cout << "Starting loop..." << std::endl;
  jjj->Loop();
  std::cout << "Exiting loop..." << std::endl;
  std::cout << "Exiting Run Code..." << std::endl;
  return 0;
}
