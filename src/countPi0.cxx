#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>
#include <TH2.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);

int    gOptNEvt;
string gOptInpFilename;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
    TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
    (int) tree->GetEntries();

  TObjArray * Hlist = new TObjArray(0);

  TH1D *h_npi0       = new TH1D("h_npi0","Number of Neutral Pions;Number;Counts",10,0,10); 
  TH1D *h_npipi0only = new TH1D("h_npipi0only","Number of Neutral Pions (no Charged Pions);Number;Counts",10,0,10); 

  Hlist->Add(h_npi0);
  Hlist->Add(h_npipi0only);

  //
  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);
    Interaction *in = event.Summary();

    //
    // Loop over all particles in this event
    //
    GHepParticle * p = 0;
    TIter event_iter(&event);

    int npi0 = 0;
    int nchgdpi = 0;

    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
      if (p->Status() == kIStStableFinalState ) {
        if (p->Pdg() == kPdgPi0) {
          npi0++;
        }
        if (p->Pdg() == kPdgPiP ||
            p->Pdg() == kPdgPiM) {
          nchgdpi++;
        }
      }

    }// end loop over particles	


    h_npi0->Fill(npi0);
    if (nchgdpi==0) {
      h_npipi0only->Fill(npi0);
    }

    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  std::string histname = "$ANA_HIST_DIR/countPi0.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();

  delete h_npi0;
  delete h_npipi0only;
  delete Hlist;

  // close input GHEP event file
  file.Close();

  LOG("myAnalysis", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("myAnalysis", pINFO) << "Parsing commad line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {
    LOG("myAnalysis", pINFO) 
      << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("myAnalysis", pFATAL) 
      << "Unspecified input filename - Exiting";
    exit(1);
  }

  // number of events to analyse
  if( parser.OptionExists('n') ) {
    LOG("myAnalysis", pINFO) 
      << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("myAnalysis", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }
}
//_________________________________________________________________________________
