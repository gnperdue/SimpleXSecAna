//____________________________________________________________________________
/*!
  Example based on `gtestEventLoop` by Costas Andreopoulos.
  */
//____________________________________________________________________________

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

  TH1D *h_gammaE = new TH1D("h_gammaE", "Photon Energy", 50, 0, 0.05); 

  Hlist->Add(h_gammaE);
  
  
  //
  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);

    //
    // Put your event analysis code here 
    //
    Interaction *in = event.Summary();
    const ProcessInfo &proc = in->ProcInfo();
    const XclsTag &xclsv = in->ExclTag();
    bool qelcc = proc.IsQuasiElastic() && proc.IsWeakCC();
    bool charm = xclsv.IsCharmEvent();

    if (qelcc && !charm) {
      LOG("myAnalysis", pNOTICE) << event.XSec();
    }


    //
    // Loop over all particles in this event
    //

    GHepParticle * p = 0;
    TIter event_iter(&event);

    while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {

      if (p->Status() == kIStStableFinalState ) {
        if (p->Pdg() == 22) {
          LOG("myAnalysis", pNOTICE)  
            << "Got a : " << p->Name() << " with E = " << p->E() << " GeV";
          h_gammaE->Fill(p->E());
        }
      }

    }// end loop over particles	

    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  std::string histname = "$ANA_HIST_DIR/supernu.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();
  
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
