//____________________________________________________________________________
#include <string>
#include <fstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TIterator.h>
#include <TH1.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Conventions/Constants.h"

using std::string;
using namespace genie;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //
  // Prepare and load the chain of GHEP files
  //
  char *s = getenv( "ANA_LIST_DIR" );
  std::string filelist = s + std::string("/diffractive_events_file_list.txt");

  std::cout << "List file = " << filelist << std::endl;
  std::ifstream filestream( filelist.c_str() );
  TChain* chain = new TChain("gtree");
  Char_t filename[200];
  while ( filestream >> filename ){
    std::cout << "  Adding " << filename << " to list..." << std::endl;
    chain->Add(filename);
  }

  Long64_t nEntries = chain->GetEntries(); 
  std::cout << "There are " << nEntries << " in the chain." << std::endl;

  NtpMCEventRecord * mcrec = 0;
  chain->SetBranchAddress("gmcrec", &mcrec);

  //
  // Set up Histograms
  //
  TObjArray * Hlist = new TObjArray(0);

  TH1D *h_proton_e = new TH1D("h_proton_e",
      "Proton Energy;GeV;Counts per 25 MeV",40,0.8,1.8); 
  TH1D *h_proton_ke = new TH1D("h_proton_ke",
      "Proton Energy;MeV;Counts per 4 MeV",100,0,400); 
  TH1D *h_pi0_e = new TH1D("h_pi0_e",
      "Neutral Pion Energy;GeV;Counts per 100 MeV",50,0,5); 
  TH1D *h_t = new TH1D("h_t",
      "t;GeV^{2};Counts per 0.02 GeV^{2}",50,0,1); 
  TH1D *h_t_over_proton_ke = new TH1D("h_t_over_proton_ke",
      "t/T_{p};GeV;Counts per 40 MeV",50,0,2); 

  Hlist->Add(h_proton_e);
  Hlist->Add(h_proton_ke);
  Hlist->Add(h_pi0_e);
  Hlist->Add(h_t);
  Hlist->Add(h_t_over_proton_ke);

  //
  // Loop over events
  //
  for(Long64_t i = 0; i < nEntries; i++) {

    // get next tree entry
    chain->GetEntry(i);

    if (i % 1000 == 0) {
      std::cout << "Event " << i << std::endl;
    }

    // get the GENIE event
    EventRecord &event = *(mcrec->event);
    Interaction *in = event.Summary();
    GHepParticle *nu = event.Probe();
    const ProcessInfo &proc = in->ProcInfo();

    // all the events are weak NC Coherent, DFR
    // we want proton energy for nu_e events

    int pdg = nu->Pdg();

    if (pdg == 14 || pdg == -14) continue;  // muon-type 
    if (pdg == 16 || pdg == -16) continue;  // tau-type 

    GHepParticle * p = 0;
    TIter event_iter(&event);

    h_t->Fill(in->KinePtr()->t(true));

    while ((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {
      if (p->Status() == kIStStableFinalState ) {
        int ppdg = p->Pdg();
        double pe = p->E();
        if (ppdg == kPdgPi0) {
          h_pi0_e->Fill(pe);
        }
        if (ppdg == kPdgProton) {
          h_proton_e->Fill(pe);
          h_proton_ke->Fill((pe - constants::kProtonMass) * 1000.0);
          h_t_over_proton_ke->Fill( in->KinePtr()->t(true) / 
              (pe - constants::kProtonMass) );
        }
        if (ppdg == kPdgAntiProton) {
          std::cout << "AntiProton E = " << p->E() << std::endl;
        }
      }
    }


    // clear current mc event record
    mcrec->Clear();

  } //end loop over events

  //
  // Close histogram file and clean up
  //
  std::string histname = "$ANA_HIST_DIR/dfr_gstana.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();

  delete h_proton_e;

  delete chain;

  return 0;
}
//_________________________________________________________________________________
