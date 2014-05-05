#define grooana_cxx
#include "grooana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void grooana::Loop()
{
  if (fChain == 0) return;

  TObjArray * Hlist = new TObjArray(0);

  // Histograms
  TH1D *h_evtXSec      = new TH1D("h_evtXSec","Cross Section;Cross Section;Counts per 50 bins",50,0,0.8); 
  TH1D *h_evtDXSec     = new TH1D("h_evtDXSec","Differential Cross Section;Differential Cross Section;Counts per 50 bins",50,0,250); 

  Hlist->Add(h_evtXSec);
  Hlist->Add(h_evtDXSec);

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    h_evtXSec->Fill(EvtXSec);
    h_evtDXSec->Fill(EvtDXSec);
  }

  std::string histname = "$ANA_HIST_DIR/qel_grooana.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();

  delete h_evtXSec;
  delete h_evtDXSec;
}
