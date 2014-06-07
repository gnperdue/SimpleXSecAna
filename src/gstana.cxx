#define gstana_cxx
#include "gstana.h"
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void gstana::Loop()
{
  if (fChain == 0) return;

  TObjArray * Hlist = new TObjArray(0);

  // Histograms
  TH1D *h_ev     = new TH1D("h_ev","Neutrino Energy;GeV;Counts per 500 MeV",20,0,10); 
  TH1D *h_neu    = new TH1D("h_neu","PDG Code;PDG Code;Counts per Code",40,-20,20); 
  TH1D *h_Q2s    = new TH1D("h_Q2s","Q^{2}s;GeV^{2};Counts per 50 bins",50,0,4); 
  TH1D *h_Q2     = new TH1D("h_Q2","Q^{2};GeV^{2};Counts per 50 bins",50,0,4); 

  Hlist->Add(h_ev);
  Hlist->Add(h_neu);
  Hlist->Add(h_Q2s);
  Hlist->Add(h_Q2);

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    h_ev->Fill(Ev);
    h_neu->Fill(neu);
    h_Q2s->Fill(Q2s);
    h_Q2->Fill(Q2);

  } // end event for loop

  std::string histname = "$ANA_HIST_DIR/qel_gstana.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();

  delete h_ev;
  delete h_neu;
  delete h_Q2s;
  delete h_Q2;

} // end of Loop() function
