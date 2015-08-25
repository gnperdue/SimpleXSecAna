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
  TH1D *h_proton_e = new TH1D("h_proton_e",
      "Proton Energy;GeV;Counts per 100 MeV",20,0,2); 

  Hlist->Add(h_proton_e);

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    h_proton_e->Fill(1.0);

  } // end event for loop

  std::string histname = "$ANA_HIST_DIR/dfr_gstana.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();

  delete h_proton_e;

} // end of Loop() function
