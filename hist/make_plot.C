make_plot() {
  gROOT->Reset();

  TCanvas *c1 = new TCanvas("c1","c1",2000,1200);
  TFile *f = new TFile("$ANA_HIST_DIR/dfr_gstana.root");
  TH1D *hw2 = (TH1D*)f->Get("h_dsigmadW2");

  c1->cd();
  c1->SetLogy();
  c1->SetLogx();
  hw2->Draw();
}
