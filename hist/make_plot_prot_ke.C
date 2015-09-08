make_plot_prot_ke() {
  gROOT->Reset();

  TCanvas *c1 = new TCanvas("c1","c1",2000,1200);
  TFile *f = new TFile("$ANA_HIST_DIR/dfr_gstana.root");
  TH1D *hpke = (TH1D*)f->Get("h_proton_ke");

  c1->cd();
  hpke->Draw();
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0100);
  hpke->Fit("expo","","",5,40);
}
