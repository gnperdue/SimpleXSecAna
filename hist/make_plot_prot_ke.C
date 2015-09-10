make_plot_prot_ke() {
  gROOT->Reset();

  TCanvas *c1 = new TCanvas("c1","c1",2000,1200);
  TFile *f = new TFile("$ANA_HIST_DIR/dfr_gstana.root");
  TH1D *hpke = (TH1D*)f->Get("h_proton_ke");

  c1->cd();
  hpke->Draw();

  // set stat and fit data reported
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1110);

  hpke->Fit("expo","","",50,350);

  // set histo and fit line widths
  hpke->SetLineWidth(6);
  TF1 *efunc = hpke->GetFunction("expo");
  efunc->SetLineWidth(6);


  // Update title and axis title
  hpke->SetTitle("Proton Kinetic Energy");
  hpke->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  hpke->GetXaxis()->SetNoExponent(kTRUE);
  hpke->GetYaxis()->SetNoExponent(kTRUE);

  // re-position stats box
  TPaveStats *st = (TPaveStats*)hpke->FindObject("stats");
  // std::cout << st->GetX1NDC() << std::endl; 
  // std::cout << st->GetX2NDC() << std::endl; 
  st->SetX1NDC(0.60); //new x start position
  st->SetX2NDC(0.85); //new x end position
  st->SetY1NDC(0.60); //new y start position
  st->SetY2NDC(0.85); //new y end position

  // save a PDF
  c1->Print("GENIE_diffractive_model_NC_proton_ke.pdf");
}
