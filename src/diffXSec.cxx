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
#include "Conventions/Units.h"

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

  int nbins = 32;
  double elow = 0.0;
  double ehigh = 3.0;
  double ebinwid = (ehigh - elow) / nbins;
  double ctlow = -1.0;
  double cthigh = 1.0;
  double ctbinwid = (cthigh - ctlow) / nbins;
  char axes1[200];
  sprintf(axes1, "Differential Cross Section;Muon Energy (GeV);XSec 1e38 cm^2 / %f GeV", ebinwid);
  TH1D *h_dsigmadT    = new TH1D("h_dsigmadT", axes1, nbins, elow, ehigh); 
  sprintf(axes1, "Differential Cross Section;Cosine Theta;XSec 1e38 cm^2 / %f Rad", ctbinwid);
  TH1D *h_dsigmadcosT = new TH1D("h_dsigmadcosT", axes1, nbins, ctlow, cthigh); 

  Hlist->Add(h_dsigmadT);
  Hlist->Add(h_dsigmadcosT);

  double total_xsec = 0.0;
  double gt_point8_fraction = 0.0;

  //
  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);
    Interaction *in = event.Summary();

    GHepParticle * fsl = event.FinalStatePrimaryLepton();
    double weight = event.XSec() / (1E-38 * units::cm2) / double(nev);
    total_xsec += weight;
    double lE = fsl->Energy();
    double lPx = fsl->Px();
    double lPy = fsl->Py();
    double lPz = fsl->Pz();
    double lP = sqrt( lPx * lPx + lPy * lPy + lPz * lPz );
    double cost = lPz / lP;
    if (cost > 0.8) {
      gt_point8_fraction += weight;
    }

 #if DEBUG
    LOG("myAnalysis", pNOTICE) << "---------------";
    LOG("myAnalysis", pNOTICE) << "total xsec = " << event.XSec() / (1E-38 * units::cm2);
    LOG("myAnalysis", pNOTICE) << "weight = " << weight;
    LOG("myAnalysis", pNOTICE) << "lE = " << lE;
    LOG("myAnalysis", pNOTICE) << "lPx = " << lPx;
    LOG("myAnalysis", pNOTICE) << "lPy = " << lPy;
    LOG("myAnalysis", pNOTICE) << "lPz = " << lPz;
    LOG("myAnalysis", pNOTICE) << "cost = " << cost;
#endif

    h_dsigmadT->Fill(lE, weight);
    h_dsigmadcosT->Fill(cost, weight);

    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  gt_point8_fraction /= total_xsec;

  std::string histname = "$ANA_HIST_DIR/diffXSec.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();

  LOG("myAnalysis", pNOTICE) << "total xsec = " << total_xsec;
  LOG("myAnalysis", pNOTICE) << "fraction at cos(t) > 0.8 = " << gt_point8_fraction;

  delete h_dsigmadT;
  delete h_dsigmadcosT;
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
