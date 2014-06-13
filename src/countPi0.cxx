#include <map>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>
#include <TH2.h>
#include <TGraph.h>

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

  double e[nev], xs[nev];
  std::map<double, double> exs;

  TObjArray * Hlist = new TObjArray(0);

  TH1D *h_npi0       = new TH1D("h_npi0","Number of Neutral Pions;Number;Counts",10,0,10); 
  TH1D *h_npipi0only = new TH1D("h_npipi0only","Number of Neutral Pions (no Charged Pions);Number;Counts",10,0,10); 
  TH1D *h_dsigdE     = new TH1D("h_dsigdE","XSec(E);Energy",100,1,3);

  Hlist->Add(h_npi0);
  Hlist->Add(h_npipi0only);
  Hlist->Add(h_dsigdE);

  double total_xsec = 0;
  double total_ccpi0_xsec = 0;
  double smallest_xsec = 1000.0;  // in units of 1e-38 cm^2
  double largest_xsec = 0.0;  // in units of 1e-38 cm^2

  //
  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);
    Interaction *in = event.Summary();
    GHepParticle * nu = event.Probe();
    const ProcessInfo &proc = in->ProcInfo();

    // check to see if the event is charged current
    bool iscc = proc.IsWeakCC();

    // get the energy and cross section in units of 1e-38 cm^2
    double nue = nu->Energy();
    double nxsec = event.XSec() / (1E-38 * units::cm2);
    // save the energy and cross section in a map keyed by energy
    exs.insert(std::pair<double,double>(nue, nxsec));
    // compute the total cross section by adding events and dividing by "flux"
    double weight = nxsec / double(nev);
    total_xsec += weight;
    // fill a histogram of energy with events weighted by cross section / flux
    h_dsigdE->Fill(nue, weight);
#if DEBUG
    LOG("myAnalysis", pNOTICE)  << "neutrino e = " << nue;
    LOG("myAnalysis", pNOTICE)  << "      xsec = " << nxsec;
    LOG("myAnalysis", pNOTICE)  << "   xsec/ev = " << weight;
#endif
    if (nxsec < smallest_xsec) smallest_xsec = nxsec;
    if (nxsec > largest_xsec) largest_xsec = nxsec;

    //
    // Loop over all particles in this event
    //
    GHepParticle * p = 0;
    TIter event_iter(&event);

    // count pions
    int npi0 = 0;
    int nchgdpi = 0;
    bool havePi0 = false;

    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
      if (p->Status() == kIStStableFinalState ) {
        if (p->Pdg() == kPdgPi0) {
          npi0++;
          if (iscc) havePi0 = true;
        }
        if (p->Pdg() == kPdgPiP ||
            p->Pdg() == kPdgPiM) {
          nchgdpi++;
        }
      }

    }// end loop over particles	

    if (havePi0) total_ccpi0_xsec += weight;

    h_npi0->Fill(npi0);
    if (nchgdpi==0) {
      h_npipi0only->Fill(npi0);
    }

    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  int indx = 0;
  for (std::map<double,double>::iterator it=exs.begin(); it!=exs.end(); ++it) {
    e[indx] = it->first;
    xs[indx] = it->second;
    ++indx;
  }

  TGraph* exs_gr = new TGraph(nev, e, xs);
  Hlist->Add(exs_gr);

  LOG("myAnalysis", pNOTICE) << "total xsec = " << total_xsec;
  LOG("myAnalysis", pNOTICE) << "total xsec by integral = " << h_dsigdE->Integral();
  LOG("myAnalysis", pNOTICE) << "smallest observed xsec = " << smallest_xsec;
  LOG("myAnalysis", pNOTICE) << "largest observed xsec = " << largest_xsec;
  LOG("myAnalysis", pNOTICE) << "total ccpi0 = " << total_ccpi0_xsec;

  // save histograms
  std::string histname = "$ANA_HIST_DIR/countPi0.root";
  TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
  Hlist->Write();
  outputfile->Close();

  // clean up memory
  delete h_npi0;
  delete h_npipi0only;
  delete h_dsigdE;
  delete exs_gr;
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
