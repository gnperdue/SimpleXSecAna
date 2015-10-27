//____________________________________________________________________________
#include <string>
#include <fstream>
#include <iostream>

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
    std::string filelist = s + std::string("/ghep_file_list.txt");

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
    char axes1[200];

    int nbins = 100;
    double hist_low = 0.5;
    double hist_high = 20.5;
    double wbinwid = (hist_high - hist_low) / nbins;
    sprintf(axes1, "#nu Differential Cross Section;W^{2} (GeV^{2});#frac{d#sigma}{dW^{2}} #times 10^{40} cm^{2} / GeV^{2}");
    TH1D *h_nu_dsigmadW2 = new TH1D("h_nu_dsigmadW2", axes1, nbins, hist_low, hist_high); 
    sprintf(axes1, "#bar{#nu} Differential Cross Section;W^{2} (GeV^{2});#frac{d#sigma}{dW^{2}} #times 10^{40} cm^{2} / GeV^{2}");
    TH1D *h_antinu_dsigmadW2 = new TH1D("h_antinu_dsigmadW2", axes1, nbins, hist_low, hist_high); 

    nbins = 100;
    hist_low = 0.0;
    hist_high = 10.0;
    double q2binwid = (hist_high - hist_low) / nbins;
    sprintf(axes1, "#nu Differential Cross Section;Q^{2} (GeV^{2});#frac{d#sigma}{dQ^{2}} #times 10^{40} cm^{2} / GeV^{2}");
    TH1D *h_nu_dsigmadQ2 = new TH1D("h_nu_dsigmadQ2", axes1, nbins, hist_low, hist_high); 
    sprintf(axes1, "#bar{#nu} Differential Cross Section;Q^{2} (GeV^{2});#frac{d#sigma}{dQ^{2}} #times 10^{40} cm^{2} / GeV^{2}");
    TH1D *h_antinu_dsigmadQ2 = new TH1D("h_antinu_dsigmadQ2", axes1, nbins, hist_low, hist_high); 

    nbins = 100;
    hist_low = 0.0;
    hist_high = 1.0;
    double xbinwid = (hist_high - hist_low) / nbins;
    sprintf(axes1, "#nu Differential Cross Section;x_{Bj};#frac{d#sigma}{dx} #times 10^{40} cm^{2}");
    TH1D *h_nu_dsigmadx = new TH1D("h_nu_dsigmadx", axes1, nbins, hist_low, hist_high); 
    TH1D *h_antinu_dsigmadx = new TH1D("h_antinu_dsigmadx", axes1, nbins, hist_low, hist_high); 

    nbins = 10;
    hist_low = 0.0;
    hist_high = 10.0;
    double npibinwid = (hist_high - hist_low) / nbins;
    sprintf(axes1, "#nu Differential Cross Section;N_{#pi};#frac{d#sigma}{dN_{#pi}} #times 10^{40} cm^{2}");
    TH1D *h_nu_dsigmadnpi = new TH1D("h_nu_dsigmadnpi", axes1, nbins, hist_low, hist_high); 
    sprintf(axes1, "#bar{#nu} Differential Cross Section;N_{#pi};#frac{d#sigma}{dN_{#pi}} #times 10^{40} cm^{2}");
    TH1D *h_antinu_dsigmadnpi = new TH1D("h_antinu_dsigmadnpi", axes1, nbins, hist_low, hist_high); 

    nbins = 10;
    hist_low = 0.0;
    hist_high = 10.0;
    double npi0binwid = (hist_high - hist_low) / nbins;
    sprintf(axes1, "#nu Differential Cross Section;N_{#pi^{0}};#frac{d#sigma}{dN_{#pi^{0}}} #times 10^{40} cm^{2}");
    TH1D *h_nu_dsigmadnpi0 = new TH1D("h_nu_dsigmadnpi0", axes1, nbins, hist_low, hist_high); 
    sprintf(axes1, "#bar{#nu} Differential Cross Section;N_{#pi^{0}};#frac{d#sigma}{dN_{#pi^{0}}} #times 10^{40} cm^{2}");
    TH1D *h_antinu_dsigmadnpi0 = new TH1D("h_antinu_dsigmadnpi0", axes1, nbins, hist_low, hist_high); 

    Hlist->Add(h_nu_dsigmadW2);
    Hlist->Add(h_nu_dsigmadQ2);
    Hlist->Add(h_nu_dsigmadx);
    Hlist->Add(h_nu_dsigmadnpi);
    Hlist->Add(h_nu_dsigmadnpi0);

    Hlist->Add(h_antinu_dsigmadW2);
    Hlist->Add(h_antinu_dsigmadQ2);
    Hlist->Add(h_antinu_dsigmadx);
    Hlist->Add(h_antinu_dsigmadnpi);
    Hlist->Add(h_antinu_dsigmadnpi0);

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
        // all the events in the files are weak NC Coherent, DFR
        EventRecord &event = *(mcrec->event);
        Interaction *in = event.Summary();
        GHepParticle *nu = event.Probe();
        const ProcessInfo &proc = in->ProcInfo();

        int pdg = nu->Pdg();
        // std::cout << pdg << std::endl;
        bool is_antinu = pdg < 0 ? true : false;
        double nue = nu->Energy();

        if (nue < 2.0) continue;
        if (nue > 6.0) continue;

        GHepParticle * p = 0;
        TIter event_iter(&event);

        // get xsec in 1e-40 cm^2 and divide by flux
        double weight = 
            event.XSec() / (1E-40 * units::cm2) / double(nEntries) / wbinwid;
        if (is_antinu) {
            h_antinu_dsigmadW2->Fill(in->KinePtr()->W(true)*in->KinePtr()->W(true), 
                    weight);
        }
        else {
            h_nu_dsigmadW2->Fill(in->KinePtr()->W(true)*in->KinePtr()->W(true), 
                    weight);
        }

        weight = event.XSec() / (1E-40 * units::cm2) / double(nEntries) / q2binwid;
        if (is_antinu) {
            h_antinu_dsigmadQ2->Fill(in->KinePtr()->Q2(true), weight);
        }
        else {
            h_nu_dsigmadQ2->Fill(in->KinePtr()->Q2(true), weight);
        }

        weight = event.XSec() / (1E-40 * units::cm2) / double(nEntries) / xbinwid;
        if (is_antinu) {
            h_antinu_dsigmadx->Fill(in->KinePtr()->x(true), weight);
        }
        else {
            h_nu_dsigmadx->Fill(in->KinePtr()->x(true), weight);
        }

        int npi0 = 0;
        int nchgdpi = 0;
        while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {
            if (p->Status() == kIStStableFinalState ) {
                if (p->Pdg() == kPdgPi0) {
                    npi0++;
                }
                if (p->Pdg() == kPdgPiP || p->Pdg() == kPdgPiM) {
                    nchgdpi++;
                }
            }
        } // end loop over particles   

        weight = event.XSec() / (1E-40 * units::cm2) / double(nEntries) / npibinwid;
        if (is_antinu) {
            h_antinu_dsigmadnpi->Fill(nchgdpi, weight);
        }
        else {
            h_nu_dsigmadnpi->Fill(nchgdpi, weight);
        }

        weight = event.XSec() / (1E-40 * units::cm2) / double(nEntries) / npi0binwid;
        if (is_antinu) {
            h_antinu_dsigmadnpi0->Fill(npi0, weight);
        }
        else {
            h_nu_dsigmadnpi0->Fill(npi0, weight);
        }

        // clear current mc event record
        mcrec->Clear();

    } //end loop over events

    //
    // Close histogram file and clean up
    //
    std::string histname = "$ANA_HIST_DIR/trans_gstana.root";
    TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
    Hlist->Write();
    outputfile->Close();

    delete chain;

    return 0;
}
//_________________________________________________________________________________
