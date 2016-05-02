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
    int optind = 1;
    std::string filename = "gntp.10005.ghep.root";
    while ((optind < argc) && (argv[optind][0]=='-')) {
        std::string sw = argv[optind];
        if (sw=="-f") {
            optind++;
            filename = argv[optind];
        }
        optind++;
    }
    //
    // Prepare and load the chain of GHEP files
    //
    char *s = getenv( "ANA_LIST_DIR" );
    std::string filelist = s + std::string("/") + filename;

    TChain* chain = new TChain("gtree");
    std::cout << "  Adding " << filename << " to chain..." << std::endl;
    chain->Add(filename.c_str());

    Long64_t nEntries = chain->GetEntries(); 
    std::cout << "There are " << nEntries << " in the chain." << std::endl;

    NtpMCEventRecord * mcrec = 0;
    chain->SetBranchAddress("gmcrec", &mcrec);

    //
    // Set up Histograms
    //
    TObjArray * Hlist = new TObjArray(0);
    char axes1[200];

    TH1D *h_pion_e_15deg = new TH1D("h_pion_e_15deg",
            "Pion Energy at 15 degrees;GeV;Counts per 25 MeV",50,0.0,1.5); 
    TH1D *h_pion_e_30deg = new TH1D("h_pion_e_30deg",
            "Pion Energy at 15 degrees;GeV;Counts per 25 MeV",50,0.0,1.5); 

    int nbins = 50;
    double ke_low = 0.0;
    double ke_high = 1.5;
    double kebinwid = (ke_high - ke_low) / nbins;
    sprintf(axes1, "Differential Cross Section;Pion KE (GeV);#frac{d#sigma}{dT_{#pi}} #times 10^{40} cm^{2} / GeV");
    TH1D *h_dsigmadTpi = new TH1D("h_disgmaTpi", axes1, nbins, ke_low, ke_high);
    sprintf(axes1, "Differential XSec at 15 deg;Pion KE (GeV);#frac{d#sigma}{dT_{#pi}} #times 10^{40} cm^{2} / GeV");
    TH1D *h_dsigmadTpi_15deg = new TH1D("h_disgmaTpi_15deg", axes1, nbins, ke_low, ke_high);
    sprintf(axes1, "Differential XSec at 30 deg;Pion KE (GeV);#frac{d#sigma}{dT_{#pi}} #times 10^{40} cm^{2} / GeV");
    TH1D *h_dsigmadTpi_30deg = new TH1D("h_disgmaTpi_30deg", axes1, nbins, ke_low, ke_high);

    Hlist->Add(h_dsigmadTpi);
    Hlist->Add(h_dsigmadTpi_15deg);
    Hlist->Add(h_dsigmadTpi_30deg);

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
        // all the events in the files are weak Coherent (some CC, some NC)
        EventRecord &event = *(mcrec->event);
        Interaction *in = event.Summary();
        GHepParticle *nu = event.Probe();
        const ProcessInfo &proc = in->ProcInfo();

        GHepParticle * p = 0;
        TIter event_iter(&event);

        while ((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {
            if (p->Status() == kIStStableFinalState ) {
                int ppdg = p->Pdg();
                double ke = p->KinE();
                double angle = p->GetP4()->Theta()/units::degree;
                // we're only looking at coherent events, so we can be sloppy
                if (ppdg == kPdgPiP || ppdg == kPdgPi0 || ppdg == kPdgPiM) {
                    // get xsec in 1e-40 cm^2 and divide by flux
                    double weight = 
                        event.XSec() / (1E-40 * units::cm2) / double(nEntries) / kebinwid;
                    h_dsigmadTpi->Fill(ke, weight);
                    if (angle > 14.0 && angle < 16.0) {
                        h_dsigmadTpi_15deg->Fill(ke, weight);
                    }
                    if (angle > 29.0 && angle < 31.0) {
                        h_dsigmadTpi_15deg->Fill(ke, weight);
                    }
                }
            }
        }

        // clear current mc event record
        mcrec->Clear();

    } //end loop over events

    //
    // Close histogram file and clean up
    //
    std::string histname = "$ANA_HIST_DIR/coh_valid.root";
    TFile * outputfile = new TFile(histname.c_str(),"RECREATE");
    Hlist->Write();
    outputfile->Close();

    delete h_dsigmadTpi;
    delete h_dsigmadTpi_15deg;
    delete h_dsigmadTpi_30deg;

    delete chain;

    return 0;
}
//_________________________________________________________________________________
