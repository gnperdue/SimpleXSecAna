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
    Long64_t max_events = 100000;
    std::string list_file_name = "default_nu_qe_like_carbon.txt";

    //
    // process command line options
    //
    int optind = 1;
    while ((optind < argc) && argv[optind][0] == '-') {
        std::string sw = argv[optind];
        if (sw == "-f" || sw == "--file") {
            optind++;
            list_file_name = argv[optind];
        }
        if (sw == "-m" || sw == "--max_events") {
            optind++;
            max_events = Long64_t(atoi(argv[optind]));
        }
        optind++;
    }

    //
    // Prepare and load the chain of GHEP files
    //
    char *s = getenv("ANA_LIST_DIR");
    std::string filelist = s + std::string("/") + list_file_name;
    std::cout << "List file = " << filelist << std::endl;
    std::ifstream filestream(filelist.c_str());
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

    int nbins = 80;
    double e_low = 0.0;
    double e_high = 20.0;
    double ebinwid = (e_high - e_low) / double(nbins);
    sprintf(axes1, "Cross Section;E_{#nu} (GeV);#frac{d#sigma}{dE} #times 10^{39} cm^{2} / GeV^{2}");
    TH1D *h_dsigmadE = new TH1D("h_dsigmadE", axes1, nbins, e_low, e_high); 
    TH1D *h_flux = new TH1D("h_flux", "E;n", nbins, e_low, e_high);

    Hlist->Add(h_dsigmadE);
    Hlist->Add(h_flux);

    //
    // Loop over events
    //
    for(Long64_t i = 0; i < nEntries; i++) {

        // get next tree entry
        chain->GetEntry(i);
        if (i % 20000 == 0) {
            std::cout << "Event " << i << std::endl;
        }

        // get the GENIE event
        EventRecord &event = *(mcrec->event);
        Interaction *in = event.Summary();
        GHepParticle *nu = event.Probe();
        const ProcessInfo &proc = in->ProcInfo();
        // const InitialState &init = in->InitState();

        // double nu_e = init.ProbeE();
        double nu_e = nu->Energy();
        int pdg = nu->Pdg();

        // is this a signal event?
        bool is_signal = false;
        if (proc.IsWeakCC() && proc.IsQuasiElastic()) {
            is_signal = true;
        }
        // loop over particles - define signal event, etc.
        /*
        GHepParticle * p = 0;
        TIter event_iter(&event);
        while ((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {
            if (p->Status() == kIStStableFinalState ) {
                int ppdg = p->Pdg();
                double pe = p->E();
                if (ppdg == kPdgProton) {
                    // blah, blah
                }
            }
        }
        */

        // get xsec in 1e-39 cm^2 and prep flux
        // (not a scalar n-entries normalization here) 
        if (is_signal) {
            double weight = 
                event.XSec() / (1E-39 * units::cm2) / ebinwid;
            h_dsigmadE->Fill(nu_e, weight);
            h_flux->Fill(nu_e, 1.0 / ebinwid);
        }

        // clear current mc event record
        mcrec->Clear();

        if (i > max_events) break;

    } //end loop over events

    //
    // divide by flux
    //
    h_dsigmadE->Divide(h_flux);

    //
    // Close histogram file and clean up
    //
    std::string histname = "$ANA_HIST_DIR/qetotalxsec_ghepana.root";
    TFile * outputfile = new TFile(histname.c_str(), "RECREATE");
    Hlist->Write();
    outputfile->Close();

    delete h_dsigmadE;
    delete h_flux;

    delete chain;

    return 0;
}
//_________________________________________________________________________________