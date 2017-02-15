//____________________________________________________________________________
#include <string>
#include <iostream>
#include <fstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TIterator.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

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




//____________________________________________________________________________
double get_bin_width(
        double value, const double bins[], const unsigned int nbins
        ) {
    if (nbins < 2) return -1.0;
    unsigned int index_check = 1;
    while (index_check <= nbins) {
        if (value < bins[index_check]) {
            return bins[index_check] - bins[index_check - 1];
        }
        index_check++;
    }
    return -1.0;
}

//____________________________________________________________________________
// check signal pdg
bool is_cc_numu_numubar(const EventRecord &event, int signal_pdg) {
    const ProcessInfo &proc = event.Summary()->ProcInfo();
    int pdg = event.Probe()->Pdg();
    return (pdg == signal_pdg) && proc.IsWeakCC();
}

//____________________________________________________________________________
// check if the event is ccqe-true == ccqe or mec
bool is_ccqe_true(const EventRecord &event, int signal_pdg) {
    const ProcessInfo &proc = event.Summary()->ProcInfo();
    if (is_cc_numu_numubar(event, signal_pdg)) {
        return proc.IsQuasiElastic() || proc.IsMEC();
    }
    return false;
}

//____________________________________________________________________________
// check if the event is ccqe-like
bool is_ccqe_like(
        const EventRecord &event, int signal_pdg, const double Tk_cut = 120.0
        ) {
    if (is_cc_numu_numubar(event, signal_pdg)) {

        int charged_muon_count = 0;
        int nu_pdg = event.Probe()->Pdg();

        GHepParticle * p = 0;
        TIter event_iter(&event);
        while ((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {
            if (p->Status() == kIStStableFinalState ) {
                int ppdg = p->Pdg();
                double pe = p->E();
                if (ppdg == kPdgProton) {
                    if (nu_pdg == kPdgAntiNuMu) {
                        double pe = p->E();
                        if (pe > Tk_cut) return false;
                    }
                }
                // absolutely no pions in the final state
                if (ppdg == kPdgPiP || ppdg == kPdgPiM || ppdg == kPdgPi0) {
                    return false;
                }
                // one and only one chaged "muon" in the final state
                if (ppdg == kPdgMuon || ppdg == kPdgAntiMuon) {
                    charged_muon_count++;
                    if (charged_muon_count > 1) {
                        return false;
                    }
                }
            }
        }
        return true;

    }
    return false;
}

//____________________________________________________________________________
// muon angle cut (rotated into beam coordinates)
bool angle_cut(const EventRecord &event, double angle_cut = 20.0) {
    angle_cut *= (constants::kPi / 180.);
    const TLorentzVector &p4fsl = *(event.FinalStatePrimaryLepton()->P4());
    TVector3 beam(0.0, 0.0, 1.0);
    // here, in pure generator case, angle is 0
    // beam.RotateX(3.0 * constants::kPi / 180.0);
    double angle = p4fsl.Angle(beam);
    return angle < angle_cut;
}

//___________________________________________________________________
int main(int argc, char ** argv)
{
    Long64_t max_events = 100000;
    std::string list_file_name = "default_nubar_qe_like_scint.txt";
    std::string out_file_name = "default_nubar_qe_like_scint.root";
    int signal_pdg = -14;

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
        if (sw == "-s" || sw == "--signal_pdg") {
            optind++;
            signal_pdg = atoi(argv[optind]);
        }
        if (sw == "-o" || sw == "--out_file") {
            optind++;
            out_file_name = argv[optind];
        }
        optind++;
    }

    if (signal_pdg != -14 && signal_pdg != 14) {
        std::cout << "This code is designed for muon neutrinos and muon";
        std::cout << "antineutrinos only. Please check your signal argument."
            << std::endl;
        return 0;
    }

    const double Ebinding = 30;

    //
    // set up histogram bins - this is ugly, but, okay, whatever
    //
    double enu_bins[] = {1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
    double q2_bins[] = {0, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 2.0};
    double nubar_pt_bins[] = {0, 0.15, 0.25, 0.4, 0.7, 1.0, 1.5};
    double nu_pt_bins[] = {0, 0.15, 0.25, 0.4, 0.7, 1.0, 1.5, 2.5, 4.0};
    double nubar_pl_bins[] = { 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0};
    double nu_pl_bins[] = {0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0};

    int nptelems = sizeof(nubar_pt_bins)/sizeof(nubar_pt_bins[0]);
    if (signal_pdg == 14) {
        nptelems = sizeof(nu_pt_bins)/sizeof(nu_pt_bins[0]);
    }
    double* pt_bins = new double[nptelems]();
    if (signal_pdg == -14) {
        for (int i = 0; i < nptelems; ++i)
            pt_bins[i] = nubar_pt_bins[i];
    }
    if (signal_pdg == 14) {
        for (int i = 0; i < nptelems; ++i)
            pt_bins[i] = nu_pt_bins[i];
    }

    int nplelems = sizeof(nubar_pl_bins)/sizeof(nubar_pl_bins[0]);
    if (signal_pdg == 14) {
        nplelems = sizeof(nu_pl_bins)/sizeof(nu_pl_bins[0]);
    }
    double* pl_bins = new double[nplelems]();
    if (signal_pdg == -14) {
        for (int i = 0; i < nplelems; ++i)
            pl_bins[i] = nubar_pl_bins[i];
    }
    if (signal_pdg == 14) {
        for (int i = 0; i < nplelems; ++i)
            pl_bins[i] = nu_pl_bins[i];
    }

    const unsigned int enu_nbins = sizeof(enu_bins) / sizeof(enu_bins[0]) - 1;
    const unsigned int q2_nbins = sizeof(q2_bins) / sizeof(q2_bins[0]) - 1;
    const unsigned int pt_nbins = sizeof(pt_bins) / sizeof(pt_bins[0]) - 1;
    const unsigned int pl_nbins = sizeof(pl_bins) / sizeof(pl_bins[0]) - 1;

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
    Long64_t nev = (max_events > 0) ?
        TMath::Min(max_events, nEntries) : nEntries;
    std::cout << "There are " << nEntries << " in the chain." << std::endl;
    std::cout << "We will process " << nev << " events." << std::endl;

    NtpMCEventRecord * mcrec = 0;
    chain->SetBranchAddress("gmcrec", &mcrec);

    //
    // Set up Histograms
    //
    TObjArray * Hlist = new TObjArray(0);

    char axes1[200];

    sprintf(axes1, "Differential cross section - true CCQE or MEC;Q^{2} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_true_nocut = new TH1D("q2_true_nocut", axes1, q2_nbins, q2_bins);

    sprintf(axes1, "Differential cross section - true CCQE or MEC (angle cut);Q^{2} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_true_cut = new TH1D("q2_true_cut", axes1, q2_nbins, q2_bins);

    sprintf(axes1, "Differential cross section - QE-like;Q^{2} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_like_nocut = new TH1D("q2_like_nocut", axes1, q2_nbins, q2_bins);

    sprintf(axes1, "Differential cross section - QE-like (angle cut);Q^{2} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_like_cut = new TH1D("q2_like_cut", axes1, q2_nbins, q2_bins);

    Hlist->Add(q2_true_nocut);
    Hlist->Add(q2_true_cut);
    Hlist->Add(q2_like_nocut);
    Hlist->Add(q2_like_cut);

    //
    // Loop over events
    //
    for(Long64_t i = 0; i < nev; i++) {

        // get next tree entry
        chain->GetEntry(i);
        if (i % 20000 == 0) {
            std::cout << "Event " << i << std::endl;
        }

        // get the GENIE event
        EventRecord &event = *(mcrec->event);
        const Kinematics &kine = event.Summary()->Kine();

        double ccqe_true_total_xsec = 0.0;
        double ccqe_like_total_xsec = 0.0;

        // is this a signal event?
        bool is_ccqe_true_event = is_ccqe_true(event, signal_pdg);
        bool is_ccqe_like_event = is_ccqe_like(event, signal_pdg);

        // get xsec in 1e-39 cm^2 
        if (is_ccqe_true_event || is_ccqe_like_event) {
            // TODO: get Q2 from lepton vars
            double Q2 = kine.Q2(true);   // think we want selected == true
            double q2_bin_wid = get_bin_width(Q2, q2_bins, q2_nbins);
            double q2_weight = 
                event.XSec() / (1E-39 * units::cm2) / q2_bin_wid / double(nev);

            if (is_ccqe_true_event) {
                ccqe_true_total_xsec += q2_weight;
                q2_true_nocut->Fill(Q2, q2_weight);
                if (angle_cut(event)) {
                    q2_true_cut->Fill(Q2, q2_weight);
                }
            }

            if (is_ccqe_like_event) {
                ccqe_like_total_xsec += q2_weight;
                q2_like_nocut->Fill(Q2, q2_weight);
                if (angle_cut(event)) {
                    q2_like_cut->Fill(Q2, q2_weight);
                }
            }
        }

        // clear current mc event record
        mcrec->Clear();

    } //end loop over events

    //
    // division by flux handled in event weight (flux integrated)
    //

    //
    // Close histogram file and clean up
    //
    std::string histname = "$ANA_HIST_DIR/" + out_file_name;
    TFile * outputfile = new TFile(histname.c_str(), "RECREATE");
    Hlist->Write();
    outputfile->Close();

    delete q2_true_nocut;
    delete q2_true_cut;
    delete q2_like_nocut;
    delete q2_like_cut;

    delete[] pt_bins;
    delete[] pl_bins;

    delete chain;

    return 0;
}

/*
// filename - input file w/o .root; assuming in root/
void nuwro_ccqe_extractor(string filename)
{

    TH2D* enuq2_true_nocut = new TH2D("enuq2_true_nocut",
            "d^2#sigma/denudq2 (true)",
            q2_nbins, q2_bins, enu_nbins, enu_bins);

    TH2D* enuq2_true_cut = new TH2D("enuq2_true_cut",
            "d^2#sigma/denudq2 (true+cut)",
            q2_nbins, q2_bins, enu_nbins, enu_bins);

    TH2D* enuq2_like_nocut = new TH2D("enuq2_like_nocut",
            "d^2#sigma/denudq2 (like)",
            q2_nbins, q2_bins, enu_nbins, enu_bins);

    TH2D* enuq2_like_cut = new TH2D("enuq2_like_cut",
            "d^2#sigma/denudq2 (like+cut)",
            q2_nbins, q2_bins, enu_nbins, enu_bins);

    TH2D* ptpl_true_nocut = new TH2D("ptpl_true_nocut",
            "d^2#sigma/dpTdpL (true)",
            pl_nbins, pl_bins, pt_nbins, pt_bins);

    TH2D* ptpl_true_cut = new TH2D("ptpl_true_cut",
            "d^2#sigma/dpTdpL (true+cut)",
            pl_nbins, pl_bins, pt_nbins, pt_bins);

    TH2D* ptpl_like_nocut = new TH2D("ptpl_like_nocut",
            "d^2#sigma/dpTdpL (like)",
            pl_nbins, pl_bins, pt_nbins, pt_bins);

    TH2D* ptpl_like_cut = new TH2D("ptpl_like_cut",
            "d^2#sigma/dpTdpL (like+cut)",
            pl_nbins, pl_bins, pt_nbins, pt_bins);

    for (unsigned int i = 0; i < nof_events; i++)
    {
        tree->GetEntry(i);

        const double Q2 = -e->q2() / 1000000.0; // in GeV^2
        const double pl = e->out[0].p().z / 1000.0; // in GeV
        const double pt = sqrt(e->out[0].p().x * e->out[0].p().x
                + e->out[0].p().y * e->out[0].p().y) / 1000.0;
        const double p =  sqrt(e->out[0].p().x * e->out[0].p().x
                + e->out[0].p().y * e->out[0].p().y
                + e->out[0].p().z * e->out[0].p().z) ;
        const double E =  e->out[0].E();
        const double theta = acos(pl*1000.0/p);//pl was in GeV
        const double enu = e->in[0].E()/1000.0; //in GeV

        const double enuQE = (pow(constants::kProtonMass,2)-pow((constants::kNeutronMass-Ebinding),2)-pow(constants::kMuonMass,2)+2*(constants::kNeutronMass-Ebinding)*E)/(2*(constants::kNeutronMass-Ebinding-E+p*cos(theta)));
        const double Q2QE = 2*enuQE*(E-p*cos(theta))-pow(constants::kMuonMass,2);
        //	cout << enuQE << "\t" << Q2QE <<"\t" << p << "\t" << E << "\t" << theta << endl;	

        if (is_ccqe_true(e))
        {
            q2_true_nocut->Fill(Q2);
            ptpl_true_nocut->Fill(pl, pt);
            enuq2_true_nocut->Fill(enuQE, Q2QE);

            if (angle_cut(e))
            {
                q2_true_cut->Fill(Q2);
                ptpl_true_cut->Fill(pl, pt);
                enuq2_true_cut->Fill(enuQE, Q2QE);
            }
        }

        if (is_ccqe_like(e))
        {
            q2_like_nocut->Fill(Q2);
            ptpl_like_nocut->Fill(pl, pt);
            enuq2_like_nocut->Fill(enuQE, Q2QE);

            if (angle_cut(e))
            {
                q2_like_cut->Fill(Q2);
                ptpl_like_cut->Fill(pl, pt); 
                enuq2_like_cut->Fill(enuQE, Q2QE);
            }
        }

        if (100 * i % nof_events == 0)
            cout << 100 * i / nof_events << "%     \r" << flush;
    }

    normalize (q2_true_nocut, xsec / nof_events);
    normalize (q2_true_cut, xsec / nof_events);
    normalize (q2_like_nocut, xsec / nof_events);
    normalize (q2_like_cut, xsec / nof_events);

    normalize2D (ptpl_true_nocut, xsec / nof_events);
    normalize2D (ptpl_true_cut, xsec / nof_events);
    normalize2D (ptpl_like_nocut, xsec / nof_events);
    normalize2D (ptpl_like_cut, xsec / nof_events);

    normalize2D (enuq2_true_nocut, xsec / nof_events);
    normalize2D (enuq2_true_cut, xsec / nof_events);
    normalize2D (enuq2_like_nocut, xsec / nof_events);
    normalize2D (enuq2_like_cut, xsec / nof_events);

    output->Write();
}
*/
