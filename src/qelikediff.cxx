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
// muon angle
double get_angle(const EventRecord &event) {
    const TLorentzVector &p4fsl = *(event.FinalStatePrimaryLepton()->P4());
    TVector3 beam(0.0, 0.0, 1.0);
    // here, in pure generator case, angle is 0
    // beam.RotateX(3.0 * constants::kPi / 180.0);
    double angle = p4fsl.Angle(beam);
    return angle;
}

//____________________________________________________________________________
// muon angle cut (rotated into beam coordinates)
bool angle_cut(const EventRecord &event, double angle_cut = 20.0) {
    angle_cut *= (constants::kPi / 180.);
    double angle = get_angle(event);
    return angle < angle_cut;
}

//____________________________________________________________________________
// normalize to cross section per nucleon
void normalize(TH1* h, const TH1* flux, const double factor) {
    h->Sumw2();
    h->Divide(flux);
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        h->SetBinContent(i, h->GetBinContent(i) * factor);
        // h->SetBinError(i, h->GetBinError(i) * factor);
        h->SetBinError(i, h->GetBinError(i) * factor);
    }
}

//____________________________________________________________________________
// normalize to cross section per nucleon
void normalize(TH1* h, const double factor) {
    h->Sumw2();
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        h->SetBinContent(i, h->GetBinContent(i) * factor);
        // h->SetBinError(i, h->GetBinError(i) * factor);
        h->SetBinError(i, 0.0);
    }
}

//____________________________________________________________________________
// normalize to cross section per nucleon
void normalize(TH1* h, const double factor,
        const double bins[], const unsigned int nbins) {
    // note - nbins is not the size of the `bins` array - it is size-1
    /////////////
    if (nbins < 2) return;
    double bin_divisors[nbins];
    for (unsigned int i = 1; i < (nbins + 1); ++i) {
        double divisor = bins[i] - bins[i - 1];
        bin_divisors[i - 1] = divisor;
    }
    h->Sumw2();
    int n_histo_bins = h->GetNbinsX();
    if (n_histo_bins != int(nbins)) return;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double bin_cont = h->GetBinContent(i);
        bin_cont = bin_cont * factor;
        bin_cont = bin_cont / bin_divisors[i - 1];
        h->SetBinContent(i, bin_cont);
        // h->SetBinContent(i, h->GetBinContent(i) * factor / bin_divisors[i - 1]);
        // h->SetBinError(i, h->GetBinError(i) * factor);
        h->SetBinError(i, 0.0);
    }
}

// double get_bin_width(
//         double value, const double bins[], const unsigned int nbins
//         ) {
//     if (nbins < 2) return -1.0;
//     unsigned int index_check = 1;
//     while (index_check <= nbins) {
//         if (value < bins[index_check]) {
//             return bins[index_check] - bins[index_check - 1];
//         }
//         index_check++;
//     }
//     return -1.0;
// }

//____________________________________________________________________________
// normalize to cross section per nucleon
void normalize2D(TH2* h, const double factor) {
    h->Sumw2();
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            h->SetBinContent(i, j, h->GetBinContent(i, j) * factor);
            h->SetBinError(i, j, h->GetBinError(i, j) * factor);
        }
    }
}

//____________________________________________________________________________
int main(int argc, char ** argv)
{
    Long64_t max_events = 100000;
    std::string list_file_name = "default_nubar_qe_like_scint.txt";
    std::string out_file_name = "default_nubar_qe_like_scint.root";
    int signal_pdg = -14;
    int per_nucleon_correction_factor = 7;   // 7 protons in CH
    double flux_e_min = 1.5;
    double flux_e_max = 10.0;

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
        if (sw == "-c" || sw == "--nucleon_correction") {
            optind++;
            per_nucleon_correction_factor = atoi(argv[optind]);
        }
        if (sw == "-n" || sw == "--flux_e_min") {
            optind++;
            flux_e_min = atof(argv[optind]);
        }
        if (sw == "-x" || sw == "--flux_e_max") {
            optind++;
            flux_e_max = atof(argv[optind]);
        }
        optind++;
    }

    if (signal_pdg != -14 && signal_pdg != 14) {
        std::cout << "This code is designed for muon neutrinos and muon ";
        std::cout << "antineutrinos only. Please check your signal argument."
            << std::endl;
        return 0;
    }

    const double Ebinding = 0.030;

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
    const unsigned int pt_nbins = signal_pdg == 14 ? 
        sizeof(nu_pt_bins) / sizeof(nu_pt_bins[0]) - 1 :
        sizeof(nubar_pt_bins) / sizeof(nubar_pt_bins[0]) - 1;
    const unsigned int pl_nbins = signal_pdg == 14 ? 
        sizeof(nu_pl_bins) / sizeof(nu_pl_bins[0]) - 1 :
        sizeof(nubar_pl_bins) / sizeof(nubar_pl_bins[0]) - 1;

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

    int nbins = 100;
    double e_low = 0.0;
    double e_high = 50.0;
    double ebinwid = (e_high - e_low) / double(nbins);
    sprintf(axes1, "Cross Section;E_{#nu} (GeV);#frac{d#sigma}{dE} #times 10^{39} cm^{2} / GeV^{2}");
    TH1D *h_dsigmadE = new TH1D("h_dsigmadE", axes1, nbins, e_low, e_high); 
    TH1D *h_flux = new TH1D("h_flux", "E;n", nbins, e_low, e_high);

    Hlist->Add(h_dsigmadE);
    Hlist->Add(h_flux);

    // CC Inclusive

    sprintf(axes1, "Differential cross section - CC;Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_cc = new TH1D("q2_cc", axes1, q2_nbins, q2_bins);

    sprintf(axes1, "Number of CC Events;Q^{2}_{QE} (GeV^{2});N-Events");
    TH1D* q2_n_cc_events = new TH1D("q2_n_cc_events", axes1, q2_nbins, q2_bins);

    // CCQE "True"

    sprintf(axes1, "Differential cross section - true CCQE or MEC;Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_true_nocut = new TH1D("q2_true_nocut", axes1, q2_nbins, q2_bins);

    sprintf(axes1, "Number of events - true CCQE or MEC;Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_n_true_nocut_events = new TH1D("q2_n_true_nocut_events", axes1, q2_nbins, q2_bins);

    // CCQE "True" - angle cut

    sprintf(axes1, "Differential cross section - true CCQE or MEC (angle cut);Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_true_cut = new TH1D("q2_true_cut", axes1, q2_nbins, q2_bins);

    // CCQE-like

    sprintf(axes1, "Differential cross section - QE-like;Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_like_nocut = new TH1D("q2_like_nocut", axes1, q2_nbins, q2_bins);

    // CCQE-like - angle cut

    sprintf(axes1, "Differential cross section - QE-like (angle cut);Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_like_cut = new TH1D("q2_like_cut", axes1, q2_nbins, q2_bins);

    // 2D

    sprintf(axes1, "Differential cross section - true CCQE or MEC;Q^{2}_{QE} (GeV^{2});E_{#nu-QE} (GeV)");
    TH2D* enuq2_true_nocut = new TH2D("enuq2_true_nocut",
            axes1, q2_nbins, q2_bins, enu_nbins, enu_bins);

    sprintf(axes1, "Differential cross section - true CCQE or MEC (angle cut);Q^{2}_{QE} (GeV^{2});E_{#nu-QE} (GeV)");
    TH2D* enuq2_true_cut = new TH2D("enuq2_true_cut",
            axes1, q2_nbins, q2_bins, enu_nbins, enu_bins);

    sprintf(axes1, "Differential cross section - QE-like;Q^{2}_{QE} (GeV^{2});E_{#nu-QE} (GeV)");
    TH2D* enuq2_like_nocut = new TH2D("enuq2_like_nocut",
            axes1, q2_nbins, q2_bins, enu_nbins, enu_bins);

    sprintf(axes1, "Differential cross section - QE-like (angle cut);Q^{2}_{QE} (GeV^{2});E_{#nu-QE} (GeV)");
    TH2D* enuq2_like_cut = new TH2D("enuq2_like_cut",
            axes1, q2_nbins, q2_bins, enu_nbins, enu_bins);

    sprintf(axes1, "Differential cross section - true CCQE or MEC;P_{#mu-L} (GeV);P_{#mu-T} (GeV)");
    TH2D* ptpl_true_nocut = new TH2D("ptpl_true_nocut",
            axes1, pl_nbins, pl_bins, pt_nbins, pt_bins);

    sprintf(axes1, "Differential cross section - true CCQE or MEC (angle cut);P_{#mu-L} (GeV);P_{#mu-T} (GeV)");
    TH2D* ptpl_true_cut = new TH2D("ptpl_true_cut",
            axes1, pl_nbins, pl_bins, pt_nbins, pt_bins);

    sprintf(axes1, "Differential cross section - QE-like;P_{#mu-L} (GeV);P_{#mu-T} (GeV)");
    TH2D* ptpl_like_nocut = new TH2D("ptpl_like_nocut",
            axes1, pl_nbins, pl_bins, pt_nbins, pt_bins);

    sprintf(axes1, "Differential cross section - QE-like (angle cut);P_{#mu-L} (GeV);P_{#mu-T} (GeV)");
    TH2D* ptpl_like_cut = new TH2D("ptpl_like_cut",
            axes1, pl_nbins, pl_bins, pt_nbins, pt_bins);

    Hlist->Add(q2_cc);
    Hlist->Add(q2_n_cc_events);

    Hlist->Add(q2_true_nocut);
    Hlist->Add(q2_n_true_nocut_events);

    Hlist->Add(q2_true_cut);
    Hlist->Add(q2_like_nocut);
    Hlist->Add(q2_like_cut);

    Hlist->Add(enuq2_true_nocut);
    Hlist->Add(enuq2_true_cut);
    Hlist->Add(enuq2_like_nocut);
    Hlist->Add(enuq2_like_cut);

    Hlist->Add(ptpl_true_nocut);
    Hlist->Add(ptpl_true_cut);
    Hlist->Add(ptpl_like_nocut);
    Hlist->Add(ptpl_like_cut);

    //
    // Loop over events
    //
    int n_cc_events = 0;
    double cc_total_cross_section = 0.0;
    for(Long64_t i = 0; i < nev; i++) {

        // get next tree entry
        chain->GetEntry(i);
        if (i % 20000 == 0) {
            std::cout << "Event " << i << std::endl;
        }

        // get the GENIE event
        EventRecord &event = *(mcrec->event);
        // is this a signal event?
        bool is_cc_event = is_cc_numu_numubar(event, signal_pdg);
        bool is_ccqe_true_event = is_ccqe_true(event, signal_pdg);
        bool is_ccqe_like_event = is_ccqe_like(event, signal_pdg);

        if (!is_cc_event) {
            mcrec->Clear();
            continue;
        }

        GHepParticle *nu = event.Probe();
        double enu = nu->Energy();
        if ((enu < flux_e_min) || (enu > flux_e_max)) {
            mcrec->Clear();
            continue;
        } 

        // for pure GENIE, the beam axis is the z-axis; no need to rotate
        const TLorentzVector &p4fsl = *(event.FinalStatePrimaryLepton()->P4());
        TVector3 p3fsl = p4fsl.Vect();
        double pl = p4fsl.Pz();
        double pt = sqrt(p4fsl.Px() * p4fsl.Px() + p4fsl.Py() * p4fsl.Py());
        double p = p3fsl.Mag();
        double E = p4fsl.E();
        double theta = get_angle(event);


        const double enuQE = (
                pow(constants::kProtonMass, 2) -
                pow((constants::kNeutronMass - Ebinding), 2) - 
                pow(constants::kMuonMass, 2) + 
                2 * (constants::kNeutronMass - Ebinding) * E
                ) / 
            (2 * (constants::kNeutronMass - Ebinding - E + p * cos(theta)));
        const double Q2QE = 2 * enuQE * (E - p * cos(theta)) -
            pow(constants::kMuonMass, 2);

        // double q2_weight = event.XSec() / (units::cm2);
        double q2_bin_wid = get_bin_width(Q2QE, q2_bins, q2_nbins);
        double enu_bin_wid = get_bin_width(enuQE, enu_bins, enu_nbins);
        double pl_bin_wid = get_bin_width(pl, pl_bins, pl_nbins);
        double pt_bin_wid = get_bin_width(pt, pt_bins, pt_nbins);

        // double q2_wt = 1.0 / q2_bin_wid;
        // double enuq2_wt = 1.0 / q2_bin_wid / enu_bin_wid;
        // double ptpl_wt = 1.0 / pl_bin_wid / pt_bin_wid;

        // double q2_wt = 1.0;
        // double enuq2_wt = 1.0;
        // double ptpl_wt = 1.0;

        double units_scale = 1.0e-39;   // 1.0e-39;
        double q2_wt = event.XSec() / (units_scale * units::cm2);
        double enuq2_wt = event.XSec() / (units_scale * units::cm2);
        double ptpl_wt = event.XSec() / (units_scale * units::cm2);

        // double q2_wt = event.XSec() / (units::cm2) / q2_bin_wid;
        // double enuq2_wt = event.XSec() / (units::cm2) / enu_bin_wid / q2_bin_wid;
        // double ptpl_wt = event.XSec() / (units::cm2) / pt_bin_wid / pl_bin_wid;

        if (is_cc_event) {
            n_cc_events += 1;
            // was (1E-39 * units::cm2)
            double eweight = event.XSec() / (units::cm2);
            // double eweight = event.XSec() / (units::cm2) / ebinwid;
            cc_total_cross_section += event.XSec() / (units::cm2);
            q2_cc->Fill(Q2QE, eweight * q2_wt);
            q2_n_cc_events->Fill(Q2QE);
            h_dsigmadE->Fill(enu, eweight);
            // h_flux->Fill(enu, 1.0 / ebinwid);
            h_flux->Fill(enu);
        }

        if (is_ccqe_true_event) {
            q2_true_nocut->Fill(Q2QE, q2_wt);
            q2_n_true_nocut_events->Fill(Q2QE, 1.0);
            enuq2_true_nocut->Fill(Q2QE, enuQE, enuq2_wt);
            ptpl_true_nocut->Fill(pl, pt, ptpl_wt);
            if (angle_cut(event)) {
                q2_true_cut->Fill(Q2QE, q2_wt);
                enuq2_true_cut->Fill(Q2QE, enuQE, enuq2_wt);
                ptpl_true_cut->Fill(pl, pt, ptpl_wt);
            }
        }

        if (is_ccqe_like_event) {
            q2_like_nocut->Fill(Q2QE, q2_wt);
            enuq2_like_nocut->Fill(Q2QE, enuQE, enuq2_wt);
            ptpl_like_nocut->Fill(pl, pt, ptpl_wt);
            if (angle_cut(event)) {
                q2_like_cut->Fill(Q2QE, q2_wt);
                enuq2_like_cut->Fill(Q2QE, enuQE, enuq2_wt);
                ptpl_like_cut->Fill(pl, pt, ptpl_wt);
            }
        }

        // clear current mc event record
        mcrec->Clear();

    } //end loop over events

    //
    // division by "flux" 
    //
    double scale_factor = 1.0;

    double integ_xsec_of_events_naive = h_dsigmadE->Integral();
    double integ_xsec_of_events = h_dsigmadE->Integral("width");
    int n_cc_events3 = h_flux->Integral();
    h_dsigmadE->Divide(h_flux);
    double integ_xsec = h_dsigmadE->Integral("width");
    double integ_xsec_naive = h_dsigmadE->Integral();
    double integ_naive = q2_cc->Integral();
    double integ = q2_cc->Integral("width");
    double nentries = q2_cc->GetEntries();
    scale_factor = integ_xsec / per_nucleon_correction_factor;

    std::cout << "per nucleon corr = " << per_nucleon_correction_factor << std::endl;
    std::cout << "integral-E(total) cross section = " << integ_xsec << std::endl;
    std::cout << "integ-Q2(here) = " << integ
        << "; nentries = " << nentries << std::endl;

    q2_n_cc_events->Sumw2();
    // double scale_factor = cc_total_cross_section / 
    //     (n_cc_events * per_nucleon_correction_factor);
    // std::cout << "n_cc_events = " << n_cc_events << std::endl;
    // double scale_factor = (integ_xsec) /
    //     (n_cc_events * per_nucleon_correction_factor);


    // std::cout << "new scale factor = " << scale_factor << std::endl;

    // scale_factor = 1.0;

    // q2_true_nocut->Divide(q2_n_true_nocut_events);
    // normalize(q2_true_nocut, 1.0 / per_nucleon_correction_factor, q2_bins, q2_nbins);
    int n_events = q2_n_true_nocut_events->Integral();
    int n_cc_events2 = q2_n_cc_events->Integral();
    // double frac = 1.0 * n_events / n_cc_events2;
    // frac = 1.0 * n_events / nev;

    int nq2_true_nocut = q2_true_nocut->GetEntries();
    int nq2_true_cut = q2_true_cut->GetEntries();
    int nq2_like_nocut = q2_like_nocut->GetEntries();
    int nq2_like_cut = q2_like_cut->GetEntries();
    std::cout << "nq2_true_nocut = " << nq2_true_nocut << std::endl;
    std::cout << "nq2_true_cut = " << nq2_true_cut << std::endl;
    std::cout << "nq2_like_nocut = " << nq2_like_nocut << std::endl;
    std::cout << "nq2_like_cut = " << nq2_like_cut << std::endl;

    scale_factor = integ_xsec / per_nucleon_correction_factor / n_cc_events3;

    // scale_factor = q2_cc->Integral("width") / q2_n_cc_events->Integral() / per_nucleon_correction_factor;

    scale_factor = 1 / q2_true_nocut->GetEntries() / per_nucleon_correction_factor;
    normalize(q2_true_nocut, scale_factor);
    normalize(q2_true_cut, scale_factor);
    normalize(q2_like_nocut, scale_factor);
    normalize(q2_like_cut, scale_factor);

    // normalize(q2_true_nocut, q2_n_cc_events, scale_factor);
    // normalize(q2_true_cut, q2_n_cc_events, scale_factor);
    // normalize(q2_like_nocut, q2_n_cc_events, scale_factor);
    // normalize(q2_like_cut, q2_n_cc_events, scale_factor);

    normalize2D(enuq2_true_nocut, scale_factor);
    normalize2D(enuq2_true_cut, scale_factor);
    normalize2D(enuq2_like_nocut, scale_factor);
    normalize2D(enuq2_like_cut, scale_factor);

    int nptpl_true_nocut = ptpl_true_nocut->GetEntries();
    int nptpl_true_cut = ptpl_true_cut->GetEntries();
    int nptpl_like_nocut = ptpl_like_nocut->GetEntries();
    int nptpl_like_cut = ptpl_like_cut->GetEntries();
    std::cout << "nptpl_true_nocut = " << nptpl_true_nocut << std::endl;
    std::cout << "nptpl_true_cut = " << nptpl_true_cut << std::endl;
    std::cout << "nptpl_like_nocut = " << nptpl_like_nocut << std::endl;
    std::cout << "nptpl_like_cut = " << nptpl_like_cut << std::endl;

    normalize2D(ptpl_true_nocut, scale_factor);
    normalize2D(ptpl_true_cut, scale_factor);
    normalize2D(ptpl_like_nocut, scale_factor);
    normalize2D(ptpl_like_cut, scale_factor);

    //
    // Close histogram file and clean up
    //
    std::string histname = "$ANA_HIST_DIR/" + out_file_name;
    TFile * outputfile = new TFile(histname.c_str(), "RECREATE");
    Hlist->Write();
    outputfile->Close();

    delete h_dsigmadE;
    delete h_flux;

    delete q2_cc;
    delete q2_n_cc_events;

    delete q2_true_nocut;
    delete q2_n_true_nocut_events;

    delete q2_true_cut;
    delete q2_like_nocut;
    delete q2_like_cut;

    delete enuq2_true_nocut;
    delete enuq2_true_cut;
    delete enuq2_like_nocut;
    delete enuq2_like_cut;

    delete ptpl_true_nocut;
    delete ptpl_true_cut;
    delete ptpl_like_nocut;
    delete ptpl_like_cut;

    delete[] pt_bins;
    delete[] pl_bins;

    delete chain;

    return 0;
}
