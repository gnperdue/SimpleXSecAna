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
            double width = bins[index_check] - bins[index_check - 1];
            return width;
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
double flux_integral(const TH1D* flux, double emin, double emax) {
    int nbins_flux = flux->GetNbinsX();
    double flux_int = 0.0;
    for (int i = 1; i <= nbins_flux; ++i) {
        double e = flux->GetBinCenter(i);
        if (e > emin && e < emax) {
            double f = flux->GetBinContent(i);
            double w = flux->GetBinWidth(i);
            flux_int += f * w;
        }
    }
    return flux_int;
}

//____________________________________________________________________________
double flux_xsec_integral(const TH1D* xsec, const TH1D* flux) {
    // assume we have identical binning
    int nbins_xsec = xsec->GetNbinsX();
    int nbins_flux = flux->GetNbinsX();
    if (nbins_xsec != nbins_flux) {
        return -1.0;
    }

    double conv_int = 0.0;
    for (int i = 1; i <= nbins_xsec; ++i) {
        double x = xsec->GetBinContent(i);
        double f = flux->GetBinContent(i);
        double w = xsec->GetBinWidth(i);
        conv_int += x * f * w;
    }
    return conv_int;
}

//____________________________________________________________________________
int main(int argc, char ** argv)
{
    Long64_t max_events = 100000;
    std::string list_file_name = "default_nubar_qe_like_scint.txt";
    std::string out_file_name = "default_nubar_qe_like_scint.root";
    int signal_pdg = -14;
    int number_of_nucleons = 13;             // 13 nucleons in CH
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
        else if (sw == "-m" || sw == "--max_events") {
            optind++;
            max_events = Long64_t(atoi(argv[optind]));
        }
        else if (sw == "-s" || sw == "--signal_pdg") {
            optind++;
            signal_pdg = atoi(argv[optind]);
        }
        else if (sw == "-o" || sw == "--out_file") {
            optind++;
            out_file_name = argv[optind];
        }
        else if (sw == "-c" || sw == "--nucleon_correction") {
            optind++;
            per_nucleon_correction_factor = atoi(argv[optind]);
        }
        else if (sw == "-u" || sw == "--number_of_nucleons") {
            optind++;
            number_of_nucleons = atoi(argv[optind]);
        }
        else if (sw == "-n" || sw == "--flux_e_min") {
            optind++;
            flux_e_min = atof(argv[optind]);
        }
        else if (sw == "-x" || sw == "--flux_e_max") {
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

    //
    // Get the input flux
    //
    std::string numubar_fluxfile = s + std::string("/numubarRHC.root");
    std::string numu_fluxfile = s + std::string("/numuFHC.root");
    TFile numubar_f(numubar_fluxfile.c_str());
    TH1D *numubar_flux = (TH1D*)numubar_f.Get("newnumubarRHC_CV_WithStatErr");
    TFile numu_f(numu_fluxfile.c_str());
    TH1D *numu_flux = (TH1D*)numu_f.Get("newnumuFHC_CV_WithStatErr");

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

    // keep a copy of the flux used to produce the samples
    Hlist->Add(numubar_flux);
    Hlist->Add(numu_flux);

    // copy the flux files to get the binning we want - we can be lazy since
    // both nu and nu-bar flux histograms have the same binning, so just use
    // whichever one
    TH1D *h_dsigmadE = (TH1D*)numubar_flux->Clone("h_dsigmadE");
    h_dsigmadE->Reset();
    h_dsigmadE->SetName("h_dsigmadE");
    h_dsigmadE->SetTitle("Total Cross Section;E_{#nu} (GeV);#frac{d#sigma}{dE} #times 10^{39} cm^{2} / GeV");
    TH1D *h_selectedsamp_flux = (TH1D*)numubar_flux->Clone("h_selectedsamp_flux");
    h_selectedsamp_flux->Reset();
    h_selectedsamp_flux->SetName("h_selectedsamp_flux");
    h_selectedsamp_flux->SetTitle("Selected events;E;n");

    Hlist->Add(h_dsigmadE);
    Hlist->Add(h_selectedsamp_flux);

    char axes1[200];

    // CC Inclusive

    sprintf(axes1, "Differential cross section - CC;Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_cc = new TH1D("q2_cc", axes1, q2_nbins, q2_bins);

    sprintf(axes1, "Number of CC Events;Q^{2}_{QE} (GeV^{2});N-Events");
    TH1D* q2_n_cc_events = new TH1D("q2_n_cc_events", axes1, q2_nbins, q2_bins);

    // CCQE "True"

    sprintf(axes1, "Differential cross section - true CCQE or MEC;Q^{2}_{QE} (GeV^{2});Cross section times 10^{39} cm^{2}");
    TH1D* q2_true_nocut = new TH1D("q2_true_nocut", axes1, q2_nbins, q2_bins);

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
    int n_ccqe_events = 0;
    int n_ccqelike_events = 0;
    double units_scale = 1.0;   // 1.0e-39;
    for(Long64_t i = 0; i < nev; i++) {

        // get next tree entry
        chain->GetEntry(i);
        if (i % 20000 == 0) {
            std::cout << "Event " << i << std::endl;
        }

        // get the GENIE event
        EventRecord &event = *(mcrec->event);
        GHepParticle *nu = event.Probe();
        double enu = nu->Energy();
        // is this a signal event?
        bool is_cc_event = is_cc_numu_numubar(event, signal_pdg);
        bool is_ccqe_true_event = is_ccqe_true(event, signal_pdg);
        bool is_ccqe_like_event = is_ccqe_like(event, signal_pdg);

        double evt_weight = event.XSec() / (units_scale * units::cm2);
        // bin width normalization falls out when we divide by h_selectedsamp_flux later
        h_dsigmadE->Fill(enu, evt_weight);
        h_selectedsamp_flux->Fill(enu);

        // only accept events inside the specified flux window
        if ((enu < flux_e_min) || (enu > flux_e_max)) {
            mcrec->Clear();
            continue;
        } 

        // only accept CC events
        if (!is_cc_event) {
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

        double q2_bin_wid = get_bin_width(Q2QE, q2_bins, q2_nbins);        
        double enu_bin_wid = get_bin_width(enuQE, enu_bins, enu_nbins);        
        double pl_bin_wid = get_bin_width(pl, pl_bins, pl_nbins);        
        double pt_bin_wid = get_bin_width(pt, pt_bins, pt_nbins);        

        double q2_wt = 1.0;
        double enuq2_wt = 1.0;
        double ptpl_wt = 1.0;

        if (is_cc_event) {
            n_cc_events += 1;
            q2_cc->Fill(Q2QE, q2_wt);
            q2_n_cc_events->Fill(Q2QE);
        }

        if (is_ccqe_true_event) {
            n_ccqe_events += 1;
            q2_true_nocut->Fill(Q2QE, q2_wt);
            enuq2_true_nocut->Fill(Q2QE, enuQE, enuq2_wt);
            ptpl_true_nocut->Fill(pl, pt, ptpl_wt);
            if (angle_cut(event)) {
                q2_true_cut->Fill(Q2QE, q2_wt);
                enuq2_true_cut->Fill(Q2QE, enuQE, enuq2_wt);
                ptpl_true_cut->Fill(pl, pt, ptpl_wt);
            }
        }

        if (is_ccqe_like_event) {
            n_ccqelike_events += 1;
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

    std::cout << "n_cc_events       = " << n_cc_events << std::endl;
    std::cout << "n_ccqe_events     = " << n_ccqe_events << std::endl;
    std::cout << "n_ccqelike_events = " << n_ccqelike_events << std::endl;

    //
    // division by "flux" 
    // - selected events, not bin-by-bin if not doing sigma(E)
    //
    h_dsigmadE->Divide(h_selectedsamp_flux);
    double conv_int = flux_xsec_integral(h_dsigmadE, numubar_flux);
    double flux_int = flux_integral(numubar_flux, flux_e_min, flux_e_max);
    double nucleon_scaling_factor = double(number_of_nucleons) /
        per_nucleon_correction_factor;
    double w = conv_int / (flux_int * double(nev)) * nucleon_scaling_factor;

    q2_cc->Scale(w, "width");
    q2_true_nocut->Scale(w, "width");
    q2_true_cut->Scale(w, "width");
    q2_like_nocut->Scale(w, "width");
    q2_like_cut->Scale(w, "width");

    enuq2_true_nocut->Scale(w, "width");
    enuq2_true_cut->Scale(w, "width");
    enuq2_like_nocut->Scale(w, "width");
    enuq2_like_cut->Scale(w, "width");

    ptpl_true_nocut->Scale(w, "width");
    ptpl_true_cut->Scale(w, "width");
    ptpl_like_nocut->Scale(w, "width");
    ptpl_like_cut->Scale(w, "width");

    //
    // Close histogram file and clean up
    //
    std::string histname = "$ANA_HIST_DIR/" + out_file_name;
    TFile * outputfile = new TFile(histname.c_str(), "RECREATE");
    Hlist->Write();
    outputfile->Close();

    delete h_dsigmadE;
    delete h_selectedsamp_flux;

    delete q2_cc;
    delete q2_n_cc_events;

    delete q2_true_nocut;
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
