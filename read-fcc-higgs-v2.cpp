#include <stdio.h>
#include <stdlib.h>
#include "Delphes.C"
#include <TMath.h>
#include <TTree.h>
#include <glob.h>
#include <TError.h>
#include <vector>

std::vector<std::string> glob(const char *pattern) {
    glob_t g;
    glob(pattern, GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
    std::vector<std::string> filelist;
    filelist.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        filelist.emplace_back(g.gl_pathv[i]);
    }
    globfree(&g);
    return filelist;
}

Float_t deltaPhi(Float_t phi1, Float_t phi2)
{
    Float_t dphi = phi1 - phi2;
    while (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
    return TMath::Abs(dphi);
}

vector<int> find_ele(Delphes *indelphes, double ptcut, int muon_index)
{
    vector<int> res;
    for (int e=0; e<indelphes->Electron_size; e++)
    {
        if (indelphes->Electron_PT[e] < ptcut) continue;
        if (TMath::Abs(indelphes->Electron_Eta[e]) > 6.0) continue;
        //if (TMath::Abs(indelphes->Electron_Eta[e]) > 1.44 && TMath::Abs(indelphes->Electron_Eta[e]) < 1.57) continue;
        //if (indelphes->Electron_IsolationVar[e] < 0.1) continue;
        if (muon_index != -1)
        {
            double deltaR = 0;
            deltaR += TMath::Power((indelphes->Electron_Eta[e] - indelphes->Muon_Eta[muon_index]), 2);
            deltaR += TMath::Power((indelphes->Electron_Phi[e] - indelphes->Muon_Phi[muon_index]), 2);
            deltaR = TMath::Sqrt(deltaR);
            if (deltaR < 0.3) continue;
            if (indelphes->Electron_Charge[e] == indelphes->Muon_Charge[muon_index]) continue;
        }
        res.push_back(e);
    }
    return res;
}

vector<int> find_mu(Delphes *indelphes, double ptcut, int electron_index)
{
    vector<int> res;
    for (int mu=0; mu<indelphes->Muon_size; mu++)
    {
        if (indelphes->Muon_PT[mu] < ptcut) continue;
        if (TMath::Abs(indelphes->Muon_Eta[mu]) > 6.0) continue;
        //if (indelphes->Muon_IsolationVar[mu] > 0.15) continue;

        if (electron_index != -1)
        {
            double deltaR = 0;
            deltaR += TMath::Power((indelphes->Muon_Eta[mu] - indelphes->Electron_Eta[electron_index]), 2);
            deltaR += TMath::Power((indelphes->Muon_Phi[mu] - indelphes->Electron_Phi[electron_index]), 2);
            deltaR = TMath::Sqrt(deltaR);
            if (deltaR < 0.3) continue;
            if (indelphes->Muon_Charge[mu] == indelphes->Electron_Charge[electron_index]) continue;
        }
        res.push_back(mu);
    }
    return res;
}

class PlotSet
{
    public:
        PlotSet() { }
        void AddHist(TH1D* hist)
        {
            histograms.push_back(hist);
        }
        void SaveAll(TFile *outfile)
        {
            outfile->cd();
            for (TH1* hist: histograms) hist->Write();
        }
        void PrimeFill(double *value)
        {
            primed_variable = value;
        }
        void Fill(int histnum)
        {
            if (histnum < histograms.size()) histograms[histnum]->Fill(*primed_variable);
        }
        vector<TH1*> histograms;

    private:
        double *primed_variable;
};

const int HIST_BINS = 20;
const double HIST_START = 0;
const double HIST_END   = 200;

void add_hist_shorthand(PlotSet *plots, TString histname, TString propername)
{
    plots->AddHist(new TH1D(histname, propername, HIST_BINS, HIST_START, HIST_END));
}

void read_fcc_higgs_v2(TString infilename, TString outfilename)
{
    gErrorIgnoreLevel = kFatal;

    TChain *intree = new TChain("Delphes");
    for (const auto &filename : glob(infilename.Data()))
    {
        printf("Reading %s\n", filename.c_str());
        //TFile *infile = new TFile(filename.c_str());
        intree->Add(filename.c_str());
    }

    Delphes *indelphes = new Delphes(intree);

    PlotSet plots_mutaue;
    PlotSet plots_etaumu;

    add_hist_shorthand(&plots_mutaue, "mutau_e_step01", "mutau_e no b-jets");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step02", "mutau_e 0, 1 jet"); // 1
    add_hist_shorthand(&plots_mutaue, "mutau_e_step03_0j", "mutau_e 0 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step03_1j", "mutau_e 1 jet"); // 3
    add_hist_shorthand(&plots_mutaue, "mutau_e_step04_0j", "mutau_e 1+ muon 0 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step04_1j", "mutau_e 1+ muon 1 jet"); // 5
    add_hist_shorthand(&plots_mutaue, "mutau_e_step05_0j", "mutau_e 1 muon 0 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step05_1j", "mutau_e 1 muon 1 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step06_0j", "mutau_e 1+ electron 0 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step06_1j", "mutau_e 1+ electron 1 jet"); // 9
    add_hist_shorthand(&plots_mutaue, "mutau_e_step07_0j", "mutau_e 1 electron 0 jet"); // 10
    add_hist_shorthand(&plots_mutaue, "mutau_e_step07_1j", "mutau_e 1 electron 1 jet"); // 11
    add_hist_shorthand(&plots_mutaue, "mutau_e_step08_0j", "mutau_e min pT 0 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step08_1j", "mutau_e min pT 1 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step09_0j", "mutau_e max deltaPhi e, met 0 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step09_1j", "mutau_e max deltaPhi e, met 1 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step10_0j", "mutau_e min deltaPhi e, mu 0 jet");
    add_hist_shorthand(&plots_mutaue, "mutau_e_step10_1j", "mutau_e min deltaPhi e, mu 1 jet");

    add_hist_shorthand(&plots_etaumu, "etau_mu_step01",    "etau_mu no b-jets");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step02",    "etau_mu 0, 1 jet"); // 1
    add_hist_shorthand(&plots_etaumu, "etau_mu_step03_0j", "etau_mu 0 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step03_1j", "etau_mu 1 jet"); // 3
    add_hist_shorthand(&plots_etaumu, "etau_mu_step04_0j", "etau_mu 1+ electron 0 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step04_1j", "etau_mu 1+ electron 1 jet"); // 5
    add_hist_shorthand(&plots_etaumu, "etau_mu_step05_0j", "etau_mu 1 electron 0 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step05_1j", "etau_mu 1 electron 1 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step06_0j", "etau_mu 1+ muon 0 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step06_1j", "etau_mu 1+ muon 1 jet"); // 9
    add_hist_shorthand(&plots_etaumu, "etau_mu_step07_0j", "etau_mu 1 muon 0 jet"); // 10
    add_hist_shorthand(&plots_etaumu, "etau_mu_step07_1j", "etau_mu 1 muon 1 jet"); // 11
    add_hist_shorthand(&plots_etaumu, "etau_mu_step08_0j", "etau_mu min pT 0 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step08_1j", "etau_mu min pT 1 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step09_0j", "etau_mu max deltaPhi mu, met 0 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step09_1j", "etau_mu max deltaPhi mu, met 1 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step10_0j", "etau_mu min deltaPhi e, mu 0 jet");
    add_hist_shorthand(&plots_etaumu, "etau_mu_step10_1j", "etau_mu min deltaPhi e, mu 1 jet");

    bool plotthis_mutaue[18];
    bool plotthis_etaumu[18];
    double mass_collinear_mutaue = 0;
    double mass_collinear_etaumu = 0;
    plots_mutaue.PrimeFill(&mass_collinear_mutaue);
    plots_etaumu.PrimeFill(&mass_collinear_etaumu);
    TLorentzVector p4_tau, p4_lepton;
    TLorentzVector p4_muon, p4_electron, p4_met;
    double pT_nu_est, x_vis_tau;
    double deltaPhi_e_met, deltaPhi_mu_met, deltaPhi_e_mu;

    for (Long64_t ievent=0; ievent < intree->GetEntries(); ievent++)
    {
        if (ievent % 10000 == 0) printf("Reading event %lld\n", ievent);
        indelphes->GetEntry(ievent);

        ///////////////////////////////////////
        // mu + tau_e
        ///////////////////////////////////////

        for (int i=0; i<18; i++) plotthis_mutaue[i] = false;

        vector<int> passed_jets;
        vector<int> passed_b_jets;
        vector<int> muon_vec;
        vector<int> electron_vec;
        int only_mu  = -1;
        int only_ele = -1;

        for (int j=0; j<indelphes->Jet_size; j++)
        {
            if (indelphes->Jet_PT[j] < 30) continue;
            if (TMath::Abs(indelphes->Jet_Eta[j]) > 6.0) continue;
            passed_jets.push_back(j);
            if (indelphes->Jet_BTag[j] & 0b111) passed_b_jets.push_back(j);
        }
        //if (passed_b_jets.size() == 0) plotthis_mutaue[0] = true;
        plotthis_mutaue[0] = passed_b_jets.size() == 0;
        plotthis_mutaue[1] = plotthis_mutaue[0] and passed_jets.size() <= 1;
        if (plotthis_mutaue[1])
        {
            plotthis_mutaue[2] = passed_jets.size() == 0;
            plotthis_mutaue[3] = passed_jets.size() == 1;
            plotthis_mutaue[4] = plotthis_mutaue[2];
            plotthis_mutaue[5] = plotthis_mutaue[3];
        }
        if (plotthis_mutaue[2] or plotthis_mutaue[3])
        {
            muon_vec = find_mu(indelphes, 53, -1);
            if (muon_vec.size() == 0)
            {
                plotthis_mutaue[4] = false;
                plotthis_mutaue[5] = false;
            }
            plotthis_mutaue[6] = plotthis_mutaue[4];
            plotthis_mutaue[7] = plotthis_mutaue[5];
            if (muon_vec.size() > 1)
            {
                plotthis_mutaue[6] = false;
                plotthis_mutaue[7] = false;
            }
            plotthis_mutaue[8] = plotthis_mutaue[6];
            plotthis_mutaue[9] = plotthis_mutaue[7];
        }
        if (plotthis_mutaue[6] or plotthis_mutaue[7])
        {
            electron_vec = find_ele(indelphes, 10, plotthis_mutaue[6] or plotthis_mutaue[7] ? muon_vec[0] : -1);
            if (electron_vec.size() == 0)
            {
                plotthis_mutaue[8] = false;
                plotthis_mutaue[9] = false;
            }
            plotthis_mutaue[10] = plotthis_mutaue[8];
            plotthis_mutaue[11] = plotthis_mutaue[9];
            if (electron_vec.size() > 1)
            {
                plotthis_mutaue[10] = false;
                plotthis_mutaue[11] = false;
            }
            plotthis_mutaue[12] = plotthis_mutaue[10];
            plotthis_mutaue[13] = plotthis_mutaue[11];
        }
        if (plotthis_mutaue[10] or plotthis_mutaue[11])
        {
            only_mu = muon_vec[0];
            only_ele = electron_vec[0];
            if (indelphes->Muon_PT[only_mu] < 60)
            {
                plotthis_mutaue[12] = false;
                plotthis_mutaue[13] = false;
            }
            plotthis_mutaue[14] = plotthis_mutaue[12];
            plotthis_mutaue[15] = plotthis_mutaue[13];
        }
        if (plotthis_mutaue[12] or plotthis_mutaue[13])
        {
            deltaPhi_e_met = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->MissingET_Phi[0]);
            if (deltaPhi_e_met > 0.7)
            {
                plotthis_mutaue[14] = false;
                plotthis_mutaue[15] = false;
            }
            plotthis_mutaue[16] = plotthis_mutaue[14];
            plotthis_mutaue[17] = plotthis_mutaue[15];
        }
        if (plotthis_mutaue[14] or plotthis_mutaue[15])
        {
            deltaPhi_e_mu = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->Muon_Phi[only_mu]);
            if (deltaPhi_e_mu < 2.2)
            {
                plotthis_mutaue[16] = false;
                plotthis_mutaue[17] = false;
            }
        }

        //TLorentzVector p4_tau;
        //TLorentzVector p4_lepton;
        if (plotthis_mutaue[16] or plotthis_mutaue[17])
        {

            p4_tau.SetPtEtaPhiM(indelphes->Electron_PT[only_ele], indelphes->Electron_Eta[only_ele], indelphes->Electron_Phi[only_ele], 0.000511);
            p4_lepton.SetPtEtaPhiM(indelphes->Muon_PT[only_mu], indelphes->Muon_Eta[only_mu], indelphes->Muon_Phi[only_mu], 0.10566);
        }
        else
        {
            //TLorentzVector p4_muon, p4_electron, p4_met;
            if (indelphes->Muon_size > 0) p4_muon.SetPtEtaPhiM(indelphes->Muon_PT[0], indelphes->Muon_Eta[0], indelphes->Muon_Phi[0], 0.10566);
            else p4_muon.SetPtEtaPhiM(0, 0, 0, 0.10566);
            if (indelphes->Electron_size > 0) p4_electron.SetPtEtaPhiM(indelphes->Electron_PT[0], indelphes->Electron_Eta[0], indelphes->Electron_Phi[0], 0.000511);
            else p4_electron.SetPtEtaPhiM(0, 0, 0, 0.000511);
            p4_met.SetPtEtaPhiM(indelphes->MissingET_MET[0], 0, indelphes->MissingET_Phi[0], 0);
            if (p4_muon.DeltaR(p4_met) < p4_electron.DeltaR(p4_met))
            {
                p4_tau.SetPtEtaPhiM(p4_muon.Pt(), p4_muon.Eta(), p4_muon.Phi(), 0.10566);
                p4_lepton.SetPtEtaPhiM(p4_electron.Pt(), p4_electron.Eta(), p4_electron.Phi(), 0.000511);
            }
            else
            {
                p4_lepton.SetPtEtaPhiM(p4_muon.Pt(), p4_muon.Eta(), p4_muon.Phi(), 0.10566);
                p4_tau.SetPtEtaPhiM(p4_electron.Pt(), p4_electron.Eta(), p4_electron.Phi(), 0.000511);
            }
        }
        pT_nu_est = indelphes->MissingET_MET[0] * TMath::Cos(deltaPhi(indelphes->MissingET_Phi[0], p4_tau.Phi()));
        x_vis_tau = p4_tau.Pt() / (p4_tau.Pt() + pT_nu_est);
        mass_collinear_mutaue = (p4_tau+p4_lepton).M() / TMath::Sqrt(x_vis_tau);

        for (int i=0; i<18; i++) if (plotthis_mutaue[i]) plots_mutaue.Fill(i);

        ///////////////////////////////////////
        // e + tau_mu
        ///////////////////////////////////////

        for (int i=0; i<18; i++) plotthis_etaumu[i] = false;

        passed_jets.clear();
        passed_b_jets.clear();
        muon_vec.clear();
        electron_vec.clear();
        only_mu  = -1;
        only_ele = -1;

        for (int j=0; j<indelphes->Jet_size; j++)
        {
            if (indelphes->Jet_PT[j] < 30) continue;
            if (TMath::Abs(indelphes->Jet_Eta[j]) > 6.0) continue;
            passed_jets.push_back(j);
            if (indelphes->Jet_BTag[j] & 0b111) passed_b_jets.push_back(j);
        }
        //if (passed_b_jets.size() == 0) plotthis_etaumu[0] = true;
        plotthis_etaumu[0] = passed_b_jets.size() == 0;
        plotthis_etaumu[1] = plotthis_etaumu[0] and passed_jets.size() <= 1;
        if (plotthis_etaumu[1])
        {
            plotthis_etaumu[2] = passed_jets.size() == 0;
            plotthis_etaumu[3] = passed_jets.size() == 1;
            plotthis_etaumu[4] = plotthis_etaumu[2];
            plotthis_etaumu[5] = plotthis_etaumu[3];
        }
        if (plotthis_etaumu[2] or plotthis_etaumu[3])
        {
            electron_vec = find_ele(indelphes, 26, -1);
            if (electron_vec.size() == 0)
            {
                plotthis_etaumu[4] = false;
                plotthis_etaumu[5] = false;
            }
            plotthis_etaumu[6] = plotthis_etaumu[4];
            plotthis_etaumu[7] = plotthis_etaumu[5];
            if (electron_vec.size() > 1)
            {
                plotthis_etaumu[6] = false;
                plotthis_etaumu[7] = false;
            }
            plotthis_etaumu[8] = plotthis_etaumu[6];
            plotthis_etaumu[9] = plotthis_etaumu[7];
        }
        if (plotthis_etaumu[6] or plotthis_etaumu[7])
        {
            muon_vec = find_mu(indelphes, 10, plotthis_etaumu[6] or plotthis_etaumu[7] ? electron_vec[0] : -1);
            if (muon_vec.size() == 0)
            {
                plotthis_etaumu[8] = false;
                plotthis_etaumu[9] = false;
            }
            plotthis_etaumu[10] = plotthis_etaumu[8];
            plotthis_etaumu[11] = plotthis_etaumu[9];
            if (muon_vec.size() > 1)
            {
                plotthis_etaumu[10] = false;
                plotthis_etaumu[11] = false;
            }
            plotthis_etaumu[12] = plotthis_etaumu[10];
            plotthis_etaumu[13] = plotthis_etaumu[11];
        }
        if (plotthis_etaumu[10] or plotthis_etaumu[11])
        {
            only_mu = muon_vec[0];
            only_ele = electron_vec[0];
            if (indelphes->Electron_PT[only_ele] < 60)
            {
                plotthis_etaumu[12] = false;
                plotthis_etaumu[13] = false;
            }
            plotthis_etaumu[14] = plotthis_etaumu[12];
            plotthis_etaumu[15] = plotthis_etaumu[13];
        }
        if (plotthis_etaumu[12] or plotthis_etaumu[13])
        {
            deltaPhi_mu_met = deltaPhi(indelphes->Muon_Phi[only_mu], indelphes->MissingET_Phi[0]);
            if (deltaPhi_mu_met > 0.7)
            {
                plotthis_etaumu[14] = false;
                plotthis_etaumu[15] = false;
            }
            plotthis_etaumu[16] = plotthis_etaumu[14];
            plotthis_etaumu[17] = plotthis_etaumu[15];
        }
        if (plotthis_etaumu[14] or plotthis_etaumu[15])
        {
            deltaPhi_e_mu = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->Muon_Phi[only_mu]);
            if (deltaPhi_e_mu < 2.2)
            {
                plotthis_etaumu[16] = false;
                plotthis_etaumu[17] = false;
            }
        }

        //TLorentzVector p4_tau;
        //TLorentzVector p4_lepton;
        if (plotthis_etaumu[16] or plotthis_etaumu[17])
        {

            p4_lepton.SetPtEtaPhiM(indelphes->Electron_PT[only_ele], indelphes->Electron_Eta[only_ele], indelphes->Electron_Phi[only_ele], 0.000511);
            p4_tau.SetPtEtaPhiM(indelphes->Muon_PT[only_mu], indelphes->Muon_Eta[only_mu], indelphes->Muon_Phi[only_mu], 0.10566);
        }
        else
        {
            //TLorentzVector p4_muon, p4_electron, p4_met;
            if (indelphes->Muon_size > 0) p4_muon.SetPtEtaPhiM(indelphes->Muon_PT[0], indelphes->Muon_Eta[0], indelphes->Muon_Phi[0], 0.10566);
            else p4_muon.SetPtEtaPhiM(0, 0, 0, 0.10566);
            if (indelphes->Electron_size > 0) p4_electron.SetPtEtaPhiM(indelphes->Electron_PT[0], indelphes->Electron_Eta[0], indelphes->Electron_Phi[0], 0.000511);
            else p4_electron.SetPtEtaPhiM(0, 0, 0, 0.000511);
            p4_met.SetPtEtaPhiM(indelphes->MissingET_MET[0], 0, indelphes->MissingET_Phi[0], 0);
            if (p4_muon.DeltaR(p4_met) < p4_electron.DeltaR(p4_met))
            {
                p4_tau.SetPtEtaPhiM(p4_muon.Pt(), p4_muon.Eta(), p4_muon.Phi(), 0.10566);
                p4_lepton.SetPtEtaPhiM(p4_electron.Pt(), p4_electron.Eta(), p4_electron.Phi(), 0.000511);
            }
            else
            {
                p4_lepton.SetPtEtaPhiM(p4_muon.Pt(), p4_muon.Eta(), p4_muon.Phi(), 0.10566);
                p4_tau.SetPtEtaPhiM(p4_electron.Pt(), p4_electron.Eta(), p4_electron.Phi(), 0.000511);
            }
        }
        pT_nu_est = indelphes->MissingET_MET[0] * TMath::Cos(deltaPhi(indelphes->MissingET_Phi[0], p4_tau.Phi()));
        x_vis_tau = p4_tau.Pt() / (p4_tau.Pt() + pT_nu_est);
        mass_collinear_etaumu = (p4_tau+p4_lepton).M() / TMath::Sqrt(x_vis_tau);

        for (int i=0; i<18; i++) if (plotthis_etaumu[i]) plots_etaumu.Fill(i);
    }

    TFile *outfile = new TFile(outfilename, "RECREATE");
    plots_mutaue.SaveAll(outfile);
    plots_etaumu.SaveAll(outfile);
    outfile->Close();
}