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

vector<int> find_ele(Delphes *indelphes, double ptcut)
{
    vector<int> res;
    for (int e=0; e<indelphes->Electron_size; e++)
    {
        if (indelphes->Electron_PT[e] < ptcut) continue;
        if (TMath::Abs(indelphes->Electron_Eta[e]) > 2.5) continue;
        //if (TMath::Abs(indelphes->Electron_Eta[e]) > 1.44 && TMath::Abs(indelphes->Electron_Eta[e]) < 1.57) continue;
        //if (indelphes->Electron_IsolationVar[e] < 0.1) continue;
        res.push_back(e);
    }
    return res;
}

vector<int> find_mu(Delphes *indelphes, double ptcut)
{
    vector<int> res;
    for (int mu=0; mu<indelphes->Muon_size; mu++)
    {
        if (indelphes->Muon_PT[mu] < ptcut) continue;
        if (TMath::Abs(indelphes->Muon_Eta[mu]) > 2.4) continue;
        //if (indelphes->Muon_IsolationVar[mu] > 0.15) continue;
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

void read_fcc_higgs_v2(TString infilename, TString outfilename, bool is_mutaue)
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

    PlotSet plots;

    add_hist_shorthand(&plots, "mutau_e_step01", "mutau_e no b-jets");
    add_hist_shorthand(&plots, "mutau_e_step02", "mutau_e 0, 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step03_0j", "mutau_e 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step03_1j", "mutau_e 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step04_0j", "mutau_e 1+ muon 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step04_1j", "mutau_e 1+ muon 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step05_0j", "mutau_e 1 muon 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step05_1j", "mutau_e 1 muon 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step06_0j", "mutau_e 1+ electron 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step06_1j", "mutau_e 1+ electron 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step07_0j", "mutau_e 1 electron 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step07_1j", "mutau_e 1 electron 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step08_0j", "mutau_e min pT 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step08_1j", "mutau_e min pT 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step09_0j", "mutau_e max deltaPhi e, met 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step09_1j", "mutau_e max deltaPhi e, met 1 jet");
    add_hist_shorthand(&plots, "mutau_e_step10_0j", "mutau_e min deltaPhi e, mu 0 jet");
    add_hist_shorthand(&plots, "mutau_e_step10_1j", "mutau_e min deltaPhi e, mu 1 jet");

    bool plotthis[18];
    double mass_collinear = 0;
    plots.PrimeFill(&mass_collinear);

    for (Long64_t ievent=0; ievent < intree->GetEntries(); ievent++)
    {
        if (ievent % 10000 == 0) printf("Reading event %lld\n", ievent);
        indelphes->GetEntry(ievent);

        for (int i=0; i<18; i++) plotthis[i] = false;

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
        //if (passed_b_jets.size() == 0) plotthis[0] = true;
        plotthis[0] = passed_b_jets.size() == 0;
        plotthis[1] = plotthis[0] and passed_jets.size() <= 1;
        if (plotthis[1])
        {
            plotthis[2] = passed_jets.size() == 0;
            plotthis[3] = passed_jets.size() == 1;
            plotthis[4] = plotthis[2];
            plotthis[5] = plotthis[3];
        }
        if (plotthis[2] or plotthis[3])
        {
            muon_vec = find_mu(indelphes, 53);
            if (muon_vec.size() == 0)
            {
                plotthis[4] = false;
                plotthis[5] = false;
            }
            plotthis[6] = plotthis[4];
            plotthis[7] = plotthis[5];
            if (muon_vec.size() > 1)
            {
                plotthis[6] = false;
                plotthis[7] = false;
            }
            plotthis[8] = plotthis[6];
            plotthis[9] = plotthis[7];
        }
        if (plotthis[6] or plotthis[7])
        {
            electron_vec = find_ele(indelphes, 10);
            if (electron_vec.size() == 0)
            {
                plotthis[8] = false;
                plotthis[9] = false;
            }
            plotthis[10] = plotthis[8];
            plotthis[11] = plotthis[9];
            if (electron_vec.size() > 1)
            {
                plotthis[10] = false;
                plotthis[11] = false;
            }
            plotthis[12] = plotthis[10];
            plotthis[13] = plotthis[11];
        }
        if (plotthis[10] or plotthis[11])
        {
            only_mu = muon_vec[0];
            only_ele = electron_vec[0];
            if (indelphes->Muon_PT[only_mu] < 60)
            {
                plotthis[12] = false;
                plotthis[13] = false;
            }
            plotthis[14] = plotthis[12];
            plotthis[15] = plotthis[13];
        }
        if (plotthis[12] or plotthis[13])
        {
            double deltaPhi_e_met = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->MissingET_Phi[0]);
            if (deltaPhi_e_met > 0.7)
            {
                plotthis[14] = false;
                plotthis[15] = false;
            }
            plotthis[16] = plotthis[14];
            plotthis[17] = plotthis[15];
        }
        if (plotthis[14] or plotthis[15])
        {
            double deltaPhi_e_mu = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->Muon_Phi[only_mu]);
            if (deltaPhi_e_mu < 2.2)
            {
                plotthis[16] = false;
                plotthis[17] = false;
            }
        }

        TLorentzVector p4_tau;
        TLorentzVector p4_lepton;
        if (plotthis[16] or plotthis[17])
        {

            p4_tau.SetPtEtaPhiM(indelphes->Electron_PT[only_ele], indelphes->Electron_Eta[only_ele], indelphes->Electron_Phi[only_ele], 0.000511);
            p4_lepton.SetPtEtaPhiM(indelphes->Muon_PT[only_mu], indelphes->Muon_Eta[only_mu], indelphes->Muon_Phi[only_mu], 0.10566);
        }
        else
        {
            TLorentzVector p4_muon, p4_electron, p4_met;
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
        double pT_nu_est = indelphes->MissingET_MET[0] * TMath::Cos(deltaPhi(indelphes->MissingET_Phi[0], p4_tau.Phi()));
        double x_vis_tau = p4_tau.Pt() / (p4_tau.Pt() + pT_nu_est);
        mass_collinear = (p4_tau+p4_lepton).M() / TMath::Sqrt(x_vis_tau);

        for (int i=0; i<18; i++) if (plotthis[i]) plots.Fill(i);
    }

    TFile *outfile = new TFile(outfilename, "RECREATE");
    plots.SaveAll(outfile);
    outfile->Close();
}