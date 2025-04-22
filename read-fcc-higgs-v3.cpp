#include <stdio.h>
#include <stdlib.h>
#include "Delphes.C"
#include <TMath.h>
#include <TTree.h>
#include <glob.h>
#include <TError.h>
#include <vector>
#include <TString.h>

using namespace std;

const int HIST_BINS = 200;
const double HIST_START = 0;
const double HIST_END   = 2000;
const int MAX_JETS = 2;

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
        int AddHist(TH1D* hist)
        {
            histograms.push_back(hist);
            return histograms.size() - 1;
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

void add_hist_shorthand(PlotSet *plots, TString histname, TString propername, vector<int> *histvector)
{
    int num = plots->AddHist(new TH1D(histname, propername, HIST_BINS, HIST_START, HIST_END));
    histvector->push_back(num);
}

void read_fcc_higgs_v3(TString infilename, TString outfilename)
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

    vector<int> histogram_numbers_mutaue_inclusive;
    vector<int> histogram_numbers_etaumu_inclusive;

    vector<vector<int>> histogram_numbers_mutaue;
    vector<vector<int>> histogram_numbers_etaumu;
    for (int jet = 0; jet <= MAX_JETS; jet++) 
    {
        vector<int> vjet;
        histogram_numbers_mutaue.push_back(vjet);
        histogram_numbers_etaumu.push_back(vjet);
    }

    /*
    histogram_numbers contain the histogram number reference 
    in PlotSet object for each jet requirement and step, 
    where the first index (referenced to vector<int>) 
    is the vector for each jet, and the second index 
    (referenced to int) represents the actual number in PlotSet object.

    Steps are as follows per each jet number
    0: separated jet
    1: 1+ muon
    2: 1 muon
    3: 1+ electron
    4: 1 electron
    5: min pT
    6: max delTaPhi e, met
    7: min deltaPhi e, mu
    8: high mass
    9: low mass
    */

    add_hist_shorthand(&plots_mutaue, "mutau_e_step00", "mutau_e no cuts", &histogram_numbers_mutaue_inclusive);
    add_hist_shorthand(&plots_mutaue, "mutau_e_step01", "mutau_e no b-jets", &histogram_numbers_mutaue_inclusive);
    add_hist_shorthand(&plots_mutaue, "mutau_e_step02", "mutau_e 0, 1 jet", &histogram_numbers_mutaue_inclusive);

    add_hist_shorthand(&plots_etaumu, "etau_mu_step00", "etau_mu no cuts", &histogram_numbers_etaumu_inclusive);
    add_hist_shorthand(&plots_etaumu, "etau_mu_step01", "etau_mu no b-jets", &histogram_numbers_etaumu_inclusive);
    add_hist_shorthand(&plots_etaumu, "etau_mu_step02", "etau_mu 0, 1 jet", &histogram_numbers_etaumu_inclusive);

    for (int jet=0; jet<=MAX_JETS; jet++)
    {
        vector<int> *mutaue_this = &(histogram_numbers_mutaue[jet]);
        vector<int> *etaumu_this = &(histogram_numbers_etaumu[jet]);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step03_%dj", jet), Form("mutau_e %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step04_%dj", jet), Form("mutau_e 1+ muon %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step05_%dj", jet), Form("mutau_e 1 muon %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step06_%dj", jet), Form("mutau_e 1+ electron %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step07_%dj", jet), Form("mutau_e 1 electron %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step08_%dj", jet), Form("mutau_e min pT %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step09_%dj", jet), Form("mutau_e max deltaPhi e, met %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_step10_%dj", jet), Form("mutau_e min deltaPhi e, mu %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_highmass_%dj", jet), Form("mutau_e high mass %d jet", jet), mutaue_this);
        add_hist_shorthand(&plots_mutaue, Form("mutau_e_lowmass_%dj", jet), Form("mutau_e low mass %d jet", jet), mutaue_this);

        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step03_%dj", jet), Form("etau_mu %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step04_%dj", jet), Form("etau_mu 1+ muon %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step05_%dj", jet), Form("etau_mu 1 muon %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step06_%dj", jet), Form("etau_mu 1+ electron %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step07_%dj", jet), Form("etau_mu 1 electron %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step08_%dj", jet), Form("etau_mu min pT %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step09_%dj", jet), Form("etau_mu max deltaPhi e, met %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_step10_%dj", jet), Form("etau_mu min deltaPhi e, mu %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_highmass_%dj", jet), Form("etau_mu high mass %d jet", jet), etaumu_this);
        add_hist_shorthand(&plots_etaumu, Form("etau_mu_lowmass_%dj", jet), Form("etau_mu low mass %d jet", jet), etaumu_this);
    }

    bool plotthis_mutaue_inclusive[3];
    bool plotthis_etaumu_inclusive[3];
    
    vector<vector<bool>> plotthis_mutaue;
    vector<vector<bool>> plotthis_etaumu;
    for (int jet=0; jet<=MAX_JETS; jet++)
    {
        vector<bool> p1, p2;
        for (int s=0; s<10; s++) 
        {
            p1.push_back(false);
            p2.push_back(false);
        }
        plotthis_mutaue.push_back(p1);
        plotthis_etaumu.push_back(p2);
    }

    double mass_collinear_mutaue = 0;
    double mass_collinear_etaumu = 0;
    plots_mutaue.PrimeFill(&mass_collinear_mutaue);
    plots_etaumu.PrimeFill(&mass_collinear_etaumu);

    bool full_calculate;
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

        plotthis_mutaue_inclusive[0] = true;
        plotthis_mutaue_inclusive[1] = false;
        plotthis_mutaue_inclusive[2] = false;
        for (int j=0; j<=MAX_JETS; j++)
        {
            for (int s=0; s<10; s++) plotthis_mutaue[j][s] = false;
        }

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
        //plotthis_mutaue[0] = passed_b_jets.size() == 0;
        //plotthis_mutaue[1] = plotthis_mutaue[0] and passed_jets.size() <= 1;

        plotthis_mutaue_inclusive[1] = passed_b_jets.size() == 0;
        plotthis_mutaue_inclusive[2] = plotthis_mutaue_inclusive[1] and passed_jets.size() <= MAX_JETS;

        /*
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
        if (plotthis_mutaue[16]) // zero jets
        {
            plotthis_mutaue[18] = indelphes->Muon_PT[only_mu] > 150 and deltaPhi_e_met < 0.3;
            plotthis_mutaue[20] = indelphes->Muon_PT[only_mu] > 60  and deltaPhi_e_met < 0.7;
        }
        if (plotthis_mutaue[17]) // one jet
        {
            plotthis_mutaue[19] = indelphes->Muon_PT[only_mu] > 150 and deltaPhi_e_met < 0.3;
            plotthis_mutaue[21] = indelphes->Muon_PT[only_mu] > 60  and deltaPhi_e_met < 0.7;
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

        for (int i=0; i<23; i++) if (plotthis_mutaue[i]) plots_mutaue.Fill(i);
        */

        for (int njet=0; njet<=MAX_JETS; njet++)
        {
            if (!plotthis_mutaue_inclusive[2]) continue;
            plotthis_mutaue[njet][0] = passed_jets.size() == njet;
            if (plotthis_mutaue[njet][0])
            {
                muon_vec = find_mu(indelphes, 53, -1);
                plotthis_mutaue[njet][1] = muon_vec.size() > 0;
                plotthis_mutaue[njet][2] = muon_vec.size() == 1;
            }
            if (plotthis_mutaue[njet][2])
            {
                electron_vec = find_ele(indelphes, 10, plotthis_mutaue[njet][2] ? muon_vec[0] : -1);
                plotthis_mutaue[njet][3] = electron_vec.size() > 0;
                plotthis_mutaue[njet][4] = electron_vec.size() == 1;
            }
            if (plotthis_mutaue[njet][4])
            {
                only_mu = muon_vec[0];
                only_ele = electron_vec[0];
                plotthis_mutaue[njet][5] = indelphes->Muon_PT[only_mu] > 60;
            }
            if (plotthis_mutaue[njet][5])
            {
                deltaPhi_e_met = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->MissingET_Phi[0]);
                plotthis_mutaue[njet][6] = deltaPhi_e_met < 0.7;
            }
            if (plotthis_mutaue[njet][6])
            {
                deltaPhi_e_mu = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->Muon_Phi[only_mu]);
                plotthis_mutaue[njet][7] = deltaPhi_e_mu > 2.2;
            }
            if (plotthis_mutaue[njet][7])
            {
                plotthis_mutaue[njet][8] = indelphes->Muon_PT[only_mu] > 150 and deltaPhi_e_met < 0.3;
                plotthis_mutaue[njet][9] = indelphes->Muon_PT[only_mu] > 60 and deltaPhi_e_met < 0.7;
            }
        }

        full_calculate = false;
        for (int njet=0; njet<=MAX_JETS; njet++) full_calculate = full_calculate or plotthis_mutaue[njet][7];

        if (full_calculate)
        {
            p4_tau.SetPtEtaPhiM(indelphes->Electron_PT[only_ele], indelphes->Electron_Eta[only_ele], indelphes->Electron_Phi[only_ele], 0.000511);
            p4_lepton.SetPtEtaPhiM(indelphes->Muon_PT[only_mu], indelphes->Muon_Eta[only_mu], indelphes->Muon_Phi[only_mu], 0.10566);
        }
        else
        {
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

        for (int i=0; i<3; i++) if (plotthis_mutaue_inclusive[i]) plots_mutaue.Fill(histogram_numbers_mutaue_inclusive[i]);
        for (int j=0; j<=MAX_JETS; j++)
        {
            for (int i=0; i<10; i++) if (plotthis_mutaue[j][i]) plots_mutaue.Fill(histogram_numbers_mutaue[j][i]);
        }

        ///////////////////////////////////////
        // e + tau_mu
        ///////////////////////////////////////

        plotthis_etaumu_inclusive[0] = true;
        plotthis_etaumu_inclusive[1] = false;
        plotthis_etaumu_inclusive[2] = false;
        for (int j=0; j<=MAX_JETS; j++)
        {
            for (int s=0; s<10; s++) plotthis_etaumu[j][s] = false;
        }

        passed_jets.clear();
        passed_b_jets.clear();
        muon_vec.clear();
        electron_vec.clear();
        only_mu  = -1;
        only_ele = -1;

        /*
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
        if (plotthis_etaumu[16]) // zero jets
        {
            plotthis_etaumu[18] = indelphes->Electron_PT[only_ele] > 150 and deltaPhi_mu_met < 0.3;
            plotthis_etaumu[20] = indelphes->Electron_PT[only_ele] > 60  and deltaPhi_mu_met < 0.7;
        }
        if (plotthis_etaumu[17]) // one jet
        {
            plotthis_etaumu[19] = indelphes->Electron_PT[only_ele] > 150 and deltaPhi_mu_met < 0.3;
            plotthis_etaumu[21] = indelphes->Electron_PT[only_ele] > 60  and deltaPhi_mu_met < 0.7;
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

        for (int i=0; i<23; i++) if (plotthis_etaumu[i]) plots_etaumu.Fill(i);
        */

        for (int j=0; j<indelphes->Jet_size; j++)
        {
            if (indelphes->Jet_PT[j] < 30) continue;
            if (TMath::Abs(indelphes->Jet_Eta[j]) > 6.0) continue;
            passed_jets.push_back(j);
            if (indelphes->Jet_BTag[j] & 0b111) passed_b_jets.push_back(j);
        }
        //if (passed_b_jets.size() == 0) plotthis_etaumu[0] = true;
        plotthis_etaumu_inclusive[1] = passed_b_jets.size() == 0;
        plotthis_etaumu_inclusive[2] = plotthis_etaumu_inclusive[1] and passed_jets.size() <= MAX_JETS;

        for (int njet=0; njet<=MAX_JETS; njet++)
        {
            if (!plotthis_etaumu_inclusive[2]) continue;
            plotthis_etaumu[njet][0] = passed_jets.size() == njet;
            if (plotthis_etaumu[njet][0])
            {
                electron_vec = find_ele(indelphes, 26, -1);
                plotthis_etaumu[njet][1] = electron_vec.size() > 0;
                plotthis_etaumu[njet][2] = electron_vec.size() == 1;
            }
            if (plotthis_etaumu[njet][2])
            {
                muon_vec = find_mu(indelphes, 10, plotthis_etaumu[njet][6] ? electron_vec[0] : -1);
                plotthis_etaumu[njet][3] = muon_vec.size() > 0;
                plotthis_etaumu[njet][4] = muon_vec.size() == 1;
            }
            if (plotthis_etaumu[njet][4])
            {
                only_mu = muon_vec[0];
                only_ele = electron_vec[0];
                plotthis_etaumu[njet][5] = indelphes->Electron_PT[only_ele] > 60;
            }
            if (plotthis_etaumu[njet][5])
            {
                deltaPhi_mu_met = deltaPhi(indelphes->Muon_Phi[only_mu], indelphes->MissingET_Phi[0]);
                plotthis_etaumu[njet][6] = deltaPhi_mu_met < 0.7;
            }
            if (plotthis_etaumu[njet][6])
            {
                deltaPhi_e_mu = deltaPhi(indelphes->Electron_Phi[only_ele], indelphes->Muon_Phi[only_mu]);
                plotthis_etaumu[njet][7] = deltaPhi_mu_met > 2.2;
            }
            if (plotthis_etaumu[njet][7])
            {
                plotthis_etaumu[njet][8] = indelphes->Electron_PT[only_ele] > 150 and deltaPhi_mu_met < 0.3;
                plotthis_etaumu[njet][9] = indelphes->Electron_PT[only_ele] > 60 and deltaPhi_mu_met < 0.7;
            }
        }

        full_calculate = false;
        for (int njet=0; njet<=MAX_JETS; njet++) full_calculate = full_calculate or plotthis_etaumu[njet][7];

        if (full_calculate)
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

        for (int i=0; i<3; i++) if (plotthis_etaumu_inclusive[i]) plots_etaumu.Fill(histogram_numbers_etaumu_inclusive[i]);
        for (int j=0; j<=MAX_JETS; j++)
        {
            for (int i=0; i<10; i++) if (plotthis_etaumu[j][i]) plots_etaumu.Fill(histogram_numbers_etaumu[j][i]);
        }
    }

    TFile *outfile = new TFile(outfilename, "RECREATE");
    plots_mutaue.SaveAll(outfile);
    plots_etaumu.SaveAll(outfile);
    outfile->Close();
}