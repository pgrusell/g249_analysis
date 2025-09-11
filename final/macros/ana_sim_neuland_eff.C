#include "/nucl_lustre/martinaff/s509/Antoine/Reduced_Data_study/mass.hh"

void ana_sim_neuland_eff(const Float_t nev = -1)
{

    gStyle->SetOptStat(0);
    const UInt_t NRGBs = 5;
    Double_t stops[NRGBs] = {0.00, 0.12, 0.32, 0.55, 1.0};
    Double_t red[NRGBs] = {87. / 255, 109. / 255, 181. / 255, 229. / 255, 242. / 255};
    Double_t green[NRGBs] = {100. / 255, 89. / 255, 101. / 255, 107. / 255, 216. / 255};
    Double_t blue[NRGBs] = {161. / 255, 122. / 255, 118. / 255, 111. / 255, 143. / 255};
    Int_t NCont = 60;
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    Int_t ci0[NRGBs];
    TColor *color0[NRGBs];

    for (int i = 0; i < 5; i++)
    {
        ci0[i] = TColor::GetFreeColorIndex();
        color0[i] = new TColor(ci0[i], red[i], green[i], blue[i]);
    }

    // for lines (1D histos and TGraphs)

    const int n = 6;
    Double_t reds[n] = {21. / 255, 111. / 255, 229. / 255, 242. / 255, 119. / 255, 198. / 255};
    Double_t greens[n] = {102. / 255, 52. / 255, 107. / 255, 216. / 255, 163. / 255, 160. / 255};
    Double_t blues[n] = {105. / 255, 79. / 255, 111. / 255, 143. / 255, 80. / 255, 56. / 255};
    Int_t color[n];

    for (int i = 0; i < n; i++)
    {
        color[i] = TColor::GetColor((Float_t)reds[i], (Float_t)greens[i], (Float_t)blues[i]);
    }

    TStopwatch timer;
    timer.Start();

    // TString filename = "/nucl_lustre/s509/sim/sim_c17.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17.geo.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17_withglad_inclxx_202508_13dp.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_cf17_withglad_inclxx_202508.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17_withglad_inclxx_hp_202508_13dp.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17_withglad_qgsp_inclxx_202508_13dp.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17_withglad_qgsp_inclxx_hp_emv_202508_13dp.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17_withoutglad_qgsp_inclxx_hp_emv_202508_13dp.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17_withoutglad_hitdata_withgeo_qgsp_inclxx_hp_emv_202508_13dp.root";
    //  TString filename = "/nucl_lustre/s509/sim/sim_c17_withglad_qgsp_inclxx_hp_emv_202508_13dp.root";
    // TString filename = "/nucl_lustre/s509/sim/sim_c17_withglad_qgsp_inclxx_hp_emv_20250906_13dp.root";
    TString filename = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/sim/sim_eff_study.root";

    TFile *file = TFile::Open(filename);
    TTree *tree = (TTree *)file->Get("evt");
    Int_t nevents = tree->GetEntries();
    if (nev > -1)
        nevents = nev;
    std::cout << "Events: " << nevents << std::endl;

    TH1F *h1_erel = new TH1F("h1_erel", "h1_erel", 300, 0, 15);
    TH1F *h1_erel_neuland = new TH1F("h1_erel_neuland", "h1_erel_neuland", 300, 0, 15);
    TH1F *h1_erel_neuland_p = new TH1F("h1_erel_neuland_p", "h1_erel_neuland_p", 300, 0, 15);
    TH1F *h1_erel_neuland_plane1 = new TH1F("h1_erel_neuland_plane1", "h1_erel_neuland_plane1", 300, 0, 15);
    TH2F *h2_erel_theta = new TH2F("h2_erel_theta", "h2_erel_theta", 300, 0, 15, 300, 0, 9);
    TH1F *h1_neulandZ = new TH1F("h1_neulandZ", "h1_neulandZ", 1000, 1550, 1750);
    TGraph *g_eff = new TGraph();
    TGraph *g_eff_p = new TGraph();
    TGraph *g_eff_plane1 = new TGraph();

    const int res_n_points = 2;
    float res_points[res_n_points] = {0.045, 0.6};
    TH1F *h1_erel_res[res_n_points];
    for (int i = 0; i < res_n_points; i++)
    {
        h1_erel_res[i] = new TH1F(Form("h1_erel_res_%d", i), Form("h1_erel_res_%d", i), 300, 0, 5);
    }
    // R3BMCTrack** track;
    TClonesArray *MCTrackCA = new TClonesArray("R3BMCTrack");
    TBranch *branchMCTrack = tree->GetBranch("MCTrack");
    branchMCTrack->SetAddress(&MCTrackCA);
    TClonesArray *neulandData = new TClonesArray("R3BNeulandPoint");
    TBranch *branchNeulandData = tree->GetBranch("NeulandPoints");
    branchNeulandData->SetAddress(&neulandData);

    TLorentzVector lvbeam;
    // double KE_beam = 530 * 19;
    double KE_beam = 650 * 23;
    const double UNIT = 931.4940954;

    int Z_frag = 8;
    int A_frag = 22;
    Int_t pdgfrag = 10000 * (Z_frag + 100000) + A_frag * 10;
    float mfrag = Nuke_Mass_Tab[22][8];
    float m_neut = 939.565379;
    const double Mbeam = Nuke_Mass_Tab[23][8];
    double E_beam = KE_beam + Mbeam;
    double P_beam = sqrt(E_beam * E_beam - Mbeam * Mbeam);
    lvbeam.SetPxPyPzE(0., 0., P_beam, E_beam); // 4-momentum beam
    TLorentzVector lvtarget;
    TLorentzVector lvstart;
    const Double_t Mtarget = 938.27208816;
    lvtarget.SetPxPyPzE(0, 0, 0, Mtarget); // 4-momentum target
    lvstart = lvbeam + lvtarget;

    float gamma_neu, gamma_frag;
    float beta_neu, beta_frag;
    bool good_frag = false;
    int initfrag = 0;
    TLorentzVector lvfrag;
    TLorentzVector lvneutron;
    TLorentzVector lvneuland;
    TVector3 pfrag, pneutron, pneuland;

    double px, py, pz, p_2, etot;
    float theta_neuland_frag, theta_n_frag;
    float Erel;
    float gamma_neuland, beta_neuland;
    int Z_track, A_track;
    float x_point, y_point, z_point;
    float ang_xz, ang_yz;
    int count_n = 0;
    int count_n_geo = 0;

    float mean_tof = 79.262e-9;
    float s_tof = 150e-12; // 150 ps (from neuland nimA)
    float dist = 16.07;    // meters from target to center of neuland
    float s_beta = s_tof * dist / (mean_tof * mean_tof * TMath::C());

    float rdm_beta, rdm_gamma, rdm_Erel;
    float Erel_tolerance = 0.002;

    for (int ie = 0; ie < nevents; ie++)
    {
        if (ie % 100 == 0)
        {
            printf("Processed: %d of %d (%.2f of 100) \r", ie, nevents, 100. * (ie) / nevents);
            fflush(stdout);
        }

        good_frag = false;
        MCTrackCA->Clear();
        neulandData->Clear();
        lvfrag.Clear();
        lvneutron.Clear();
        pfrag.Clear();
        pneutron.Clear();
        lvneuland.Clear();

        gamma_neu = 0;
        gamma_frag = 0;
        beta_frag = 0;
        beta_neu = 0;
        Erel = 0;

        tree->GetEntry(ie);

        Int_t MCtracksPerEvent = MCTrackCA->GetEntriesFast();
        if (MCtracksPerEvent > 0)
        {
            R3BMCTrack **track = new R3BMCTrack *[MCtracksPerEvent];
            for (Int_t j = 0; j < MCtracksPerEvent; j++)
            {
                // check for good final fragment
                track[j] = (R3BMCTrack *)MCTrackCA->At(j);
                if (track[j]->GetMotherId() == -1)
                {

                    auto mom = track[j]->GetMomentumMass() * 1000.;
                    if (track[j]->GetPdgCode() > 1000000)
                    {

                        // Charge and mass are obtained from PDG Code
                        Z_track = int(track[j]->GetPdgCode() / 10000) - 100000.;
                        A_track = 0.1 * (track[j]->GetPdgCode() - (100000 + Z_track) * 10000.);

                        if (Z_track == Z_frag && A_track == A_frag)
                        {
                            // std::cout << "hola" << std::endl;

                            good_frag = true;
                            initfrag++;
                            mom.SetPx(mom.X() + gRandom->Gaus(0, 0.005) * mom.X());
                            mom.SetPy(mom.Y() + gRandom->Gaus(0, 0.005) * mom.Y());
                            mom.SetPz(mom.Z() + gRandom->Gaus(0, 0.005) * mom.Z());
                            double new_E = sqrt(mom.X() * mom.X() + mom.Y() * mom.Y() + mom.Z() * mom.Z() + (track[j]->GetMass() * track[j]->GetMass()) * 1e6);
                            mom.SetPxPyPzE(mom.X(), mom.Y(), mom.Z(), new_E);
                            lvfrag.SetPxPyPzE(mom.X(), mom.Y(), mom.Z(), mom.E());
                            pfrag.SetXYZ(mom.X(), mom.Y(), mom.Z());
                            gamma_frag = mom.E() / mfrag;
                            beta_frag = mom.Beta();
                        }
                    }

                    if (track[j]->GetPdgCode() == 2112)
                    {
                        auto mom = track[j]->GetMomentumMass() * 1000.;
                        lvneutron.SetPxPyPzE(mom.X(), mom.Y(), mom.Z(), mom.E());
                        pneutron.SetXYZ(mom.X(), mom.Y(), mom.Z());
                        if (pneutron.Theta() * TMath::RadToDeg() < 15)
                        {
                            count_n++;
                            gamma_neu = mom.E() / m_neut;
                            beta_neu = mom.Beta();
                            rdm_beta = beta_neu + gRandom->Gaus(0, s_beta);
                            rdm_gamma = 1. / sqrt(1. - rdm_beta * rdm_beta);
                            ang_xz = atan(pneutron.X() / pneutron.Z());
                            ang_yz = atan(pneutron.Y() / pneutron.Z());
                            x_point = tan(ang_xz) * 1590;
                            y_point = tan(ang_yz) * 1590;
                            if (abs(x_point) < 125 && abs(y_point) < 125)
                                count_n_geo++;
                        }
                    }
                }
            }

            if (track)
                delete[] track;

            theta_n_frag = pneutron.Angle(pfrag);
            Erel = sqrt(mfrag * mfrag + m_neut * m_neut + 2 * gamma_neu * gamma_frag * mfrag * m_neut * (1 - beta_neu * beta_frag * TMath::Cos(theta_n_frag))) - mfrag - m_neut;
            rdm_Erel = sqrt(mfrag * mfrag + m_neut * m_neut + 2 * rdm_gamma * gamma_frag * mfrag * m_neut * (1 - rdm_beta * beta_frag * TMath::Cos(theta_n_frag))) - mfrag - m_neut;
            // std::cout << Erel << " " << rdm_Erel << std::endl;
            if (Erel > 0)
            {
                h1_erel->Fill(Erel);
                h2_erel_theta->Fill(Erel, pneutron.Theta() * TMath::RadToDeg());
            }

            for (int i = 0; i < res_n_points; i++)
            {
                if (abs(Erel - res_points[i]) < Erel_tolerance)
                {
                    h1_erel_res[i]->Fill(rdm_Erel);
                }
            }
        }

        bool neuland_p = false;
        bool first = true;
        std::vector<int> neuland_pdgs;
        std::vector<int> neuland_energies;
        Int_t neulandPerEvent = neulandData->GetEntries();
        if (neulandPerEvent > 0 && Erel > 0)
            h1_erel_neuland->Fill(Erel);
        R3BNeulandPoint **neulandPointData = new R3BNeulandPoint *[neulandPerEvent];
        for (Int_t j = 0; j < neulandPerEvent; j++)
        {
            neulandPointData[j] = (R3BNeulandPoint *)(neulandData->At(j));
            Int_t TrackId = neulandPointData[j]->GetTrackID();
            R3BMCTrack *Track = (R3BMCTrack *)MCTrackCA->At(TrackId);
            Int_t MotherId = Track->GetMotherId();
            Int_t pdg_neuland = Track->GetPdgCode();

            if (abs(neulandPointData[j]->GetZ() - 1590) < 2.5 && first)
            {
                first = false;
                h1_erel_neuland_plane1->Fill(Erel);
            }

            // neuland_pdgs.push_back(pdg_neuland);
            float px = neulandPointData[j]->GetPx() * 1000;
            float py = neulandPointData[j]->GetPy() * 1000;
            float pz = neulandPointData[j]->GetPz() * 1000;
            // neuland_energies.push_back(sqrt(px*px+py*py+pz*pz));
            // if (MotherId==-1 && pdg_neuland == 2112)
            if (pdg_neuland == 2212 && sqrt(px * px + py * py + pz * pz) > 60) //&& neulandPointData[j]->GetEnergyLoss()>0.001)
            {
                neuland_p = true;
                h1_neulandZ->Fill(neulandPointData[j]->GetZ());
                // px = neulandPointData[j]->GetPx()*1000;
                // py = neulandPointData[j]->GetPy()*1000;
                // pz = neulandPointData[j]->GetPz()*1000;
                // p_2 = px*px+py*py+pz*pz;
                // pneuland.SetXYZ(px,py,pz);
                // etot = sqrt(p_2 + m_neut*m_neut);
                // lvneuland.SetPxPyPzE(px,py,pz,etot);
                // theta_neuland_frag = pneuland.Angle(pfrag);
                // gamma_neuland = lvneuland.E()/m_neut;
                // beta_neuland = lvneuland.Beta();
                // Erel = sqrt(mfrag*mfrag + m_neut*m_neut + 2*gamma_neuland*gamma_frag*mfrag*m_neut*(1-beta_neuland*beta_frag*TMath::Cos(theta_neuland_frag))) - mfrag - m_neut;
                // //std::cout << Erel << std::endl;
                // if (Erel>0)
                // {
                //     h1_erel_neuland->Fill(Erel);
                // }
            }
        }
        if (neuland_p && Erel)
            h1_erel_neuland_p->Fill(Erel);
        // if (neuland_p)
        // {
        //   //std::cout << std::endl;
        //   for (int i=0;i<neuland_pdgs.size();i++)
        //   { std::cout << neuland_pdgs[i] << " " <<  neuland_energies[i] << std::endl; }
        // }

        // if (!good_fragment_found) continue;
    }

    TCanvas *c_erel = new TCanvas("c_erel", "c_erel", 700, 500);
    c_erel->SetMargin(0.15, 0.03, 0.15, 0.02);
    h1_erel->SetLineColor(color[2]);
    h1_erel->SetLineWidth(2);
    h1_erel->SetTitle("");
    h1_erel->SetLineStyle(3);
    h1_erel->GetXaxis()->SetRangeUser(0.001, 12);
    // h1_erel->GetYaxis()->SetRangeUser(0.001,1300);
    h1_erel->GetXaxis()->SetTickLength(0.02);
    h1_erel->GetYaxis()->SetTickLength(0.02);
    h1_erel->GetXaxis()->SetNdivisions(505);
    h1_erel->GetYaxis()->SetNdivisions(505);
    h1_erel->GetXaxis()->SetTitle("E_{rel} [MeV]");
    h1_erel->GetYaxis()->SetTitle("Counts");
    h1_erel->GetXaxis()->SetTitleSize(0.07);
    h1_erel->GetYaxis()->SetTitleSize(0.07);
    h1_erel->GetXaxis()->SetTitleOffset(0.85);
    h1_erel->GetYaxis()->SetTitleOffset(1.05);
    h1_erel->GetXaxis()->SetLabelSize(0.06);
    h1_erel->GetYaxis()->SetLabelSize(0.06);
    h1_erel->GetXaxis()->SetNdivisions(209);
    h1_erel->GetYaxis()->SetNdivisions(205);
    h1_erel->GetXaxis()->CenterTitle();
    h1_erel->GetYaxis()->CenterTitle();
    h1_erel->Draw();
    h1_erel_neuland->SetLineColor(color[1]);
    h1_erel_neuland->SetLineWidth(2);
    h1_erel_neuland->SetLineStyle(1);
    h1_erel_neuland->Draw("same");
    h1_erel_neuland_p->SetLineColor(color[0]);
    h1_erel_neuland_p->SetLineStyle(9);
    h1_erel_neuland_p->SetLineWidth(2);
    h1_erel_neuland_p->Draw("same");
    h1_erel->SetName("MCTrack");
    h1_erel_neuland->SetName("Neuland");
    h1_erel_neuland->SetName("Neuland proton");
    h1_erel_neuland_plane1->SetLineColor(kOrange);
    // h1_erel_neuland_plane1->Draw("same");
    TLatex latex1;
    latex1.SetTextSize(0.05);
    latex1.SetTextFont(42);
    // latex1.DrawLatex(-0.2, -0.2, "0");
    TLegend *leg = new TLegend(0.5, 0.5, 0.85, 0.8);
    leg->AddEntry(h1_erel, "MCTrack", "l");
    leg->AddEntry(h1_erel_neuland, "Neuland all", "l");
    leg->AddEntry(h1_erel_neuland_p, "Neuland proton", "l");
    // leg->AddEntry(h1_erel_neuland_plane1,"Neuland plane 1","l");
    leg->SetBorderSize(0);
    leg->Draw();
    c_erel->SaveAs("c_erel_eff_eval_withglad.pdf");
    h1_erel->SetName("h1_erel_incl");
    h1_erel->SaveAs("h1_erel_incl_withglad.C");

    TH1F *h1_eff = new TH1F("h1_eff", "h1_eff", 300, 0, 15);
    for (int i = 0; i < h1_erel->GetNbinsX(); i++)
    {
        float eff = h1_erel_neuland->GetBinContent(i + 1) / h1_erel->GetBinContent(i + 1);
        h1_eff->SetBinContent(i + 1, eff);
        g_eff->SetPoint(i, h1_erel->GetBinCenter(i + 1), eff * 100);
        g_eff_p->SetPoint(i, h1_erel->GetBinCenter(i + 1), h1_erel_neuland_p->GetBinContent(i + 1) / h1_erel->GetBinContent(i + 1) * 100);
        g_eff_plane1->SetPoint(i, h1_erel->GetBinCenter(i + 1), h1_erel_neuland_plane1->GetBinContent(i + 1) / h1_erel->GetBinContent(i + 1));
    }

    TCanvas *c_eff = new TCanvas("c_eff", "c_eff", 700, 500);
    c_eff->SetMargin(0.15, 0.03, 0.15, 0.02);
    // h1_eff->Draw();
    g_eff->SetMarkerStyle(20);
    g_eff->SetMarkerSize(1);
    g_eff->SetMarkerColor(color[1]);
    g_eff->SetLineColor(color[1]);
    g_eff->SetLineWidth(2);
    g_eff->SetLineStyle(2);
    g_eff->GetXaxis()->SetTickLength(0.02);
    g_eff->GetYaxis()->SetTickLength(0.02);
    g_eff->GetXaxis()->SetTitle("E_{rel} [MeV]");
    g_eff->GetYaxis()->SetTitle("Efficiency (%)");
    g_eff->GetXaxis()->SetTitleSize(0.07);
    g_eff->GetYaxis()->SetTitleSize(0.07);
    g_eff->GetXaxis()->SetTitleOffset(0.85);
    g_eff->GetYaxis()->SetTitleOffset(0.8);
    g_eff->GetXaxis()->SetLabelSize(0.06);
    g_eff->GetYaxis()->SetLabelSize(0.06);
    g_eff->GetXaxis()->SetNdivisions(206);
    g_eff->GetYaxis()->SetNdivisions(206);
    g_eff->GetXaxis()->CenterTitle();
    g_eff->GetYaxis()->CenterTitle();
    g_eff->GetXaxis()->SetRangeUser(0.05, 9.8);
    g_eff->GetYaxis()->SetRangeUser(0.001 * 100, 1.05 * 100);
    TLatex latex;
    latex.SetTextSize(0.05);
    latex.SetTextFont(42);
    // latex.DrawLatex(-0.25, -0.25, "0");
    g_eff->Draw("AL");
    g_eff_p->SetMarkerStyle(21);
    g_eff_p->SetMarkerSize(1);
    g_eff_p->SetMarkerColor(color[0]);
    g_eff_p->SetLineColor(color[0]);
    g_eff_p->SetLineWidth(2);
    g_eff_p->SetLineStyle(1);
    g_eff_p->Draw("LSAME");
    g_eff_plane1->SetMarkerStyle(21);
    g_eff_plane1->SetMarkerSize(1.5);
    g_eff_plane1->SetMarkerColor(kOrange);
    TLine *leff = new TLine(0, float(count_n_geo) / count_n, 12, float(count_n_geo) / count_n);
    leff->SetLineStyle(kDashed);
    leff->SetLineWidth(4);
    leff->SetLineColor(color[2]);
    // leff->Draw("same l");
    TLegend *leg_eff = new TLegend(0.16, 0.18, 0.6, 0.45);
    leg_eff->AddEntry(g_eff, "Neuland all", "l");
    leg_eff->AddEntry(g_eff_p, "Neuland proton", "l");
    // leg_eff->AddEntry(leff,"<#epsilon_{geo}>","l");
    leg_eff->SetBorderSize(0);
    leg_eff->Draw("same");
    // g_eff_plane1->Draw("PSAME");
    c_eff->SaveAs("c_eff_withglad.pdf");
    c_eff->SaveAs("c_eff_withglad.C");

    TCanvas *c_theta_erel = new TCanvas("c_theta_erel", "c_theta_erel", 700, 500);
    c_theta_erel->SetMargin(0.15, 0.03, 0.15, 0.02);
    h2_erel_theta->SetTitle("");
    h2_erel_theta->GetXaxis()->SetRangeUser(0.0, 12);
    h2_erel_theta->GetYaxis()->SetRangeUser(0.05, 9);
    h2_erel_theta->GetXaxis()->SetTickLength(0.02);
    h2_erel_theta->GetYaxis()->SetTickLength(0.02);
    h2_erel_theta->GetXaxis()->SetTitle("E_{rel} [MeV]");
    h2_erel_theta->GetYaxis()->SetTitle("polar angle #theta [deg]");
    h2_erel_theta->GetXaxis()->SetTitleSize(0.07);
    h2_erel_theta->GetYaxis()->SetTitleSize(0.07);
    h2_erel_theta->GetXaxis()->SetTitleOffset(0.85);
    h2_erel_theta->GetYaxis()->SetTitleOffset(0.8);
    h2_erel_theta->GetXaxis()->SetLabelSize(0.06);
    h2_erel_theta->GetYaxis()->SetLabelSize(0.06);
    h2_erel_theta->GetXaxis()->SetNdivisions(208);
    h2_erel_theta->GetYaxis()->SetNdivisions(206);
    h2_erel_theta->GetXaxis()->CenterTitle();
    h2_erel_theta->GetYaxis()->CenterTitle();
    h2_erel_theta->Draw("col");
    TLatex latex3;
    latex3.SetTextSize(0.05);
    latex3.SetTextFont(42);
    // latex3.DrawLatex(-0.25, -0.25, "0");
    TLine *l = new TLine(0, 4.63, 12, 4.63);
    l->SetLineStyle(kDashed);
    l->SetLineWidth(4);
    l->SetLineColor(color[0]);
    l->Draw("same l");
    // c_theta_erel->cd(2);
    // h1_neulandZ->Draw();
    c_theta_erel->SaveAs("c_theta_erel_withglad.pdf");

    TCanvas *c_erel_res = new TCanvas("c_erel_res", "c_erel_res", 700, 500);
    c_erel_res->SetMargin(0.15, 0.03, 0.15, 0.02);
    for (int i = 0; i < res_n_points; i++)
    {
        h1_erel_res[i]->SetLineColor(color[i]);
        h1_erel_res[i]->SetLineWidth(2);
        h1_erel_res[i]->SetTitle("");
        h1_erel_res[i]->GetXaxis()->SetTickLength(0.02);
        h1_erel_res[i]->GetYaxis()->SetTickLength(0.02);
        h1_erel_res[i]->GetXaxis()->SetTitleSize(0.07);
        h1_erel_res[i]->GetYaxis()->SetTitleSize(0.07);
        h1_erel_res[i]->GetXaxis()->SetTitleOffset(0.85);
        h1_erel_res[i]->GetYaxis()->SetTitleOffset(0.85);
        h1_erel_res[i]->GetXaxis()->SetLabelSize(0.06);
        h1_erel_res[i]->GetYaxis()->SetLabelSize(0.06);
        h1_erel_res[i]->GetXaxis()->SetNdivisions(206);
        h1_erel_res[i]->SetName(Form("E_{rel} = %.1f MeV", res_points[i]));

        if (i == 0)
        {
            h1_erel_res[i]->Draw("hist");
        }
        else
        {
            h1_erel_res[i]->Draw("same");
        }

        h1_erel_res[i]->Fit("gaus", "Q0");
        h1_erel_res[i]->GetFunction("gaus")->SetLineColor(color[i]);
        h1_erel_res[i]->GetFunction("gaus")->SetLineStyle(2);
        h1_erel_res[i]->GetFunction("gaus")->Draw("same");
    }
    h1_erel_res[0]->GetXaxis()->SetTitle("E_{rel} [MeV]");
    h1_erel_res[0]->GetYaxis()->SetTitle("Counts");
    h1_erel_res[0]->GetXaxis()->CenterTitle();
    h1_erel_res[0]->GetYaxis()->CenterTitle();

    TCanvas *c_res_extra = new TCanvas("c_res_extra", "c_res_extra", 700, 500);
    c_res_extra->SetMargin(0.15, 0.03, 0.15, 0.02);
    TGraphErrors *g_res = new TGraphErrors();
    TGraphErrors *g_sigma_o_mean = new TGraphErrors();
    float sigma, mean, ratio_err;
    for (int i = 0; i < res_n_points; i++)
    {
        mean = h1_erel_res[i]->GetFunction("gaus")->GetParameter(1);
        sigma = h1_erel_res[i]->GetFunction("gaus")->GetParameter(2);
        g_res->SetPoint(i, res_points[i], sigma);
        g_res->SetPointError(i, 0, h1_erel_res[i]->GetFunction("gaus")->GetParError(2));
        g_sigma_o_mean->SetPoint(i, res_points[i], sigma / mean);
        ratio_err = sqrt(pow(sigma / mean / mean, 2) * pow(h1_erel_res[i]->GetFunction("gaus")->GetParError(1), 2) + pow(mean, -2) * pow(h1_erel_res[i]->GetFunction("gaus")->GetParError(2), 2));
        g_sigma_o_mean->SetPointError(i, 0, ratio_err);
    }
    g_res->SetMarkerStyle(20);
    g_res->SetMarkerSize(2);
    g_res->SetMarkerColor(color[1]);
    g_res->GetXaxis()->SetTickLength(0.02);
    g_res->GetYaxis()->SetTickLength(0.02);
    g_res->GetXaxis()->SetTitle("E_{rel} [MeV]");
    g_res->GetYaxis()->SetTitle("#sigma_{E_{rel}} [MeV]");
    g_res->GetXaxis()->SetTitleSize(0.07);
    g_res->GetYaxis()->SetTitleSize(0.07);
    g_res->GetXaxis()->SetTitleOffset(0.85);
    g_res->GetYaxis()->SetTitleOffset(0.85);
    g_res->GetXaxis()->SetLabelSize(0.06);
    g_res->GetYaxis()->SetLabelSize(0.06);
    g_res->GetXaxis()->SetNdivisions(206);
    g_res->GetYaxis()->SetNdivisions(206);
    g_res->GetXaxis()->CenterTitle();
    g_res->GetYaxis()->CenterTitle();
    g_res->SetLineColor(color[1]);
    g_res->Draw("AP");
    TF1 *f_res = new TF1("f_res", "[0]*x^[1]", 0, 4.1);
    g_res->Fit("f_res", "Q0");
    f_res->SetLineColor(color[1]);
    f_res->SetLineStyle(2);
    f_res->Draw("same");
    std::cout << "Erel \u03c3 at 1MeV: " << f_res->Eval(1) << std::endl;
    g_sigma_o_mean->SetMarkerStyle(21);
    g_sigma_o_mean->SetMarkerSize(2);
    g_sigma_o_mean->SetMarkerColor(color[0]);
    g_sigma_o_mean->Draw("PSAME");

    std::cout << "Number of neutrons: " << count_n << std::endl;
    std::cout << "Number of neutrons in geometric acceptance: " << count_n_geo << std::endl;
    std::cout << "Geo eff: " << float(count_n_geo) / count_n << std::endl;
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    Float_t cpuUsage = ctime / rtime;
    cout << endl;
    cout << "CPU used: " << cpuUsage << endl;
    cout << "Real time " << rtime / 60. << " min, CPU time " << ctime / 60.
         << "min" << endl;
    cout << "Macro finished successfully." << endl;
}
