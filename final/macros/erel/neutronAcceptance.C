void neutronAcceptance(const char *inName =
                           "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/data_23O.root",
                       const char *outName = "23O_neutron_acceptance.root")
{
    const double m_neut = 0.939565; // GeV

    // offsets (23O1n.txt): fragX, fragY, neuX, neuY, betaMatch
    const double fFragOffX = -0.0026286137;
    const double fFragOffY = -0.016158624;
    const double fNeuOffX = -0.0012187827;
    const double fNeuOffY = -0.018529642;
    const double fBetaMatch = -0.00140415;

    const double fxOff = -fFragOffX, fyOff = -fFragOffY;
    const double dxOff = -fNeuOffX, dyOff = -fNeuOffY;

    TFile *f = TFile::Open(inName, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: no se pudo abrir " << inName << std::endl;
        return;
    }
    TTree *t = dynamic_cast<TTree *>(f->Get("FilterDataTree"));
    if (!t)
    {
        std::cerr << "ERROR: no se encontro FilterDataTree.\n";
        f->Close();
        return;
    }

    Double_t M_frag, beta_frag, beta_neu, p_neu;
    Double_t px_frag, py_frag, pz_frag;
    Double_t px_neu, py_neu, pz_neu;
    Double_t x_neu_hit, y_neu_hit, z_neu_hit;

    t->SetBranchAddress("M_frag", &M_frag);
    t->SetBranchAddress("beta_frag", &beta_frag);
    t->SetBranchAddress("px_frag", &px_frag);
    t->SetBranchAddress("py_frag", &py_frag);
    t->SetBranchAddress("pz_frag", &pz_frag);
    t->SetBranchAddress("beta_neu", &beta_neu);
    t->SetBranchAddress("p_neu", &p_neu);
    t->SetBranchAddress("px_neu", &px_neu);
    t->SetBranchAddress("py_neu", &py_neu);
    t->SetBranchAddress("pz_neu", &pz_neu);
    t->SetBranchAddress("x_neu_hit", &x_neu_hit);
    t->SetBranchAddress("y_neu_hit", &y_neu_hit);
    t->SetBranchAddress("z_neu_hit", &z_neu_hit);

    // --- histogramas 1D (una particula) --------------------------------------
    TH1F *hBetaNeu = new TH1F("hBetaNeu", "#beta_{neu};#beta;counts", 100, 0.3, 0.9);
    TH1F *hPNeu = new TH1F("hPNeu", "p_{neu};p [GeV/c];counts", 100, 0, 5);
    TH1F *hCosAng = new TH1F("hCosAng", "cos(#theta_{nf});cos;counts", 100, 0.9, 1.001);
    TH1F *hOpaLab = new TH1F("hOpaLab", "opening angle n-frag (lab);#theta [rad];counts", 100, 0, 0.2);
    TH1F *hXneu = new TH1F("hXneu", "x_{neu} hit;x [cm];counts", 120, -150, 150);
    TH1F *hYneu = new TH1F("hYneu", "y_{neu} hit;y [cm];counts", 120, -150, 150);

    // --- 2D: posicion de impacto en NeuLAND ----------------------------------
    TH2F *hXY = new TH2F("hXY", "NeuLAND hit;x [cm];y [cm]", 120, -150, 150, 120, -150, 150);

    // --- 2D: Erel vs variable de 1 particula (LO DIAGNOSTICO) ----------------
    TH2F *hErelVsCos = new TH2F("hErelVsCos", "E_{rel} vs cos(#theta_{nf});cos;E_{rel} [MeV]",
                                100, 0.9, 1.001, 100, -1, 10);
    TH2F *hErelVsOpa = new TH2F("hErelVsOpa", "E_{rel} vs opening angle;#theta [rad];E_{rel} [MeV]",
                                100, 0, 0.2, 100, -1, 10);
    TH2F *hErelVsX = new TH2F("hErelVsX", "E_{rel} vs x_{neu};x [cm];E_{rel} [MeV]",
                              120, -150, 150, 100, -1, 10);
    TH2F *hErelVsY = new TH2F("hErelVsY", "E_{rel} vs y_{neu};y [cm];E_{rel} [MeV]",
                              120, -150, 150, 100, -1, 10);

    const Long64_t N = t->GetEntries();
    std::cout << "Eventos: " << N << std::endl;

    for (Long64_t i = 0; i < N; ++i)
    {
        t->GetEntry(i);

        double bf = beta_frag + fBetaMatch;
        double fx = px_frag / pz_frag + fxOff;
        double fy = py_frag / pz_frag + fyOff;
        double dx = px_neu / pz_neu + dxOff;
        double dy = py_neu / pz_neu + dyOff;

        double cos_ang =
            (dx * fx + dy * fy + 1.0) /
            (std::sqrt(dx * dx + dy * dy + 1.0) * std::sqrt(fx * fx + fy * fy + 1.0));
        double opa = std::acos(std::min(1.0, std::max(-1.0, cos_ang)));

        double g_neu = 1.0 / std::sqrt(1.0 - beta_neu * beta_neu);
        double g_frag = 1.0 / std::sqrt(1.0 - bf * bf);
        double m_f = M_frag;
        double Erel = 1000.0 * (std::sqrt(m_f * m_f + m_neut * m_neut +
                                          2.0 * g_neu * g_frag * m_f * m_neut * (1.0 - beta_neu * bf * cos_ang)) -
                                m_f - m_neut);

        hBetaNeu->Fill(beta_neu);
        hPNeu->Fill(p_neu);
        hCosAng->Fill(cos_ang);
        hOpaLab->Fill(opa);
        hXneu->Fill(x_neu_hit);
        hYneu->Fill(y_neu_hit);
        hXY->Fill(x_neu_hit, y_neu_hit);

        hErelVsCos->Fill(cos_ang, Erel);
        hErelVsOpa->Fill(opa, Erel);
        hErelVsX->Fill(x_neu_hit, Erel);
        hErelVsY->Fill(y_neu_hit, Erel);
    }
    // f->Close();

    // --- dibujar -------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "1-particle distributions", 1200, 800);
    c1->Divide(3, 2);
    c1->cd(1);
    hBetaNeu->Draw("hist");
    c1->cd(2);
    hPNeu->Draw("hist");
    c1->cd(3);
    hOpaLab->Draw("hist");
    c1->cd(4);
    hXneu->Draw("hist");
    c1->cd(5);
    hYneu->Draw("hist");
    c1->cd(6);
    hXY->Draw("colz");

    TCanvas *c2 = new TCanvas("c2", "Erel vs 1-particle vars", 1200, 800);
    c2->Divide(2, 2);
    c2->cd(1);
    hErelVsOpa->Draw("colz");
    c2->cd(2);
    hErelVsCos->Draw("colz");
    c2->cd(3);
    hErelVsX->Draw("colz");
    c2->cd(4);
    hErelVsY->Draw("colz");

    // --- guardar -------------------------------------------------------------
    TFile *fout = new TFile(outName, "RECREATE");
    hBetaNeu->Write();
    hPNeu->Write();
    hCosAng->Write();
    hOpaLab->Write();
    hXneu->Write();
    hYneu->Write();
    hXY->Write();
    hErelVsCos->Write();
    hErelVsOpa->Write();
    hErelVsX->Write();
    hErelVsY->Write();
    fout->Close();

    std::cout << "Guardado en: " << outName << std::endl;
}