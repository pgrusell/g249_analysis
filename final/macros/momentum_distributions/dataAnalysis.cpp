#include "dataAnalysis.hh"

DataAnalysis::DataAnalysis(std::string dataFilePath, std::string outFileName, bool hasNeutrons)
{
    fDataFilePath = dataFilePath;
    fHasNeutrons = hasNeutrons;
    fOutFileName = outFileName;

    fErel = new TH1F("fErel", "Erel", 100, -5, 10);
    fDeltaBeta = new TH1F("fDeltaBeta", "#Delta#beta", 100, -0.3, 0.3);
}

DataAnalysis::~DataAnalysis()
{
    if (fErel)
    {
        delete fErel;
        fErel = nullptr;
    }

    if (fDeltaBeta)
    {
        delete fDeltaBeta;
        fDeltaBeta = nullptr;
    }
}

void DataAnalysis::setOffsetsFromTxt(std::string offFile)
{
    const std::string txtPath =
        std::string(getenv("repopath")) + "/final/settings/" + offFile;

    std::ifstream in(txtPath);
    std::vector<double> offsets;

    std::string line;
    while (std::getline(in, line))
    {
        offsets.push_back(std::atof(line.c_str()));
    }

    fNeutronPOffsets[0] = offsets[2];
    fNeutronPOffsets[1] = offsets[3];
    fFragmentPOffsets[0] = offsets[0];
    fFragmentPOffsets[1] = offsets[1];
    fBetaMatchValue = offsets[4];
}

void DataAnalysis::getData(bool called)
{
    // Open dataFile
    TFile *f = TFile::Open(fDataFilePath.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[anaErelTree] ERROR: Couldn't open " << fDataFilePath << std::endl;
        return;
    }

    // Get TTree
    TTree *t = dynamic_cast<TTree *>(f->Get("FilterDataTree"));
    if (!t)
    {
        std::cerr << "[FilterDataTree] ERROR: Couldn't find TTree \"FilterDataTree\".\n";
        f->Close();
        return;
    }

    Long64_t nentries = t->GetEntries();

    if (called)
        std::cout << "Number of entries: " << nentries << "\n";

    // Read values from branches
    Double_t AoQ_frag;
    Double_t Z_frag_est;
    Double_t A_frag;
    Double_t M_frag;
    Double_t beta_frag;
    Double_t beta_neu;
    Double_t p_frag;
    Double_t p_neu;
    Double_t califa_opa;
    Double_t px_frag, py_frag, pz_frag;
    Double_t px_neu, py_neu, pz_neu;
    Double_t x_neu_hit, y_neu_hit, z_neu_hit;
    Double_t tof;
    Double_t beta_proj;
    Double_t px_proj, py_proj, pz_proj;

    t->SetBranchAddress("AoQ_frag", &AoQ_frag);
    t->SetBranchAddress("Z_frag_est", &Z_frag_est);
    t->SetBranchAddress("A_frag", &A_frag);
    t->SetBranchAddress("M_frag", &M_frag);
    t->SetBranchAddress("beta_frag", &beta_frag);
    t->SetBranchAddress("p_frag", &p_frag);
    t->SetBranchAddress("califa_opa", &califa_opa);
    t->SetBranchAddress("px_frag", &px_frag);
    t->SetBranchAddress("py_frag", &py_frag);
    t->SetBranchAddress("pz_frag", &pz_frag);
    t->SetBranchAddress("beta_proj", &beta_proj);
    t->SetBranchAddress("px_in", &px_proj);
    t->SetBranchAddress("py_in", &py_proj);
    t->SetBranchAddress("pz_in", &pz_proj);

    if (fHasNeutrons)
    {
        t->SetBranchAddress("beta_neu", &beta_neu);
        t->SetBranchAddress("p_neu", &p_neu);
        t->SetBranchAddress("px_neu", &px_neu);
        t->SetBranchAddress("py_neu", &py_neu);
        t->SetBranchAddress("pz_neu", &pz_neu);
        t->SetBranchAddress("x_neu_hit", &x_neu_hit);
        t->SetBranchAddress("y_neu_hit", &y_neu_hit);
        t->SetBranchAddress("z_neu_hit", &z_neu_hit);
        t->SetBranchAddress("tof_neuland", &tof);
    }

    // Variables that will be stored
    std::string outFile = std::string(getenv("repopath")) +
                          "/results/final/" + fOutFileName;

    TFile *fout = nullptr;
    TTree *tout = nullptr;

    if (called)
    {
        fout = new TFile(outFile.c_str(), "RECREATE");
        tout = new TTree("KinTree", "Kinematics from FilterDataTree");
    }

    Double_t Erel;
    Double_t beta_n_out;
    Double_t beta_f_out;
    Double_t delta_beta_out;
    Double_t dx_neu, dy_neu, dz_neu;
    Double_t fx_frag, fy_frag, fz_frag;
    Double_t opa_out;
    Double_t neuland_tof;
    Double_t pL, pT;
    Double_t p_sys, p_sys_CM;
    Double_t px_sys, py_sys, pz_sys;
    Double_t px_sys_CM, py_sys_CM, pz_sys_CM;
    Double_t opa_lab;
    Double_t px_rf_rot, py_rf_rot, pz_rf_rot;

    if (called)
    {
        tout->Branch("AoQ_frag", &AoQ_frag);
        tout->Branch("Z_frag_est", &Z_frag_est);
        tout->Branch("beta_frag", &beta_f_out);
        tout->Branch("fx_frag", &fx_frag);
        tout->Branch("fy_frag", &fy_frag);
        tout->Branch("fz_frag", &fz_frag);
        tout->Branch("px_frag", &px_frag);
        tout->Branch("py_frag", &py_frag);
        tout->Branch("pz_frag", &pz_frag);
        tout->Branch("pL", &pL);
        tout->Branch("pT", &pT);
        tout->Branch("califa_opa", &opa_out);

        tout->Branch("p_sys", &p_sys);
        tout->Branch("p_sys_CM", &p_sys_CM);
        tout->Branch("px_sys", &px_sys);
        tout->Branch("py_sys", &py_sys);
        tout->Branch("pz_sys", &pz_sys);
        tout->Branch("px_sys_CM", &px_sys_CM);
        tout->Branch("py_sys_CM", &py_sys_CM);
        tout->Branch("pz_sys_CM", &pz_sys_CM);

        tout->Branch("px_rf_rot", &px_rf_rot);
        tout->Branch("py_rf_rot", &py_rf_rot);
        tout->Branch("pz_rf_rot", &pz_rf_rot);

        if (fHasNeutrons)
        {
            tout->Branch("Erel", &Erel);
            tout->Branch("beta_neu", &beta_n_out);
            tout->Branch("delta_beta", &delta_beta_out);
            tout->Branch("dx_neu", &dx_neu);
            tout->Branch("dy_neu", &dy_neu);
            tout->Branch("dz_neu", &dz_neu);
            tout->Branch("px_neu", &px_neu);
            tout->Branch("py_neu", &py_neu);
            tout->Branch("pz_neu", &pz_neu);
            tout->Branch("x_neu", &x_neu_hit);
            tout->Branch("y_neu", &y_neu_hit);
            tout->Branch("z_neu", &z_neu_hit);
            tout->Branch("neuland_tof", &neuland_tof);
            tout->Branch("opa_lab", &opa_lab);
        }
    }

    // Loop over events to calculate variables
    for (Long64_t i = 0; i < nentries; ++i)
    {

        t->GetEntry(i);

        opa_out = califa_opa;
        beta_frag += fBetaMatchValue;
        beta_f_out = beta_frag;

        if (fHasNeutrons)
        {

            neuland_tof = tof;

            ///////////// Relative energy spectrum //////////////
            // Center transversal components of the momentum

            Double_t nz_corr = pz_neu;
            Double_t dx_neu_off = -fNeutronPOffsets[0];
            Double_t dy_neu_off = -fNeutronPOffsets[1];
            Double_t fx_frag_off = -fFragmentPOffsets[0];
            Double_t fy_frag_off = -fFragmentPOffsets[1];

            dx_neu = px_neu / nz_corr + dx_neu_off;
            dy_neu = py_neu / nz_corr + dy_neu_off;
            dz_neu = 1.0;

            fx_frag = px_frag / pz_frag + fx_frag_off;
            fy_frag = py_frag / pz_frag + fy_frag_off;
            fz_frag = 1.0;

            Double_t cos_ang =
                (dx_neu * fx_frag + dy_neu * fy_frag + 1.0) /
                (std::sqrt(dx_neu * dx_neu + dy_neu * dy_neu + 1.0) *
                 std::sqrt(fx_frag * fx_frag + fy_frag * fy_frag + 1.0));

            Double_t gamma_neu = 1.0 / std::sqrt(1.0 - beta_neu * beta_neu);
            Double_t gamma_frag = 1.0 / std::sqrt(1.0 - beta_frag * beta_frag);

            // Double_t m_f = M_frag * 1000.;

            Double_t m_f = M_frag;

            Erel = std::sqrt(m_f * m_f + m_neut * m_neut +
                             2.0 * gamma_neu * gamma_frag * m_f * m_neut *
                                 (1.0 - beta_neu * beta_frag * cos_ang)) -
                   m_f - m_neut;

            fErel->Fill(Erel * 1000);
            fDeltaBeta->Fill(beta_frag - beta_neu);

            beta_n_out = beta_neu;
        }

        ///////// MOMENTA DISTRIBUTIONS ///////////

        // P4 of the fragment
        fx_frag = px_frag / pz_frag;
        fy_frag = py_frag / pz_frag;
        fz_frag = 1.0;

        Double_t gamma_frag = 1.0 / std::sqrt(1.0 - beta_frag * beta_frag);
        Double_t m_f = M_frag;

        TVector3 u_frag(px_frag, py_frag, pz_frag);

        u_frag = u_frag.Unit();

        Double_t p_frag_lab = gamma_frag * m_f * beta_frag;
        TVector3 pvec_frag_lab = u_frag * p_frag_lab;
        TLorentzVector P4_frag_lab;

        P4_frag_lab.SetPxPyPzE(pvec_frag_lab.X(), pvec_frag_lab.Y(), pvec_frag_lab.Z(),
                               gamma_frag * m_f);

        // P4 of the incoming
        TVector3 u_proj(px_proj, py_proj, pz_proj);
        u_proj = u_proj.Unit();
        double beta_proj_mag = beta_proj;
        double gamma_proj = 1.0 / std::sqrt(1.0 - beta_proj_mag * beta_proj_mag);
        double p_proj_mag = gamma_proj * Mproj_GeV * beta_proj_mag;

        TVector3 p_proj_lab = u_proj * p_proj_mag;
        double E_proj_lab = gamma_proj * Mproj_GeV;

        TLorentzVector P4_proj_lab(
            p_proj_lab.X(),
            p_proj_lab.Y(),
            p_proj_lab.Z(),
            E_proj_lab);

        TLorentzVector P4_neu_lab;
        P4_neu_lab.SetPxPyPzE(0., 0., 0., 0.);

        if (fHasNeutrons)
        {
            Double_t nz_corr = pz_neu;
            dx_neu = px_neu / nz_corr;
            dy_neu = py_neu / nz_corr;

            Double_t gamma_neu = 1.0 / std::sqrt(1.0 - beta_neu * beta_neu);

            // TVector3 u_neu(dx_neu, dy_neu, 1.);
            TVector3 u_neu(px_neu, py_neu, pz_neu);
            u_neu = u_neu.Unit();

            Double_t p_neu_lab = gamma_neu * m_neut * beta_neu;

            TVector3 pvec_neu_lab = u_neu * p_neu_lab;

            P4_neu_lab.SetPxPyPzE(pvec_neu_lab.X(), pvec_neu_lab.Y(), pvec_neu_lab.Z(),
                                  gamma_neu * m_neut);
        }

        // Composed P4 of the fragment + neutron (LAB)
        TLorentzVector P4_sys_lab = P4_frag_lab;
        if (fHasNeutrons)
            P4_sys_lab = P4_neu_lab + P4_frag_lab;

        // Rotate everything into the RF of the projectile
        TVector3 u_proj_lab = P4_proj_lab.Vect().Unit();

        const double phi = u_proj_lab.Phi();
        const double theta = u_proj_lab.Theta();

        TRotation R;
        R.RotateZ(-phi);
        R.RotateY(-theta);

        auto RotateP4 = [&](const TLorentzVector &in) -> TLorentzVector
        {
            TLorentzVector out = in;
            TVector3 v = in.Vect();
            v.Transform(R);
            out.SetVect(v);
            return out;
        };

        TLorentzVector P4_sys_rot = RotateP4(P4_sys_lab);
        TLorentzVector P4_proj_rot = RotateP4(P4_proj_lab);

        TLorentzVector P4_sys_CM = P4_sys_rot;
        P4_sys_CM.Boost(TVector3(0.0, 0.0, -beta_proj)); // boost opposite to +Z

        px_sys = P4_sys_lab.Px();
        py_sys = P4_sys_lab.Py();
        pz_sys = P4_sys_lab.Pz();
        p_sys = P4_sys_lab.P();

        px_sys_CM = P4_sys_CM.Px();
        py_sys_CM = P4_sys_CM.Py();
        pz_sys_CM = P4_sys_CM.Pz();
        p_sys_CM = P4_sys_CM.P();

        TVector3 p3 = P4_sys_CM.Vect();

        px_rf_rot = p3.X();
        py_rf_rot = p3.Y();
        pz_rf_rot = p3.Z();

        pL = pz_rf_rot;
        pT = std::sqrt(px_rf_rot * px_rf_rot + py_rf_rot * py_rf_rot);

        // Other variables (angle in LAB kept as you had it)
        opa_lab = P4_neu_lab.Vect().Angle(P4_frag_lab.Vect());

        if (called)
            tout->Fill();
    }

    if (called)
    {
        fout->cd();
        tout->Write();
        fout->Close();
    }

    if (called)
        std::cout << "[KinTree] TTree saved in: " << outFile << "\n";

    f->Close();
}

void DataAnalysis::matchBeta()
{
    const double number = 100.;
    const double dBeta0 = -0.1;
    double dPaso = 0.2 / number;

    auto *gr = new TGraph(number);
    gr->SetLineColor(2);
    gr->SetLineWidth(4);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerStyle(7);

    for (int i = 0; i < number; i++)
    {
        const double dBeta = dBeta0 + i * dPaso;
        fBetaMatchValue = dBeta;
        fErel->Reset();
        fDeltaBeta->Reset();
        getData(false);
        int binmax = fDeltaBeta->GetMaximumBin();

        const int nSideBins = 3;
        int binL = std::max(1, binmax - nSideBins);
        int binR = std::min(fDeltaBeta->GetNbinsX(), binmax + nSideBins);

        double xL = fDeltaBeta->GetXaxis()->GetBinLowEdge(binL);
        double xR = fDeltaBeta->GetXaxis()->GetBinUpEdge(binR);

        double amp0 = fDeltaBeta->GetBinContent(binmax);
        double mean0 = fDeltaBeta->GetXaxis()->GetBinCenter(binmax);
        double sigma0 = fDeltaBeta->GetXaxis()->GetBinWidth(binmax) * 1.5;

        TF1 gaus("gaus_local", "gaus", xL, xR);
        gaus.SetParameters(amp0, mean0, sigma0);
        gaus.SetParLimits(2, 1e-6, (xR - xL));

        int fitStatus = fDeltaBeta->Fit(&gaus, "RQ0");

        double x_peak = mean0;
        if (fitStatus == 0)
        {
            x_peak = gaus.GetParameter(1);
        }

        gr->SetPoint(i, fBetaMatchValue, x_peak);
    }

    gr->Draw("ALP");
}