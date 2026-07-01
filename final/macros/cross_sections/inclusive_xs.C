double efficiency(double E)
{
    if (E < 4.61)
        return 78.2688 / 100.;

    else
    {
        return (115.849 - E * 8.14434) / 100.;
    }
}

double calculateEfficiency(TH1F *hErel, double min, double max)
{

    double eff = 0;
    double N = 0;

    for (int i = 1; i <= hErel->GetNbinsX(); i++)
    {
        double energy = hErel->GetBinCenter(i);

        if (energy > min && energy < max)
        {

            eff += hErel->GetBinContent(i) * efficiency(energy);
            N += hErel->GetBinContent(i);
        }
    }

    return (double)eff / N;
}

void inclusive_xs()
{
    // Data files
    std::string unreactedFilePath = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/data_25F.root";
    std::string ox23FilePath = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/23O_analyzed_noNeuland.root";
    std::string ox23WithNeulandFilePath = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/23O_analyzed.root";

    // Unreacted cuts
    auto *unreacted_25f_cut = new TCutG("unreacted_25f_cut", 23);
    unreacted_25f_cut->SetPoint(0, 2.74329, 8.48631);
    unreacted_25f_cut->SetPoint(1, 2.75374, 8.55453);
    unreacted_25f_cut->SetPoint(2, 2.7594, 8.66688);
    unreacted_25f_cut->SetPoint(3, 2.76376, 8.80331);
    unreacted_25f_cut->SetPoint(4, 2.76637, 8.92368);
    unreacted_25f_cut->SetPoint(5, 2.7655, 9.04005);
    unreacted_25f_cut->SetPoint(6, 2.76289, 9.15642);
    unreacted_25f_cut->SetPoint(7, 2.75853, 9.26877);
    unreacted_25f_cut->SetPoint(8, 2.75418, 9.3731);
    unreacted_25f_cut->SetPoint(9, 2.74939, 9.4373);
    unreacted_25f_cut->SetPoint(10, 2.7446, 9.48545);
    unreacted_25f_cut->SetPoint(11, 2.73937, 9.48545);
    unreacted_25f_cut->SetPoint(12, 2.73153, 9.45736);
    unreacted_25f_cut->SetPoint(13, 2.72413, 9.36507);
    unreacted_25f_cut->SetPoint(14, 2.71934, 9.23667);
    unreacted_25f_cut->SetPoint(15, 2.71542, 9.08018);
    unreacted_25f_cut->SetPoint(16, 2.71454, 8.9598);
    unreacted_25f_cut->SetPoint(17, 2.71454, 8.83139);
    unreacted_25f_cut->SetPoint(18, 2.72064, 8.65885);
    unreacted_25f_cut->SetPoint(19, 2.72979, 8.53045);
    unreacted_25f_cut->SetPoint(20, 2.7385, 8.47427);
    unreacted_25f_cut->SetPoint(21, 2.74285, 8.47026);
    unreacted_25f_cut->SetPoint(22, 2.74329, 8.48631);
    unreacted_25f_cut->SetVarX("AoQ_frag");
    unreacted_25f_cut->SetVarY("Z_frag_est");

    // 23O cuts
    auto *reacted_23o_cutg = new TCutG("reacted_23o_cutg", 13);
    reacted_23o_cutg->SetVarX("AoQ_frag");
    reacted_23o_cutg->SetVarY("Z_frag_est");
    reacted_23o_cutg->SetTitle("Graph");
    reacted_23o_cutg->SetFillStyle(1000);
    reacted_23o_cutg->SetPoint(0, 2.78967, 8.38836);
    reacted_23o_cutg->SetPoint(1, 2.78589, 8.14281);
    reacted_23o_cutg->SetPoint(2, 2.794, 7.85634);
    reacted_23o_cutg->SetPoint(3, 2.81074, 7.68753);
    reacted_23o_cutg->SetPoint(4, 2.85233, 7.67218);
    reacted_23o_cutg->SetPoint(5, 2.87663, 7.83588);
    reacted_23o_cutg->SetPoint(6, 2.87933, 8.11212);
    reacted_23o_cutg->SetPoint(7, 2.87231, 8.30139);
    reacted_23o_cutg->SetPoint(8, 2.85503, 8.5009);
    reacted_23o_cutg->SetPoint(9, 2.83612, 8.56229);
    reacted_23o_cutg->SetPoint(10, 2.80966, 8.52648);
    reacted_23o_cutg->SetPoint(11, 2.79616, 8.44463);
    reacted_23o_cutg->SetPoint(12, 2.78967, 8.38836);
    reacted_23o_cutg->SetVarX("AoQ_frag");
    reacted_23o_cutg->SetVarY("Z_frag_est");

    double opaMin = 1.25;
    double opaMax = 1.65;
    TString condOpa = Form("califa_opa > %f && califa_opa < %f", opaMin, opaMax);

    ///////////// Unreacted fragment ///////////

    auto *f_25f = new TFile(unreactedFilePath.c_str(), "READ");
    auto *tr_25f = static_cast<TTree *>(f_25f->Get("FilterDataTree"));

    int kNumUnreacted = tr_25f->Draw("Z_frag_est:AoQ_frag", "unreacted_25f_cut", "goff");
    std::cout << "Unreacted  = " << kNumUnreacted << std::endl;

    ///////////// Outgoing fragment ///////////

    auto *f_23o = new TFile(ox23FilePath.c_str(), "READ");
    auto *tr_23o = static_cast<TTree *>(f_23o->Get("KinTree"));

    int kNumFragments = tr_23o->Draw("Z_frag_est:AoQ_frag", static_cast<TString>("reacted_23o_cutg &&") + condOpa, "goff");
    std::cout << "Fragments  = " << kNumFragments << std::endl;

    ///////////// Outgoing fragment (with neutrons) ///////////

    auto *f_23o_neu = new TFile(ox23WithNeulandFilePath.c_str(), "READ");
    auto *tr_23o_neu = static_cast<TTree *>(f_23o_neu->Get("KinTree"));

    auto *hErel = new TH1F("hErel", "Erel", 70, 0, 7.5);

    int kNumFragmentsWithNeuland = tr_23o_neu->Draw("Z_frag_est:AoQ_frag", static_cast<TString>("reacted_23o_cutg && Erel*1000 <7.5 &&") + condOpa, "goff");
    tr_23o_neu->Draw("Erel*1000>>hErel", static_cast<TString>("reacted_23o_cutg && Erel*1000 <7.5 &&") + condOpa, "");

    // Efficiencies
    double eff_p2p_tot = 0.566897;
    double eff_neu = calculateEfficiency(hErel, 0, 7.5);

    std::cout << "Fragments with one Neutron  = " << kNumFragmentsWithNeuland / eff_neu << std::endl;

    // Calculation of the cross section
    double reaction_rate = 0.024;
    double nproj = kNumUnreacted / (1. - reaction_rate);

    std::cout << "Projectiles  = " << nproj << std::endl;

    const double kMt = 1.000784;
    const double kRho = 0.0715;
    const double kZTarget = 5;
    const double kNa = 6.022E23;
    double xs = -TMath::Log(1. - kNumFragments / nproj / 1. / eff_p2p_tot) * kMt / kRho / kZTarget / kNa * 1e24 * 1e3;

    double R = kNumFragments / nproj / eff_p2p_tot;
    double C = kMt / kRho / kZTarget / kNa * 1e24 * 1e3;

    double dR = R * TMath::Sqrt(1.0 / kNumFragments + 1.0 / kNumUnreacted);
    double dSigma = C / (1.0 - R) * dR;

    std::cout << xs << " +/- " << dSigma << " [mb]\n";

    ////////////////

    xs = -TMath::Log(1. - kNumFragmentsWithNeuland / eff_neu / nproj / 1. / eff_p2p_tot) * kMt / kRho / kZTarget / kNa * 1e24 * 1e3;

    R = kNumFragmentsWithNeuland / eff_neu / nproj / eff_p2p_tot;
    C = kMt / kRho / kZTarget / kNa * 1e24 * 1e3;

    dR = R * TMath::Sqrt(1.0 / kNumFragmentsWithNeuland / eff_neu + 1.0 / kNumUnreacted);
    dSigma = C / (1.0 - R) * dR;

    std::cout << xs << " +/- " << dSigma << " [mb]\n";

    //////////////// Sigma vs Erel_max plot ////////////////

    // Histogram with extended range up to 15 MeV
    auto *hErel_wide = new TH1F("hErel_wide", "Erel wide", 150, 0, 15.);
    tr_23o_neu->Draw("Erel*1000>>hErel_wide",
                     static_cast<TString>("reacted_23o_cutg && Erel*1000 < 15. &&") + condOpa, "goff");

    // Scan Erel_max from 7.5 to 15 MeV
    const int nPoints = 100;
    const double erelMin = 7.5;
    const double erelMax = 15.;

    auto *gr = new TGraphErrors(nPoints);
    gr->SetTitle(";E_{rel}^{max} [MeV];#sigma [mb]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetLineWidth(2);

    for (int ip = 0; ip < nPoints; ip++)
    {
        double erelCut = erelMin + (erelMax - erelMin) * ip / (nPoints - 1.);

        // Integrate histogram up to erelCut for counts and weighted efficiency
        double nFrag = 0;
        double effSum = 0;

        for (int ib = 1; ib <= hErel_wide->GetNbinsX(); ib++)
        {
            double energy = hErel_wide->GetBinCenter(ib);
            if (energy > erelCut)
                break;

            double content = hErel_wide->GetBinContent(ib);
            nFrag += content;
            effSum += content * efficiency(energy);
        }

        if (nFrag < 1)
            continue;

        double effNeu = effSum / nFrag;
        double nFragCorr = nFrag / effNeu;

        double Ri = nFragCorr / nproj / eff_p2p_tot;
        double sigma_i = -TMath::Log(1. - Ri) * C;

        // Error propagation
        double dRi = Ri * TMath::Sqrt(1.0 / nFragCorr + 1.0 / kNumUnreacted);
        double dSigma_i = C / (1.0 - Ri) * dRi;

        gr->SetPoint(ip, erelCut, sigma_i);
        gr->SetPointError(ip, 0, dSigma_i);
    }

    auto *cSigma = new TCanvas("cSigma", "Sigma vs Erel max", 800, 600);
    gr->Draw("AP");
    cSigma->Update();
}