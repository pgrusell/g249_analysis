void plot(const char *fname = "momdis.txt")
{
    gStyle->SetOptStat(0);

    std::ifstream in(fname);
    if (!in.is_open())
    {
        std::cerr << "ERROR: cannot open file: " << fname << std::endl;
        return;
    }

    // Full vectors (no filter): Qz and Qy
    std::vector<double> Qi_all, dQz_all, dQy_all;

    // Filtered vectors (only y>0): Qt and Q
    std::vector<double> Qi_t, dQt_pos;
    std::vector<double> Qi_q, dQ_pos;

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        std::stringstream ss(line);

        double qi, z, t, y, q;
        if (!(ss >> qi >> z >> t >> y >> q))
            continue;

        // Keep all for Qz and Qy (as in your table)
        Qi_all.push_back(qi);
        dQz_all.push_back(z);
        dQy_all.push_back(y);

        // Only positive for Qt
        if (qi > 0.0)
        {
            Qi_t.push_back(qi);
            dQt_pos.push_back(t);
        }

        // Only positive for Q
        if (qi > 0.0)
        {
            Qi_q.push_back(qi);
            dQ_pos.push_back(q);
        }
    }

    if (Qi_all.empty())
    {
        std::cerr << "ERROR: no data parsed from file. Check formatting." << std::endl;
        return;
    }

    auto makeGraph = [&](const char *name, const char *title,
                         const std::vector<double> &x, const std::vector<double> &y)
    {
        if (x.size() != y.size() || x.empty())
        {
            std::cerr << "ERROR: graph " << name << " has invalid sizes: "
                      << x.size() << " vs " << y.size() << std::endl;
            return (TGraph *)nullptr;
        }
        auto *g = new TGraph((int)x.size());
        g->SetName(name);
        g->SetTitle(title);
        for (int i = 0; i < (int)x.size(); ++i)
            g->SetPoint(i, x[i], y[i]);
        g->SetLineWidth(2);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.6);
        return g;
    };

    TGraph *gQz = makeGraph("gQz", "dS/dQ_{z} vs Q_{i}", Qi_all, dQz_all);
    TGraph *gQy = makeGraph("gQy", "dS/dQ_{y} vs Q_{i}", Qi_all, dQy_all);
    TGraph *gQt = makeGraph("gQt", "dS/dQ_{t} vs Q_{i} (t>0)", Qi_t, dQt_pos);
    TGraph *gQ = makeGraph("gQ", "dS/dQ vs Q_{i} (Q>0)", Qi_q, dQ_pos);

    if (!gQz || !gQy || !gQt || !gQ)
        return;

    auto *c = new TCanvas("c_mom", "Momentum distributions", 1200, 900);
    c->Divide(2, 2);

    auto setupPad = [&](int ipad, TGraph *g, const char *ytitle)
    {
        c->cd(ipad);
        gPad->SetGrid(1, 1);
        g->Draw("ALP");
        g->GetXaxis()->SetTitle("Q_{i} (MeV/c)");
        g->GetYaxis()->SetTitle(ytitle);
        g->GetXaxis()->CenterTitle(true);
        g->GetYaxis()->CenterTitle(true);
    };

    setupPad(1, gQz, "dS/dQ_{z} (mb / (MeV/c))");
    setupPad(2, gQt, "dS/dQ_{t} (mb / (MeV/c))");
    setupPad(3, gQy, "dS/dQ_{y} (mb / (MeV/c))");
    setupPad(4, gQ, "dS/dQ (mb / (MeV/c))");

    c->Update();
    c->SaveAs("momentum_distributions_vs_Qi_filtered.pdf");
    c->SaveAs("momentum_distributions_vs_Qi_filtered.png");

    std::cout << "Saved: momentum_distributions_vs_Qi_filtered.pdf and .png" << std::endl;
}