#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TEfficiency.h>

#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <algorithm>
#include <cmath>

const int NBINS = 641;
const double XLOW = -48.0;
const double XHI = 48.0;
const double BWIDTH = (XHI - XLOW) / NBINS;

const int NFOOTS = 4;
const char *branchNames[NFOOTS] = {"PosFoot5", "PosFoot6", "PosFoot7", "PosFoot8"};
const char *footLabels[NFOOTS] = {"Foot5", "Foot6", "Foot7", "Foot8"};

const std::vector<std::vector<int>> deadStripsPerFoot = {
    {257, 321, 385},                                                                           // Foot5
    {257, 321, 385},                                                                           // Foot6
    {65, 129, 193, 230, 257, 303, 304, 306, 307, 309, 321, 376, 385, 394, 449, 485, 513, 577}, // Foot7
    {193, 257, 321, 385, 449}                                                                  // Foot8
};

const int N_SKIP = 2;
const int NPAIRS = 6;

struct PairInfo
{
    int pairIdx;  // index into hCorr2D[]
    bool seedIsX; // true if the seed foot is the X axis of the histogram
};

bool isDeadStrip(int footNumber, int strip)
{
    int idx = footNumber - 5;
    if (idx < 0 || idx >= NFOOTS)
        return false;
    for (int ds : deadStripsPerFoot[idx])
    {
        if (ds == strip)
            return true;
    }
    return false;
}

std::set<std::pair<int, int>> getInterpolatedBins(int footIdxX, int footIdxY)
{
    std::set<int> deadSetX(deadStripsPerFoot[footIdxX].begin(),
                           deadStripsPerFoot[footIdxX].end());
    std::set<int> deadSetY(deadStripsPerFoot[footIdxY].begin(),
                           deadStripsPerFoot[footIdxY].end());

    std::set<std::pair<int, int>> interpBins;

    // Vertical bands (dead X): all bins (bx, by) with bx in deadSetX
    for (int bx : deadSetX)
        for (int by = 1; by <= NBINS; by++)
            interpBins.insert({bx, by});

    // Horizontal bands (dead Y): all bins (bx, by) with by in deadSetY
    for (int by : deadSetY)
        for (int bx = 1; bx <= NBINS; bx++)
            interpBins.insert({bx, by});

    return interpBins;
}

TH2D *interpolate2D(TH2D *hOrig, int footIdxX, int footIdxY, const char *name)
{
    TH2D *hCorr = (TH2D *)hOrig->Clone(name);

    const auto &deadX = deadStripsPerFoot[footIdxX];
    const auto &deadY = deadStripsPerFoot[footIdxY];
    std::set<int> deadSetX(deadX.begin(), deadX.end());
    std::set<int> deadSetY(deadY.begin(), deadY.end());

    // Interpolate vertical bands (dead X)
    // Group adjacent dead strips in X
    std::vector<int> sortedX(deadX.begin(), deadX.end());
    std::sort(sortedX.begin(), sortedX.end());

    std::vector<std::pair<int, int>> groupsX;
    {
        std::set<int> vis;
        for (int ds : sortedX)
        {
            if (vis.count(ds))
                continue;
            int first = ds, last = ds;
            vis.insert(ds);
            while (deadSetX.count(last + 1))
            {
                last++;
                vis.insert(last);
            }
            groupsX.push_back({first, last});
        }
    }

    for (auto &[gFirst, gLast] : groupsX)
    {
        for (int by = 1; by <= NBINS; by++)
        {
            // Find live neighbor to the left in X, skipping N_SKIP
            double valLeft = 0, errLeft = 0;
            int binLeft = -1;
            int skipped = 0;
            for (int bx = gFirst - 1; bx >= 1; bx--)
            {
                if (!deadSetX.count(bx))
                {
                    if (skipped < N_SKIP)
                    {
                        skipped++;
                        continue;
                    }
                    valLeft = hOrig->GetBinContent(bx, by);
                    errLeft = hOrig->GetBinError(bx, by);
                    binLeft = bx;
                    break;
                }
            }

            // Live neighbor to the right in X
            double valRight = 0, errRight = 0;
            int binRight = -1;
            skipped = 0;
            for (int bx = gLast + 1; bx <= NBINS; bx++)
            {
                if (!deadSetX.count(bx))
                {
                    if (skipped < N_SKIP)
                    {
                        skipped++;
                        continue;
                    }
                    valRight = hOrig->GetBinContent(bx, by);
                    errRight = hOrig->GetBinError(bx, by);
                    binRight = bx;
                    break;
                }
            }

            // Linear interpolation for each bin in the group
            for (int bx = gFirst; bx <= gLast; bx++)
            {
                double interp = 0;
                double errInterp = 0;
                if (binLeft >= 0 && binRight >= 0)
                {
                    double frac = (double)(bx - binLeft) / (binRight - binLeft);
                    interp = valLeft * (1.0 - frac) + valRight * frac;
                    errInterp = std::sqrt((1.0 - frac) * (1.0 - frac) * errLeft * errLeft +
                                          frac * frac * errRight * errRight);
                }
                else if (binLeft >= 0)
                {
                    interp = valLeft;
                    errInterp = errLeft;
                }
                else if (binRight >= 0)
                {
                    interp = valRight;
                    errInterp = errRight;
                }
                hCorr->SetBinContent(bx, by, interp);
                hCorr->SetBinError(bx, by, errInterp);
            }
        }
    }

    // Interpolate horizontal bands (dead Y)
    std::vector<int> sortedY(deadY.begin(), deadY.end());
    std::sort(sortedY.begin(), sortedY.end());

    std::vector<std::pair<int, int>> groupsY;
    {
        std::set<int> vis;
        for (int ds : sortedY)
        {
            if (vis.count(ds))
                continue;
            int first = ds, last = ds;
            vis.insert(ds);
            while (deadSetY.count(last + 1))
            {
                last++;
                vis.insert(last);
            }
            groupsY.push_back({first, last});
        }
    }

    // Temporary copy to avoid read/write contamination
    TH2D *hTemp = (TH2D *)hCorr->Clone("hTemp");

    for (auto &[gFirst, gLast] : groupsY)
    {
        for (int bx = 1; bx <= NBINS; bx++)
        {
            bool xAlsoDead = deadSetX.count(bx);

            // Live neighbor below in Y, skipping N_SKIP
            double valDown = 0, errDown = 0;
            int binDown = -1;
            int skipped = 0;
            for (int by = gFirst - 1; by >= 1; by--)
            {
                if (!deadSetY.count(by))
                {
                    if (skipped < N_SKIP)
                    {
                        skipped++;
                        continue;
                    }
                    valDown = hCorr->GetBinContent(bx, by);
                    errDown = hCorr->GetBinError(bx, by);
                    binDown = by;
                    break;
                }
            }

            // Live neighbor above in Y
            double valUp = 0, errUp = 0;
            int binUp = -1;
            skipped = 0;
            for (int by = gLast + 1; by <= NBINS; by++)
            {
                if (!deadSetY.count(by))
                {
                    if (skipped < N_SKIP)
                    {
                        skipped++;
                        continue;
                    }
                    valUp = hCorr->GetBinContent(bx, by);
                    errUp = hCorr->GetBinError(bx, by);
                    binUp = by;
                    break;
                }
            }

            for (int by = gFirst; by <= gLast; by++)
            {
                double interpY = 0;
                double errInterpY = 0;
                if (binDown >= 0 && binUp >= 0)
                {
                    double frac = (double)(by - binDown) / (binUp - binDown);
                    interpY = valDown * (1.0 - frac) + valUp * frac;
                    errInterpY = std::sqrt((1.0 - frac) * (1.0 - frac) * errDown * errDown +
                                           frac * frac * errUp * errUp);
                }
                else if (binDown >= 0)
                {
                    interpY = valDown;
                    errInterpY = errDown;
                }
                else if (binUp >= 0)
                {
                    interpY = valUp;
                    errInterpY = errUp;
                }

                if (xAlsoDead)
                {
                    // Crossing: average X interpolation (already in hCorr) with Y
                    double interpX = hCorr->GetBinContent(bx, by);
                    double errInterpX = hCorr->GetBinError(bx, by);
                    hTemp->SetBinContent(bx, by, 0.5 * (interpX + interpY));
                    hTemp->SetBinError(bx, by, 0.5 * std::sqrt(errInterpX * errInterpX + errInterpY * errInterpY));
                }
                else
                {
                    hTemp->SetBinContent(bx, by, interpY);
                    hTemp->SetBinError(bx, by, errInterpY);
                }
            }
        }
    }

    // Copy final result from hTemp to hCorr
    for (auto &[gFirst, gLast] : groupsY)
    {
        for (int bx = 1; bx <= NBINS; bx++)
        {
            for (int by = gFirst; by <= gLast; by++)
            {
                hCorr->SetBinContent(bx, by, hTemp->GetBinContent(bx, by));
                hCorr->SetBinError(bx, by, hTemp->GetBinError(bx, by));
            }
        }
    }

    delete hTemp;
    return hCorr;
}

double runStarMC(TH2D *h2D[NPAIRS],
                 int pairI[NPAIRS], int pairJ[NPAIRS],
                 PairInfo pairLookup[NFOOTS][NFOOTS],
                 TRandom3 &rng, int nMC,
                 int *nTotalOut = nullptr,
                 int *nLostOut = nullptr,
                 int *nLostByFoot = nullptr,
                 int *nLostByMult = nullptr,
                 int *nSeedUsed = nullptr)
{
    // Build marginals from the provided histograms
    TH1D *hMarg[NFOOTS];
    hMarg[0] = h2D[0]->ProjectionX("_marg0_mc");
    hMarg[1] = h2D[0]->ProjectionY("_marg1_mc");
    hMarg[2] = h2D[1]->ProjectionY("_marg2_mc");
    hMarg[3] = h2D[2]->ProjectionY("_marg3_mc");

    int nTotal = 0, nLost = 0;
    int lostByFoot[NFOOTS] = {0};
    int lostByMult[NFOOTS + 1] = {0};
    int seedUsed[NFOOTS] = {0};

    for (int mc = 0; mc < nMC; mc++)
    {
        // choose seed foot at random
        int seed = rng.Integer(NFOOTS);
        seedUsed[seed]++;

        // generate seed position from its 1D marginal
        double seedPos = hMarg[seed]->GetRandom();
        int seedBin = hMarg[seed]->FindBin(seedPos);

        // generate other 3 foots from conditionals
        double pos[NFOOTS];
        pos[seed] = seedPos;
        bool valid = true;

        for (int other = 0; other < NFOOTS; other++)
        {
            if (other == seed)
                continue;

            PairInfo &pi = pairLookup[seed][other];
            TH2D *h = h2D[pi.pairIdx];

            TH1D *hSlice = nullptr;
            if (pi.seedIsX)
                hSlice = h->ProjectionY(Form("_sl_%d_%d_%d", mc, seed, other),
                                        seedBin, seedBin);
            else
                hSlice = h->ProjectionX(Form("_sl_%d_%d_%d", mc, seed, other),
                                        seedBin, seedBin);

            if (hSlice->Integral() <= 0)
            {
                valid = false;
                delete hSlice;
                break;
            }
            pos[other] = hSlice->GetRandom();
            delete hSlice;
        }

        if (!valid)
            continue;

        // convert positions to strips and check dead strips
        nTotal++;
        int nDeadHit = 0;
        for (int i = 0; i < NFOOTS; i++)
        {
            int strip = (int)((pos[i] - XLOW) / BWIDTH) + 1;
            if (isDeadStrip(i + 5, strip))
            {
                lostByFoot[i]++;
                nDeadHit++;
            }
        }
        if (nDeadHit > 0)
            nLost++;
        lostByMult[nDeadHit]++;
    }

    // Copy out optional counters
    if (nTotalOut)
        *nTotalOut = nTotal;
    if (nLostOut)
        *nLostOut = nLost;
    if (nLostByFoot)
        for (int i = 0; i < NFOOTS; i++)
            nLostByFoot[i] = lostByFoot[i];
    if (nLostByMult)
        for (int m = 0; m <= NFOOTS; m++)
            nLostByMult[m] = lostByMult[m];
    if (nSeedUsed)
        for (int i = 0; i < NFOOTS; i++)
            nSeedUsed[i] = seedUsed[i];

    // Clean up marginals
    for (int i = 0; i < NFOOTS; i++)
        delete hMarg[i];

    double eff = (nTotal > 0) ? (double)(nTotal - nLost) / nTotal : 0.0;
    return eff;
}

void footStripEff()
{
    gStyle->SetOptStat(0);

    // ─── Open file ────────────────────────────────────────────────
    TString path = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_25F_test.root";
    auto *f = new TFile(path, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Error opening " << path << std::endl;
        return;
    }

    auto *tr = static_cast<TTree *>(f->Get("FilterDataTree"));
    if (!tr)
    {
        std::cerr << "Error: FilterDataTree not found" << std::endl;
        return;
    }

    double posVar[NFOOTS];
    for (int i = 0; i < NFOOTS; i++)
        tr->SetBranchAddress(branchNames[i], &posVar[i]);

    Long64_t nEntries = tr->GetEntries();
    std::cout << "Events in tree: " << nEntries << std::endl;

    // Create and fill the 6 original 2D histograms
    int pairI[NPAIRS], pairJ[NPAIRS];
    {
        int p = 0;
        for (int i = 0; i < NFOOTS; i++)
            for (int j = i + 1; j < NFOOTS; j++)
            {
                pairI[p] = i;
                pairJ[p] = j;
                p++;
            }
    }

    TH2D *hOrig2D[NPAIRS];
    TH2D *hCorr2D[NPAIRS];

    for (int p = 0; p < NPAIRS; p++)
    {
        int fi = pairI[p], fj = pairJ[p];
        hOrig2D[p] = new TH2D(
            Form("hOrig2D_%s_%s", footLabels[fi], footLabels[fj]),
            Form("%s vs %s;%s [mm];%s [mm]",
                 footLabels[fi], footLabels[fj],
                 footLabels[fi], footLabels[fj]),
            NBINS, XLOW, XHI, NBINS, XLOW, XHI);
        // Enable Sumw2 so ROOT stores sqrt(N) errors from the start
        hOrig2D[p]->Sumw2();
    }

    std::cout << "Filling 2D histograms..." << std::endl;
    for (Long64_t ev = 0; ev < nEntries; ev++)
    {
        tr->GetEntry(ev);
        for (int p = 0; p < NPAIRS; p++)
            hOrig2D[p]->Fill(posVar[pairI[p]], posVar[pairJ[p]]);
    }

    // ─── 2. Interpolate dead bands in each 2D histogram ──────────
    std::cout << "\nInterpolating the 6 2D histograms..." << std::endl;
    for (int p = 0; p < NPAIRS; p++)
    {
        int fi = pairI[p], fj = pairJ[p];
        std::cout << "  " << footLabels[fi] << " vs " << footLabels[fj]
                  << " (" << deadStripsPerFoot[fi].size() << " + "
                  << deadStripsPerFoot[fj].size() << " dead strips)..."
                  << std::flush;
        hCorr2D[p] = interpolate2D(hOrig2D[p], fi, fj,
                                   Form("hCorr2D_%s_%s", footLabels[fi], footLabels[fj]));
        std::cout << " done." << std::endl;
    }

    std::vector<std::set<std::pair<int, int>>> interpBinSets(NPAIRS);
    for (int p = 0; p < NPAIRS; p++)
        interpBinSets[p] = getInterpolatedBins(pairI[p], pairJ[p]);

    auto *c1 = new TCanvas("c1", "Original", 1800, 1200);
    c1->Divide(3, 2);
    for (int p = 0; p < NPAIRS; p++)
    {
        c1->cd(p + 1);
        gPad->SetLogz();
        hOrig2D[p]->Draw("COLZ");
    }

    auto *c2 = new TCanvas("c2", "Corrected", 1800, 1200);
    c2->Divide(3, 2);
    for (int p = 0; p < NPAIRS; p++)
    {
        c2->cd(p + 1);
        gPad->SetLogz();
        hCorr2D[p]->Draw("COLZ");
    }

    std::cout << "\n6 corrected 2D histograms available on the heap:" << std::endl;
    for (int p = 0; p < NPAIRS; p++)
        std::cout << "  hCorr2D[" << p << "] = "
                  << footLabels[pairI[p]] << " vs " << footLabels[pairJ[p]]
                  << "  (" << hCorr2D[p]->GetName() << ")" << std::endl;

    PairInfo pairLookup[NFOOTS][NFOOTS];
    for (int p = 0; p < NPAIRS; p++)
    {
        int fi = pairI[p], fj = pairJ[p];
        pairLookup[fi][fj] = {p, true};
        pairLookup[fj][fi] = {p, false};
    }

    TRandom3 rng(12345);
    const int nMC = 100000;

    int nTotal = 0, nLost = 0;
    int nLostByFoot[NFOOTS] = {0};
    int nLostByMultiplicity[NFOOTS + 1] = {0};
    int nSeedUsed[NFOOTS] = {0};

    std::cout << "\nRunning nominal MC (" << nMC << " attempts)..." << std::endl;

    double effNominal = runStarMC(hCorr2D, pairI, pairJ, pairLookup, rng, nMC,
                                  &nTotal, &nLost, nLostByFoot,
                                  nLostByMultiplicity, nSeedUsed);

    // Repeat MC to obtain uncertainties
    const int N_TOYS = 100;
    std::vector<double> toyEfficiencies(N_TOYS);

    std::cout << "\nRunning " << N_TOYS << " toy MCs to estimate interpolation systematic..."
              << std::endl;

    for (int itoy = 0; itoy < N_TOYS; itoy++)
    {
        if ((itoy + 1) % 10 == 0)
            std::cout << "  Toy " << (itoy + 1) << "/" << N_TOYS << std::endl;

        // Clone the corrected histograms
        TH2D *hToy[NPAIRS];
        for (int p = 0; p < NPAIRS; p++)
        {
            hToy[p] = (TH2D *)hCorr2D[p]->Clone(Form("hToy_%d_%d", itoy, p));

            // Fluctuate only the interpolated bins
            for (auto &[bx, by] : interpBinSets[p])
            {
                double val = hCorr2D[p]->GetBinContent(bx, by);
                double err = hCorr2D[p]->GetBinError(bx, by);
                if (err <= 0)
                    continue; // nothing to fluctuate

                double fluct = rng.Gaus(val, err);
                if (fluct < 0)
                    fluct = 0; // bin content cannot be negative
                hToy[p]->SetBinContent(bx, by, fluct);
            }
        }

        // Run the star MC on the fluctuated histograms
        toyEfficiencies[itoy] = runStarMC(hToy, pairI, pairJ, pairLookup,
                                          rng, nMC);

        // Clean up toy histograms
        for (int p = 0; p < NPAIRS; p++)
            delete hToy[p];
    }

    // Compute mean and RMS of the toy efficiency distribution
    double toyMean = 0, toyRMS = 0;
    for (double e : toyEfficiencies)
        toyMean += e;
    toyMean /= N_TOYS;
    for (double e : toyEfficiencies)
        toyRMS += (e - toyMean) * (e - toyMean);
    toyRMS = std::sqrt(toyRMS / (N_TOYS - 1));

    // Extract the 68% central interval (analogous to 1-sigma)
    std::sort(toyEfficiencies.begin(), toyEfficiencies.end());
    int idxLo = (int)(0.5 * (1.0 - 0.683) * N_TOYS);
    int idxHi = (int)(0.5 * (1.0 + 0.683) * N_TOYS);
    if (idxHi >= N_TOYS)
        idxHi = N_TOYS - 1;
    double toy68Lo = toyEfficiencies[idxLo];
    double toy68Hi = toyEfficiencies[idxHi];

    // Results
    const double confLevel = 0.683;
    int nPassed = nTotal - nLost;
    double effTotalLo = TEfficiency::ClopperPearson(nTotal, nPassed, confLevel, false);
    double effTotalHi = TEfficiency::ClopperPearson(nTotal, nPassed, confLevel, true);

    std::cout << Form("  Valid events:         %d", nTotal) << std::endl;
    std::cout << Form("  Lost events:          %d", nLost) << std::endl;
    std::cout << Form("  Total efficiency:     %.4f %% [- %.4f, + %.4f] %% (stat, Clopper-Pearson)",
                      100.0 * effNominal,
                      100.0 * (effNominal - effTotalLo),
                      100.0 * (effTotalHi - effNominal))
              << std::endl;
    std::cout << "  Per foot:" << std::endl;
    for (int i = 0; i < NFOOTS; i++)
    {
        int nPassedFoot = nTotal - nLostByFoot[i];
        double effFoot = (double)nPassedFoot / nTotal;
        double effFootLo = TEfficiency::ClopperPearson(nTotal, nPassedFoot, confLevel, false);
        double effFootHi = TEfficiency::ClopperPearson(nTotal, nPassedFoot, confLevel, true);
        std::cout << Form("    %s: %2d dead strips -> lost %7d (%.4f %%)  "
                          "Eff: %.4f %% [- %.4f, + %.4f]",
                          footLabels[i],
                          (int)deadStripsPerFoot[i].size(),
                          nLostByFoot[i],
                          100.0 * nLostByFoot[i] / nTotal,
                          100.0 * effFoot,
                          100.0 * (effFoot - effFootLo),
                          100.0 * (effFootHi - effFoot))
                  << std::endl;
    }

    std::cout << "  Dead-foot multiplicity per event:" << std::endl;
    for (int m = 0; m <= NFOOTS; m++)
    {
        std::cout << Form("    %d foots hit: %8d events (%.4f %%)",
                          m, nLostByMultiplicity[m],
                          100.0 * nLostByMultiplicity[m] / nTotal)
                  << std::endl;
    }
    std::cout << "  Seed usage:" << std::endl;
    for (int i = 0; i < NFOOTS; i++)
    {
        std::cout << Form("    %s: %d times (%.1f %%)",
                          footLabels[i], nSeedUsed[i],
                          100.0 * nSeedUsed[i] / nMC)
                  << std::endl;
    }

    std::cout << Form("  Nominal efficiency:   %.4f %%", 100.0 * effNominal) << std::endl;
    std::cout << Form("  Toy mean efficiency:  %.4f %%", 100.0 * toyMean) << std::endl;
    std::cout << Form("  Toy RMS (syst):       %.4f %%", 100.0 * toyRMS) << std::endl;
    std::cout << Form("  Toy 68%% interval:     [%.4f, %.4f] %%",
                      100.0 * toy68Lo, 100.0 * toy68Hi)
              << std::endl;
    std::cout << Form("  Stat (C-P):           [- %.4f, + %.4f] %%",
                      100.0 * (effNominal - effTotalLo),
                      100.0 * (effTotalHi - effNominal))
              << std::endl;
    std::cout << Form("  Combined (quad sum):  +/- %.4f %%",
                      100.0 * std::sqrt(toyRMS * toyRMS +
                                        std::pow(0.5 * (effTotalHi - effTotalLo), 2)))
              << std::endl;
}