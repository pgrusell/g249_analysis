void fiberAcceptance()
{
    std::vector<TString> files = {
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260129_160553_160_5.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260129_160624_154_5.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260129_220002_160_1.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260129_220043_152_1.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_050528_152_2.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260129_230235_160_6.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260129_222241_154_1.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_003935_154_6.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_043553_160_2.root",

        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_055357_160_7.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_064307_154_2.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_075827_154_7.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_112817_160_3.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_123817_160_8.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_151342_154_8.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_130412_154_3.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_183646_160_4.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_192006_160_9.root",

        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_202101_154_4.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260130_201448_154_9.root",
        "/nucl_lustre/pablogrusell/g249/root_files/unpack_trans/g249_all_det_offline_0001_20260131_031444_154_10.root",
    };

    // Create a TChain to read the "evt" tree from all files
    TChain *chain = new TChain("evt");

    for (const auto &file : files)
    {
        chain->Add(file);
    }
    o

            std::cout
        << "Total entries: " << chain->GetEntries() << std::endl;

    // Event loop

    double x_fib30;
    double x_fib31;
    double x_fib33;
    double y_fib32;

    double in_aoq;
    double in_Z;

    double out_aoq;
    double out_Z;

    chain->SetBranchAddress("Fi32Hit.fY", &y_fib32);
    chain->SetBranchAddress("Fi31Hit.fX", &x_fib31);
    chain->SetBranchAddress("Fi30Hit.fX", &x_fib30);
    chain->SetBranchAddress("Fi33Hit.fX", &x_fib33);

    chain->SetBranchAddress("FrsData.fAq", &in_aoq);
    chain->SetBranchAddress("FrsData.fZ", &in_Z);
    chain->SetBranchAddress("FragmentMDFTrack.fCharge", &out_Z);
    chain->SetBranchAddress("FragmentMDFTrack.fMass", &out_aoq);

    auto *x30vsY = new TH2F("x30vsY", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);
    auto *x31vsY = new TH2F("x31vsY", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);
    auto *x33vsY = new TH2F("x33vsY", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);

    auto *x30vsY_24 = new TH2F("x30vsY_24", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);
    auto *x31vsY_24 = new TH2F("x31vsY_24", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);
    auto *x33vsY_24 = new TH2F("x33vsY_24", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);

    auto *x30vsY_23 = new TH2F("x30vsY_23", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);
    auto *x31vsY_23 = new TH2F("x31vsY_23", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);
    auto *x33vsY_23 = new TH2F("x33vsY_23", ";x [cm]; y[cm]", 50, -20, 20, 50, -20, 20);

    Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        std::cout << i << std::endl;
        chain->GetEntry(i);
    }

    std::cout << "Done!" << std::endl;
}