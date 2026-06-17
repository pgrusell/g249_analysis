#include "dataAnalysis.cpp"

void runDataAnalysis(TString mode = "ana23")
{

    if (mode == "ana23")
    {
        auto *da = new DataAnalysis("/nucl_lustre/pablogrusell/g249/g249_analysis/results/data_23O.root", "23O_analyzed_test.root", true);
        da->setOffsetsFromTxt("23O1n.txt");
        da->getData();
    }

    if (mode == "ana24")
    {
        auto *da = new DataAnalysis("/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_24O_test.root", "24O_analyzed_test.root", false);
        // da->setOffsetsFromTxt("23O1n.txt");
        da->getData();
    }

    if (mode == "ana22O")
    {
        auto *da = new DataAnalysis("/home/e12exp/ssd/PhD_analysis/g249_analysis/results/final/data_22O.root", "22O_analyzed.root", true);
        da->setOffsetsFromTxt("22O.txt");
        da->getData();
    }

    if (mode == "betaMatching")
    {
        auto *da = new DataAnalysis("/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_23O_test.root", "23O_analyzed_test.root", true);
        da->setOffsetsFromTxt("23O1n.txt");
        da->matchBeta();
    }

        if (mode == "betaMatching22O")
    {
        auto *da = new DataAnalysis("/home/e12exp/ssd/PhD_analysis/g249_analysis/results/final/data_22O.root", "22O_analyzed.root", true);
        da->setOffsetsFromTxt("22O.txt");
        da->matchBeta();
    }

}
