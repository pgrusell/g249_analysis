#include "dataAnalysis.cpp"

void runDataAnalysis()
{
    auto *da = new DataAnalysis("/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_23O_test.root", "23O_analyzed_test.root", true);
    da->setOffsetsFromTxt("23O1n.txt");
    // da->matchBeta();
    da->getData();
}