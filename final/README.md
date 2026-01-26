# Final Analysis

This section contains the macros and utils used to compute and obtain the different distributions for the final analysis. Here we explain the workflow to reproduce the analysis.

## Data retrieval

The first step is pass from the data of the unpacked `.root` files to something more usefull. In particular, the idea is to extract the data of the `.root` file for every specific ion/reaction and disregard events with other ions or without the complete information. 

For example, if we want to study only the product of the reaction $^{25}F$ (p,2p) $^{24}O$ decaying by one neutron emission (i.e. $^{25}F$ to $^{23}O$ + 1n), we will generate a filtered root file which will only contain data of this channel.

This is done using the macro `eventFilter.C` that has to be called specifying the setting, the reaction, and there is one more argument in the signature for testing.

## Data analysis

Once we have the specific filtered data, it's time to start with the data analysis itself. For this, the class dataAnalysis has been developed, and can be called from the runDataAnalysisMacro.

- ****

- **Beta matching**: In principle, the neutron and the fragment should have the same beta, however, due to energy losses and intrinsec uncertainties on the MDF, the momentum of the fragment is biased. To solve this problem, we need to match the beta of the neutron and the beta of the fragment through different methods. This is  

