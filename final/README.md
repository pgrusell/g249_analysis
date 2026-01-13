# Final Analysis

This section contains the macros and utils used to compute and obtain the different distributions for the final analysis. Here we explain the workflow to reproduce the analysis.

## Data retrieval

The first step is pass from the data of the unpacked `.root` files to something more usefull. In particular, the idea is to extract the data of the `.root` file for every specific ion/reaction and disregard events with other ions or without the complete information. 

For example, if we want to study only the product of the reaction $^{25}F$ (p,2p) $^{24}O$ decaying by one neutron emission (i.e. $^{25}F$ to $^{23}O$ + 1n), we will generate a filtered root file which will only contain data of this channel.

This is done using the macro `eventFilter.C` that has to be called specifying the setting, the reaction, and there is one more argument in the signature for testing.

