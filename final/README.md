# Final Analysis

This section contains the macros and utils used to compute and obtain the different distributions for the final analysis. Here we explain the workflow to reproduce the analysis.

## Data reduction

The first step is pass from the data of the unpacked `.root` files to something more usefull. In particular, the idea is to extract the data of the `.root` file for every specific ion/reaction and disregard events with other ions or without the complete information. 

For example, if we want to study only the product of the reaction $^{25}F$ (p,2p) $^{24}O$ decaying by one neutron emission (i.e. $^{25}F$ to $^{23}O$ + 1n), we will generate a filtered root file which will only contain data of this channel.

This is done using the macro `eventFilter.C` that has to be called specifying the setting, the reaction, and there is one more argument in the signature for testing. An example of use is the following:

```bash
root -l 'eventFilter.C("25F.txt", "25F23O", false)'
```

There, the file `settings/25F.txt` contains the absolute paths to all the unpacked lmd files for the $^{25}F$ setting. 


## Data analysis

Once we have the specific filtered data, it's time to start with the data analysis itself. For this, the class dataAnalysis has been developed, and can be called from the runDataAnalysisMacro. The steps to use this macro are the following:

- **Settings:** First of all, in the folder ´settings´ we need to define a setting for the reactions we are interested in. For example `23O1n.txt` or `24O.txt`. This file has to contain offsets for the variables dx, dy, fx and fy. To get them, first create an empty file for the reaction channel and run the macro (modifiying the corresponding paths) runDataAnalysis. This will generate a rootfile that stores dx, dy... and you can put the offsets in every line of the .txt file (order fx\nfy\ndx\ndy).

- **Beta matching**: In principle, the neutron and the fragment should have the same beta, however, due to energy losses and intrinsic uncertainties on the MDF, the momentum of the fragment is biased. To solve this problem, we need to match the beta of the neutron and the beta of the fragment through different methods. This offset should be added in the last line of the setting file and can be calculated running the macro runDataAnalysis in betaMatch mode. This will plot a line that should pass to (0,0). The x-coordinate for y = 0 is the beta offset.

Once the settings file is ready, we can run runDataAnalysis in normal mode to get a file with momentum distributions and Erel calculated. Finally all these results can be plotted with the macro plot_23O_24O.C. 

### Other macros WIP

- **fiberAcceptance.C**: Macro to study the position of the beamspot in the fibers (currently not compiling).

- **fitMomentaDist.C**: Macro to fit momentum distributions to theoretical calculations (since we don't have these calculations this macro has not been tested yet, so it's very preliminary).



