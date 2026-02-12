# g249_analysis

These are the codes to analyze the g249 experiment data. It will take as input the hit/cluster (or cal, if needed) data from the corresponding rootfile provided by [R3BRoot](https://github.com/R3BRootGroup/R3BRoot) and use it to generate the corresponding histograms and other analysis. 

## Usage

The codes are based on ROOT macros and use R3BRoot for data packaging so only a recent installation of R3BRoot is needed. 

After setting up the corresponding variables of R3BRoot, we have to generate the data directory structure and launch the path variables by executing calling the `config.sh` script.

## Structure

- doc: Information about the code, analysis, physics... (work in progress).
- cal: Macros and classes to perform calibrations of FOOTs and NeuLAND.
- sim: Classes to perform simulations of the setup (work in progress).
- theory: Theoretical calculations to compare with experimental data (wrong for the moment).
- utils: Helpers to plot stuff.
- final: Macros to obtain final physical results.

