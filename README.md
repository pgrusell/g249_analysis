# g249_analysis

These are the codes to analyze the g249 experiment data. It will take as input the hit/cluster (or cal, if needed) data from the corresponding rootfile provided by [R3BRoot](https://github.com/R3BRootGroup/R3BRoot) and use it to generate the corresponding histograms and other analysis. 

## Usage

The codes are based on ROOT macros and use R3BRoot for data packaging so only a recent installation of R3BRoot is needed. 

After setting up the corresponding variables of R3BRoot, we have to generate the data directory structure and launch the path variables by executing calling the `config.sh` script.

Once this is done all the macros can be executed using the command:

```bash
root -l macro.C
```

(with the different arguments depending on the macro itself).

## Structure

- doc: Information about the code, analysis, physics and usage.
- cal: Macros to test the calibrations of different detectors (currently only for eta correction on FOOTs).

_Under development:_

- sim: Simulation of the experiment.
- track: Reconstruction algorithms, trained with simulated data, to retrieve vertex, B$\rho$, A/Q and other physical variables for particle ID.
- final: Macros to obtain the physical observables from the data.