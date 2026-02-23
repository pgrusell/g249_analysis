#!/bin/bash

# Load env variables
export repopath=$(realpath .)

# Create the data directory (if it does not exist yet)
mkdir -p data

# Create the results directory and its structure (if it does not exists yet)
mkdir -p results
mkdir -p results/{cal,sim,final}
mkdir -p sim/gen
