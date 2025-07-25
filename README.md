# SO-ASE
Code repository: Southern Ocean and Antarctic Sea ice Evolution (SO-ASE)

## Overview

This repository contains tools and workflows for studying the Southern Ocean and Antarctic sea ice evolution. It includes Python modules for data analysis and scripts for parallelized processing of model outputs.

### A) `so_ase` Python Module

The `so_ase` Python module provides utilities for analyzing and visualizing Southern Ocean and Antarctic sea ice data. Key features include:

- **Data Processing**: Functions for handling model outputs and observational datasets.
- **Visualization**: Tools for creating plots and maps of sea ice and ocean properties.

The python package can be installed using pip.

1. Clone the repository 
2. cd ./SO-ASE/
2. conda env create -f ./environment/environment.yml or conda create -f ./environment/environment.yml
3. conda activate so_ase
4. cd ./so_ase/
5. pip install -e .
6. python -m ipykernel install --user --name=so_ase --display-name "Python (so-ase)" (optional)

The -e will install the package in editable mode, so local changes are read at import.  
The package can now be imported from the so_ase environment in a python file and is further aailable as an ipykernel in Jupyter:     
import so_ase as so  

### B) `parallelized_workflows` Subdirectory

The `parallelized_workflows` subdirectory contains scripts for efficiently processing FESOM2 model outputs. These workflows are designed to run on high-performance computing systems using SLURM. Key features include:

- **Parallelized Interpolation**: Scripts for interpolating model data onto regular grids.

- **Parallelized Unrotating**: Scripts for unrotating model vector quantities.

### Usage

Refer to the respective documentation for detailed instructions on using the `so_ase` module and the parallelized workflows.
