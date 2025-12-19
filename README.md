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
3. conda env create -f ./environment/environment.yml or conda create -f ./environment/environment.yml
4. conda activate so_ase
5. pip install -e .
6. python -m ipykernel install --user --name=so_ase --display-name "Python (so-ase)" (optional)

The -e will install the package in editable mode, so local changes are read at import.  

The package can now be imported from the so_ase environment in a python file.   
import so_ase as so  
It is further available as a kernel in Jupyter.     

### B) `workflows_*` Subdirectories

The `workflows+*` subdirectories contain scripts for various application, i.e. fesom input & mesh file generation and parallelized fesom output postprocessing.
