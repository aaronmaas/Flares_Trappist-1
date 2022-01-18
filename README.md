# Flares on Trappist-1 with K2 and MuSCAT 1 and 2 data

This repository contains the neccesarry modules and methods to excess all results from following work in Progress:

1. Maas, A. J. et al (2022, in prep.). *Lower temperatures than expected for TRAPPIST-1 flares*
2. Maas, A. J. et al (2022, in prep.). *Are late M-dwarfs flares really that cold? *


#### How to navigate this repository

**Nomenclature:**

- File names that contain `TRAPPIST_1` or `AUMic` belong to project 1. 
- File names that contain `Other_Star` belong to project 2.
- File names with neither belong to both, like `notebooks/exampple.py`.
- File names with a `_` prefix are work in progress.

**Folders:**

- `notebooks/` contains all the notebooks and scripts 
  - also contains the modules used in both, stored in `funcs/`
- `data/` contains ancillary data, such as the K2, M1, M2 transmission curves but also theoretical calculations
- `results/` contain the tables and figures for the papers

## Project 1: Lower than observed temperatures for TRAPPIST-1 flares

Notebooks that produce the figures and tables in the paper found in `notebooks/`:

- 

## Project 1: Are late M-dwarf flares really that cold? 

Notebooks that produce the figures and tables in the paper found in `notebooks/`:

- 


## Required Python packages

- numpy
- scipy
- astropy
- matplotlib
- lightkurve
- altaipony
- PyTransit
- Muscat2transitpipeline
- emcee
- corner`
