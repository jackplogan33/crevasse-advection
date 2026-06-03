# Part 2: Stress Derivation and Crevasse Advection

This directory contains the analysis and figure-generation code for [CITATION].
It does not require running Part 1 &mdash; all necessary input data can be downloaded 
from the Zenodo archive below.

## Input Data

Processed velocity fields and supporting datasets are archived on Zenodo.

> **Zenodo record:** [DOI: TBD]

Download and place the data in `analysis/data/` before running the notebooks.


## Notebooks

| Notebook | Figures generated | Key outputs |
|----------|------------------|-------------|
| `01-stress_derivation.ipynb` | Figures 1, 3 | Strain rate tensor, deviatoric stress, von Mises stress, effective stress |
| `02-parcel_tracking.ipynb` | Figures 4&ndash;9 | Lagrangian parcel trajectories, crevasse advection timeseries |
| `03-supplemental.ipynb` | Supplemental figures and movies | S1&ndash;S7, Movies 1&ndash;5 |

Notebooks 01 and 02 can be run independently of each other. Notebook 03 depends 
on outputs from both and should be run last.

> [!NOTE]
> Movie generation in notebook 03 can be time-intensive; individual 
> cells are labeled so you can generate a subset if needed.

## Dependencies

All dependencies are covered by the shared `environment.yml` in the root of this 
repository. The core analysis functions live in `blue_ice_tools.py`, also in the 
root directory.

## Notes

- The `LagrangianTracking` class in `blue_ice_tools.py` is the central object for 
  Part 2. Its docstring and `02-parcel_tracking.ipynb` together are the best 
  entry point if you want to adapt the advection approach for your own work.