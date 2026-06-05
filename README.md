# Crevasse Advection at Shirase Glacier

This repository accompanies the publication [CITATION]. It contains two independent
components: a Sentinel-1 SAR velocity processing pipeline (Part 1) and the scientific 
analysis and figure-generation code for the paper (Part 2).

## Where do you want to start?

**I want to reproduce the paper figures**
* The processed velocity data is archived on Zenodo [DOI: 10.5281/zenodo.20535986].
  Download it and head to [`02-analysis/`](./02-analysis/README.md).

**I want to run the velocity pipeline myself**
* See [`01-processing/`](./01-processing/README.md) for installation and usage.

**I want to adapt the pipeline for a different region or sensor**
* Start with the [`01-processing/`](./01-processing/README.md) README, then focus on `02-run_topsapp.ipynb`
  and `03-postprocess.ipynb` where the key parameters are documented.

---

## Repository Structure

```
crevasse-advection/
├── 01-processing/         # Part 1: Sentinel-1 pixel offset tracking with ISCE2
│   ├── offsets/           # Python scripts for automated offset processing
│   ├── 01-setup.ipynb        
│   ├── 02-run_topsapp.ipynb  
│   └── 03-postprocess.ipynb  
├── 02-analysis/           # Part 2: Stress derivation, parcel tracking, figures
│   ├── blue_ice_tools.py  # Core utility library
│   ├── 01-stress_derivation.ipynb
│   ├── 02-parcel_tracking.ipynb
│   └── 03-supplemental.ipynb
├── environment.yml
├── data/                  # Placeholder -- filled by downloads from Zenodo or Part 1
├── media/                 # Placeholder -- filled by Part 2 notebooks
├── environment.yml
├── setup.sh
└── README.md              # You are here
```

---

## Quick Start

### Prerequisites

* [Conda](https://docs.conda.io/en/latest/) (Miniconda or Anaconda)
* An [Earthdata](https://urs.earthdata.nasa.gov/) account (Part 1 only)

### Installation

A pre-built JupyterHub image with all dependencies configured is avaible in quay.io:

```
quay.io/jackplogan/isce-image:0e337edc16fc
```

For a local install, create the conda environment and run the setup script:

```bash
git clone https://github.com/jackplogan33/crevasse-advection/
cd crevasse-advection
conda env create -f environment.yml
conda activate isce2
./setup.sh
conda deactivate && conda activate isce2   # reactivate to apply PATH changes
```

See [`01-processing/README.md`](./01-processing/README.md) for full installation details and troubleshooting.

---

## Overview

### Part 1: Velocity Processing

Surface velocity fields are derived from consecutive Sentinel-1 SAR image pairs using amplitude cross correlations (pixel offset tracking) via the [ISCE2](https://github.com/isce-framework/isce2) `topsApp` pipeline

The workflow covers:
1. Automated SAFE and orbit file download via the ASF Python API
2. Per-pair topsApp configuration and execution
3. Dense offset stacking along a `mid_date` time dimension
4. Outlier masking and biharmonic inpainting
5. Conversion from pixel offsets to surface velocity in m/yr

All parameters are controlled through a single `config.json` file. 
An in-depth explanation of the pipeline and configuration of parameters is provided in [`01-processing/README.md`](./01-processing/README.md).

### Part 2: Stress Derivation and Crevasse Advection

Starting from the velocity time series, Part 2 derives the strain rate tensor and associated stress fields, then tracks ice parcels through the velocity field to study crevasse advection.

Key steps:
1. Strain rate computation via Savitzky-Golay filtering or finite differences
2. Rotation of the strain tensor into along- and across-flow components
3. Deviatoric and Cauchy stress calculation following Glen's flow law
4. Lagrangian parcel tracking with 4-th order Runge-Kutta integration
5. Ensemble advection over randomly sampled starting locations

The central objects are `calc_strain_stress()` and the `LagrangianTracking` class showcased in `02-analysis/blue_ice_tools.py`.

---

## Citation

If you use this code or data, please cite both the paper and this repository:

> **Paper:** [TBD]    
> **Data:**  [DOI: 10.5281/zenodo.20535986] (reserved, not published)    
> **Code:** [DOI: 10.5281/zenodo.20550146] (reserved, not published)    
