# Crevasse Advection at Shirase Glacier

This repository accompanies the publication [CITATION]. It contains two independent
components: a Sentinel-1 SAR velocity processing pipeline (Part 1) and the scientific 
analysis and figure-generation code for the paper (Part 2).

## Where do you want to start?

**I want to reproduce the paper figures**
* The processed velocity data is archived on Zenodo [DOI: TBD].
  Download it and head to [`analysis/`](./analysis/README.md).

**I want to run the velocity pipeline myself**
* See [`processing/`](./processing/README.md) for installation and usage.

**I want to adapt the pipeline for a different region or sensor**
* Start with the [`processing/`](./processing/README.md) README, then focus on `02-run_topsapp.ipynb`
  and `03-postprocess.ipynb` where the key parameters are documented.

## Repository Structure

```
crevasse-advection/
├── 01-processing/         # Part 1: Sentinel-1 pixel offset tracking with ISCE2
│   ├── offsets/           # Python scripts for automated offset processing
│   └── notebooks 01-03
├── 02-analysis/           # Part 2: Stress derivation, parcel tracking, figures
│   ├── blue_ice_tools.py  # Utility library
│   └── notebooks 01-03
├── environment.yml
├── data/                  # Placeholder. Filled by downloads from Zenodo or processing
├── media/                 # Placeholder. Filled by Part 2 notebooks
└── README.md              # You are here
```

## Citation

If you use this code or data, please cite both the paper and this repository:

> **Paper:** [tbd]
> **Data:**  [tbf - Zenodo DOI]
> **Code:** [tbd - Zenodo DOI]