# Part 1: Sentinel-1 Velocity Processing

This directory contains a cloud-native pixel offset tracking pipeline built on [ISCE2](https://github.com/isce-framework/isce2). 
It ingests Sentinel-1 SAR image pairs and produces surface velocity fields via amplitude tracking. 
The workflow is designed to be reproducible for Shirase Glacier specifically, but is straightforward to adapt for any region with adequate Sentinel-1 coverage.

## Notebooks

| Notebook | Purpose | Required? |
|----------|---------|-----------|
| `01-setup.ipynb` | Download DEM, set up file structure | Yes |
| `02-run_topsapp.ipynb` | Download inputs, configure and run ISCE2 for a single offset | Recommended |
| `03-postprocess.ipynb` | Convert pixel offsets to ground velocity in xy-coordinates | Recommended |

Notebooks 02 and 03 each have a corresponding Python script (`run_isce2.py`, 
`postprocess.py`) that can be run from the command line once you're familiar with 
the parameters.

## Installation

Choose the option that matches your computing environment:

| Option | Environment |
|--------|-------------|
| [Custom Image](#custom-image) | Cloud environment where you can bring your own image *(easiest)* |
| [Manual Install](#manual-install) | Cloud environment with a fixed image, or local machine |

---

### Custom Image

A pre-built JupyterHub image with all dependencies and environment variables configured
is available on quay.io:
```quay.io/jackplogan/isce-image:0e337edc16fc ```

**Setup:**

1. Clone the repository:
```bash
git clone https://github.com/jackplogan33/crevasse-advection/
```

That's it! You can now run all the notebooks.

> [!NOTE]
> This image includes `jupyter-keepalive`, which keeps your JupyterHub session
> active during long processing runs (offset processing can take ~20 hours).

---

### Manual Install

This option requires setting several environment variables that can conflict with
existing configurations. Additional debugging may be needed depending on your setup.

1. Clone the repository:
```bash
git clone https://github.com/jackplogan33/crevasse-advection/
cd crevasse-advection
```

2. Create the conda environment:
```bash
conda env create -f environment.yml
```

3. Activate the environment and run the setup script:
```bash
conda activate isce2
./setup.sh
```

> [!IMPORTANT]
> Activate the environment *before* running `setup.sh`. The script uses `$CONDA_PREFIX`
> to locate the correct environment and edit its activation/deactivation scripts,
> which makes the ISCE2 command-line tools accessible on your PATH.

4. Reactivate to apply the changes:
```bash
conda deactivate
conda activate isce2
```

**Verify the install:**
```bash
which topsApp.py
```
A returned path confirms the install succeeded. When running notebooks, make sure
the **ISCE2** kernel is selected.

## Adapting for a Different Region

The key parameters to modify are:

- **DEM extent:** defined in `01-setup.ipynb`
- **Scene selection and date range:** defined in `02-run_topsapp.ipynb`
- **Geocoding and filtering parameters:** documented in `03-postprocess.ipynb` and `postprocess.py`

The pipeline assumes Sentinel-1 IW mode data available via ASF. Other SAR sensors 
would require modifications to the topsApp XML configuration.
