# Part 1: Sentinel-1 Velocity Processing

This directory contains a cloud-native pixel offset tracking pipeline built on [ISCE2](https://github.com/isce-framework/isce2). 
It ingests Sentinel-1 SAR image pairs and produces surface velocity fields via amplitude tracking. 
The workflow is designed to be reproducible for Shirase Glacier specifically, but is straightforward to adapt for any region with adequate Sentinel-1 coverage.

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
```
quay.io/jackplogan/isce-image:0e337edc16fc
```

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
A returned path confirms the install succeeded. When running notebooks, make sure the **ISCE2** kernel is selected.

---

## Usage
 
This directory contains two scripts and three notebooks that together take you from
raw Sentinel-1 scenes to a calibrated surface velocity timeseries.
 
The intended workflow is:
 
1. **Run `01-setup.ipynb`** to download the DEM and set up the file structure.
2. **Edit `config.json`** to define your scene, processing parameters, and output paths.
3. **Run `run_isce2.py`** to download SAFE files and compute per-pair pixel offsets.
4. **Run `postprocess.py`** to stack, clean, and convert offsets to velocity.
Notebooks `02-run_topsapp.ipynb` and `03-postprocess.ipynb` walk through steps 3 and 4
interactively and are the recommended starting point before using the CLI scripts directly.
Once you're familiar with the parameters, the scripts can be run end-to-end without
opening a notebook.
 
All processing parameters for both scripts are controlled through a single `config.json`
file, documented in the [Configuration](#configuration) section.

---

## Notebooks

| Notebook | Purpose | Required? |
|----------|---------|-----------|
| `01-setup.ipynb` | Download DEM, set up file structure | Yes |
| `02-run_topsapp.ipynb` | Download inputs, configure and run ISCE2 for a single offset | Recommended |
| `03-postprocess.ipynb` | Convert pixel offsets to ground velocity in xy-coordinates | Recommended |

Notebooks 02 and 03 each have a corresponding Python script (`run_isce2.py`, 
`postprocess.py`) that can be run from the command line once you're familiar with 
the parameters.
 
---
 
## Running ISCE2

`run_isce2.py` automates the full ISCE2 offset processing pipeline from data download through
per-pair topsApp runs.

### Usage

```bash
python run_isce2.py
```

The script is interactive: it prompts for the number of SAFE files to download,
your Earthdata credentials, and a start date for the scene search. All other
parameters are read from `config.json` in the current directory.

### What it does

1. **Downloads SAFE files** via the ASF Python API, filtered by the `aoi`, `frames`,
   and `polarization` defined in config. Files are saved to `./SAFE/`.
2. **Downloads orbit files** for each SAFE file using `fetchOrbit.py` (included with ISCE2).
3. **Iterates through consecutive image pairs** (reference → secondary), and for each:
   - Creates a pair subdirectory named `YYYYMMDD-YYYYMMDD/`
   - Writes `reference.xml`, `secondary.xml`, and `topsApp.xml` from config parameters
   - Scales the Ampcor search window by the number of 12-day Sentinel-1 passes separating
     the pair (e.g. a 24-day pair gets twice the search range of a 12-day pair)
   - Runs `topsApp.py` and streams output to `topsapp_processing.log`
   - On success: crops the dense offset, SNR, and COV outputs to the extent defined
     in `crop`, saves them as netCDF, and removes large intermediates to save disk space
   - Removes the reference SAFE file after a successful pair to manage storage
4. **Removes the final SAFE file** once all pairs are complete.
 
### File structure after a run
 
```
offsets/
├── SAFE/                        # Downloaded .zip files (removed as pairs complete)
├── orbits/                      # Orbit files
├── topsapp_processing.log
├── 20200101-20200113/
│   ├── merged/
│   │   ├── dense_offsets.nc     # Cropped azimuth/range offsets
│   │   ├── dense_offsets_snr.nc
│   │   ├── dense_offsets_cov.nc
│   ├── reference.xml
│   ├── secondary.xml
│   └── topsApp.xml
├── 20200113-20200125/
│   └── ...
└── ...
```

---
 
## Offset Postprocessing

`postprocess.py` stacks per-pair offset files, removes outliers, and converts to surface velocity.
Driven entirely by the `postprocessing` block of `config.json`.

### Usage
 
```bash
python postprocess.py config.json
```
 
### What it does
 
#### Step 1 — Stacking (optional)
 
If `postprocessing.stack.enabled` is `true`, the script finds all
`20*-20*/merged/dense_offsets.nc` files in the working directory and merges them
into a single dataset along a `mid_date` time dimension. It does the same for the
SNR and COV files. Each offset is:
 
- Interpolated onto the common spatial grid defined in `stack.grid`
- Normalized by the temporal baseline (days between acquisitions) to give pixels/day
- Tagged with the mid-date of the acquisition pair
The three stacked files are written to `stack.output_dir`. The path to
`dense_offsets.nc` is then derived automatically for the next step.

If `stack.enabled` is `false`, the script reads directly from the path given
in `postprocessing.offsets`.
 
#### Step 2 — Cleaning
 
Outliers are removed using two complementary strategies applied with a logical OR:
 
- **Hard threshold:** pixels outside `az_bounds` or `rg_bounds` are masked immediately.
- **MAD filter:** for each timestep, a large spatial median is computed
  (`window_size × window_size`), and pixels where the absolute deviation exceeds
  `mad_mult × MAD` are masked.
Masked pixels are reconstructed using biharmonic inpainting (∇⁴u = 0) applied
independently to each 2D time slice. A final 3D median filter (`temporal_filter`)
is applied across the full stack. The cleaned dataset is written to
`postprocessing.output`.
 
#### Step 3 — Velocity conversion (optional)
 
If `postprocessing.velocity.enabled` is `true`, pixel offsets (pixels/day) are
converted to surface velocity (m/year) in EPSG:3031 map coordinates:
 
1. Azimuth offsets are scaled by `az_res`; range offsets are scaled by
   `rg_res / sin(inc_angle)` to project from slant range to ground range.
2. The azimuth/range vector is rotated into map x/y using the satellite `heading`.
3. m/day is converted to m/year.
Output variables are `vx`, `vy`, and `vv` (speed), written to `velocity.output`.

### Output files

| File | Contents |
|------|----------|
| `{stack.output_dir}/dense_offsets.nc` | Stacked azimuth and range offsets (pixels/day) |
| `{stack.output_dir}/dense_offsets_snr.nc` | Stacked signal-to-noise ratio |
| `{stack.output_dir}/dense_offsets_cov.nc` | Stacked covariance |
| `{postprocessing.output}` | Cleaned and filtered offsets |
| `{velocity.output}` | Surface velocity in m/year (vx, vy, vv) |

---

## Adapting for a Different Region

The key parameters to modify are:

- **DEM extent:** defined in `01-setup.ipynb`
- **Scene selection and date range:** `scene.aoi`, `scene.frames` in `config.json`, the start date prompted by `run_isce2.py`,
  and the number of scenes to download.
- **Crop extent:** `crop.xmin/xmax/ymin/ymax` in `config.json`
- **Stack grid:** `postprocessing.stack.grid` in `config.json` — should match the crop
  extent with a 50 m inset on each edge
- **Geocoding and filtering parameters:** `ampcor`, `cleaning`, and `velocity` blocks
  in `config.json`

The pipeline assumes Sentinel-1 IW mode data available via ASF. Other SAR sensors 
would require modifications to the topsApp XML configuration.

---

## Configuration
 
Both scripts are driven by a single `config.json` file. A full annotated example is
shown below, followed by a description of each block.
 
```json
{
    "scene": {
        "aoi": "POLYGON((38.0336 -69.7358,38.0336 -70.4952,39.6985 -70.4952,39.6985 -69.7358,38.0336 -69.7358))",
        "frames": [830, 834, 936, 938, 939],
        "dem_file": "./dem/REMA_10m_shirase.bil",
        "swaths": [2],
        "region_of_interest": [-70.45441464, -69.88490745, 38.29553816, 39.83817435],
        "polarization": "hh"
    },
    "insar": {
        "do_interferogram": false,
        "do_esd": false,
        "do_unwrap": false,
        "do_unwrap_2stage": false,
        "do_ionosphere": false,
        "geocode_list": []
    },
    "ampcor": {
        "do_denseoffsets": true,
        "window_width": 256,
        "window_height": 64,
        "skip_width": 44,
        "skip_height": 8,
        "search_width_per_pass": 30,
        "search_height_per_pass": 10
    },
    "crop": {
        "xmin": 1345000,
        "xmax": 1430000,
        "ymin": 1646000,
        "ymax": 1760000,
        "output_epsg": 3031
    },
    "compute": {
        "omp_num_threads": 8
    },
    "postprocessing": {
        "stack": {
            "enabled": true,
            "output_dir": "../../data",
            "grid": {
                "xmin": 1345050,
                "xmax": 1429950,
                "nx": 1699,
                "ymin": 1646050,
                "ymax": 1759950,
                "ny": 2279
            }
        },
        "offsets": "../../data/dense_offsets.nc",
        "output": "../../data/filt_dense_offsets.nc",
        "cleaning": {
            "window_size": 31,
            "mad_mult": 10.0,
            "az_bounds": [-0.125, 0.55],
            "rg_bounds": [-0.58, 1.75],
            "temporal_filter": [3, 7, 7]
        },
        "velocity": {
            "enabled": true,
            "output": "../../data/velocity/velocity-timeseries.nc",
            "az_res": 14.1,
            "rg_res": 2.3,
            "inc_angle": 38.3,
            "heading": 0.0097
        }
    }
}
```
 
### `scene`
 
Controls Sentinel-1 scene selection and acquisition geometry.

| Key | Type | Description |
|-----|------|-------------|
| `aoi` | WKT string | Area of interest in WGS84 lon/lat, used to search the ASF catalog |
| `frames` | list of int | Sentinel-1 frame numbers to restrict the search |
| `dem_file` | string | Path to the DEM (`.bil` format). Must have sidecar `.xml` and `.vrt` files |
| `swaths` | list of int | IW subswaths to process (1, 2, and/or 3) |
| `region_of_interest` | list of float | Bounding box `[lat_min, lat_max, lon_min, lon_max]` passed to topsApp |
| `polarization` | string | SAR polarization channel (e.g. `"hh"`, `"vv"`) |
 
### `insar`
 
Toggles InSAR processing steps in topsApp. For pixel offset tracking, all of these
should be `false` and `geocode_list` should be empty.
 
### `ampcor`
 
Controls the ISCE2 amplitude cross-correlation (Ampcor) dense offset estimation.
 
| Key | Type | Description |
|-----|------|-------------|
| `do_denseoffsets` | bool | Must be `true` to run offset tracking |
| `window_width` / `window_height` | int | Correlation window size in pixels. Larger = more stable but coarser |
| `skip_width` / `skip_height` | int | Step size between correlation windows. Controls output resolution |
| `search_width_per_pass` / `search_height_per_pass` | int | Per-pass search range. Scaled by the number of Sentinel-1 passes separating the image pair, so a 24-day pair searches twice as wide as a 12-day pair |
 
### `crop`
 
Defines the output spatial extent and coordinate reference system after geocoding.
Coordinates are in metres in the projection defined by `output_epsg` (EPSG:3031,
Antarctic Polar Stereographic, by default).
 
### `compute`
 
| Key | Type | Description |
|-----|------|-------------|
| `omp_num_threads` | int | Number of OpenMP threads for Ampcor. Set to the number of available cores |
 
### `postprocessing`
 
Controls stacking, cleaning, and velocity conversion. This block is used exclusively
by `postprocess.py` and is ignored by `run_isce2.py`.
 
#### `stack`
 
Determines whether raw per-pair offset files are stacked into a single dataset.
 
| Key | Type | Description |
|-----|------|-------------|
| `enabled` | bool | If `true`, stack all per-pair offsets found in the working directory |
| `output_dir` | string | Directory where the stacked `.nc` files are written. Required when `enabled: true` |
| `grid` | object | Target spatial grid for interpolation (see below). Required when `enabled: true` |
 
`grid` sub-keys:
 
| Key | Description |
|-----|-------------|
| `xmin`, `xmax`, `nx` | x-axis extent (metres) and number of grid points |
| `ymin`, `ymax`, `ny` | y-axis extent (metres) and number of grid points |
 
**Two modes of operation:**
 
*You need to stack:* set `stack.enabled: true` and configure `output_dir` and `grid`.
The stacked `dense_offsets.nc` will be written to `output_dir` automatically — you do
not need to set `offsets` in this case.
 
```json
"stack": { "enabled": true, "output_dir": "../data", "grid": { ... } },
"offsets": null
```
 
*You already have a stack:* set `stack.enabled: false` and point `offsets` at your
existing file. The `output_dir` and `grid` keys are ignored.
 
```json
"stack": { "enabled": false },
"offsets": "../data/dense_offsets.nc"
```
 
#### `offsets`
 
Path to the stacked offset dataset. Only required when `stack.enabled: false`.
When stacking is enabled, this path is derived automatically from `stack.output_dir`.
 
#### `output`
 
Path where the cleaned, filtered offset dataset is written.
 
#### `cleaning`
 
Outlier masking and inpainting parameters.
 
| Key | Type | Description |
|-----|------|-------------|
| `window_size` | int | Spatial median filter kernel size (pixels) used for local MAD computation |
| `mad_mult` | float | MAD threshold multiplier. Lower values = more aggressive masking |
| `az_bounds` | [min, max] | Hard azimuth offset bounds (pixels/day). Values outside are masked |
| `rg_bounds` | [min, max] | Hard range offset bounds (pixels/day). Values outside are masked |
| `temporal_filter` | [t, y, x] | 3D median filter kernel applied after inpainting |
 
#### `velocity`
 
Controls conversion from pixel offsets to surface velocity.
 
| Key | Type | Description |
|-----|------|-------------|
| `enabled` | bool | If `true`, write a velocity dataset after cleaning |
| `output` | string | Path for the velocity output file |
| `az_res` | float | Azimuth pixel spacing in metres (Sentinel-1 IW: 14.1 m) |
| `rg_res` | float | Slant-range pixel spacing in metres (Sentinel-1 IW: 2.3 m) |
| `inc_angle` | float | Radar incidence angle in degrees |
| `heading` | float | Satellite heading in radians relative to the EPSG:3031 y-axis |