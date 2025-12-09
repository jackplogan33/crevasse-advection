# ISCE2 Processing and Crevasse Advection

Using the ISCE2 interferometric processor, we perform pixel offset tracking to determine the flow of glacial ice from Sentinel-1 Synthetic Aperture Radar (SAR).
This workflow is cloud native and is designed to process multiple offsets in series, exploiting the amplitude of a SAR image. 
It can easily be applied to create your own surface velocity timeseries anywhere with adequate Sentinel-1 coverage.

More information will be detailed in our forthcoming publication.

## Installation

This repository requires creating two separate conda environments. 
The first includes all required packages to run ISCE2. 
The second includes all packages required for the input files and output postprocessing. 

1. Clone repository:
```bash
git clone https://github.com/jackplogan33/crevasse-advection/
cd ./crevasse-advection
```

2. ISCE Environment:
Run the following commands below in terminal:
```bash
conda create --name isce2 --file isce_requirements.txt
conda activate isce2
conda install isce2 boto3 jupyter conda-build ipykernel
```

3. Processing Environment:
```bash
conda env create -f environment.yml
```

4. Create a kernel for use in Jupyter Lab
```bash
python -m ipykernel install --user --name crevasse-advection --display-name "IPython - crevasse-advection"
```

## Usage

We have written a series of notebooks to explain the setup, run, and postprocessing steps. 

`01-ISCE2_setup.ipynb`: <br>
Download the DEM for the region of interest and ensure the file structure is properly setup to run ISCE2 in a contained file structure.

`02-run_topsapp.ipynb` <br>
Download the necessary input files for a single offset and write the XML files for topsApp. 
This notebook explains each of the steps occuring in the run script `./offsets/run_isce2.py`.

`03-postprocess_offsets` <br>
Turn the offsets from azimuth-range displacements in pixels to ground velocity in xy-coordinates.
