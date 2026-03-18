# ISCE2 Processing and Crevasse Advection

Using the ISCE2 interferometric processor, we perform pixel offset tracking to determine the flow of glacial ice from Sentinel-1 Synthetic Aperture Radar (SAR).
This workflow is cloud native and is designed to process multiple offsets in series, exploiting the amplitude of a SAR image. 
It can easily be applied to create your own surface velocity timeseries anywhere with adequate Sentinel-1 coverage.

More information will be detailed in our forthcoming publication.

## Installation

There are two installation options depending on your computing environment and permissions.
1. You are in a cloud computing environment and *can* bring your own image (easiest)

2. You are in a cloud computing environment and *cannot* bring your own image, **or** are running locally (less easy)

### Custom Image

I built a JupyterHub image hosted in quay.io.
If your cloud computing environment allows you to pull your own image, you can do so using the following link:
`quay.io/jackplogan/isce-image:b8a558f43a10`

This image has all of the environment variable set, as far as I have tested.

1. Clone repository:
```bash
git clone https://github.com/jackplogan33/crevasse-advection/
```
That's it. Now you can run all the notebooks.

### Running locally or a community image
This process is a little bit more involved. 
There are many conflicting environment variables that need to be specified for everything to run properly.
Depending on your machine's configuration, these steps may not work and will require additional debugging.

1. Clone repository:
```bash
git clone https://github.com/jackplogan33/crevasse-advection/
cd ./crevasse-advection
```

2. Install conda environment:
```bash
conda env create -f environment.yml
```

3. Set environment variables and create a kernel for use in Jupyter Lab:
```bash
conda activate isce2
./setup.sh
source ~/.bashrc
```

> Make sure you activate the environment first. This ensures that the `$CONDA_PREFIX` is pointing to the isce environment, not the default environment
>
> We have had issues with the proj variables during setup in community-maintained images.
> To ensure the right paths are set, they need to be explicity written.

In theory, this should work. To test, you can run:
```bash
which topsApp.py
```
If it returns a path, it worked. 
Now you can run the notebooks. 
When running the notebooks, make sure you have selected the "ISCE2" kernel.


## Usage
We have written a series of notebooks to explain the setup, run, and postprocessing steps. 

`01-ISCE2_setup.ipynb`: <br>
Download the DEM for the region of interest and ensure the file structure is properly setup to run ISCE2 in a contained file structure.

`02-run_topsapp.ipynb` <br>
Download the necessary input files for a single offset and write the XML files for topsApp. 
This notebook explains each of the steps occuring in the run script `./offsets/run_isce2.py`.

`03-postprocess_offsets` <br>
Turn the offsets from azimuth-range displacements in pixels to ground velocity in xy-coordinates.