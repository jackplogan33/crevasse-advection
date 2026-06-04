import numpy as np
import xarray as xr
import sys
import glob
import json
import argparse

from joblib import Parallel, delayed
from pathlib import Path
from scipy.ndimage import median_filter
from skimage.restoration.inpaint import inpaint_biharmonic

########### STACKING ##########################################################
def resave_offset_nc(filebase, xmin, xmax, nx, ymin, ymax, ny, output_dir='.', chunk_size=34):
    """
    Read per-pair netCDF offset files produced by run_isce2.py and stack them
    into a single dataset along a 'mid_date' time dimension.

    Offsets are normalised to pixels/day by dividing by the number of days
    between the two acquisition dates encoded in the directory name
    (e.g. 20200101-20200113).

    Parameters
    ----------
    filebase : str
        Filename to glob for inside each 20*-20*/merged/ directory.
        If 'dense_offsets.nc', band_data is unpacked into 'azimuth'/'range'
        and divided by the temporal baseline. Otherwise the file is passed
        through unchanged (e.g. SNR, COV).
    xmin, xmax, nx : float, float, int
        Target x-axis extent and number of grid points.
    ymin, ymax, ny : float, float, int
        Target y-axis extent and number of grid points.
    chunk_size : int
        Number of offsets to merge in each intermediate batch before the
        final merge, to avoid memory spikes on large stacks.
    """
    # . . Get all offsets in local directory
    offset_dirs = sorted(glob.glob(f'20*-20*/merged/{filebase}'))
    if not offset_dirs:
        raise FileNotFoundError(         
            f"No directories matching '20*-20*/merged/{filebase}' found. "
            "Are you running from the correct working directory?"
        )

    x_grid = np.linspace(xmin, xmax, nx)
    y_grid = np.linspace(ymin, ymax, ny)

    # . . Read each offset and append
    offsets = []
    for i, offset_dir in enumerate(sorted(offset_dirs)):
        print(f"  Opening offset {i+1}/{len(offset_dirs)}: {offset_dir}")

        # . . Calculate mid date of S1 pass
        # Parse acquisition dates from directory name (format: YYYYMMDD-YYYYMMDD)
        date_pair = Path(offset_dir).parts[0]          # e.g. '20200101-20200113'
        start_str, end_str = date_pair.split('-')
        start = np.datetime64(f'{start_str[:4]}-{start_str[4:6]}-{start_str[6:]}')
        end   = np.datetime64(f'{end_str[:4]}-{end_str[4:6]}-{end_str[6:]}')
        days  = (end - start).astype('int64')          # temporal baseline in days
        mid_date = start + (end - start) // 2

        # . . Open dataset and interpolate to grid
        ds = xr.open_dataset(offset_dir, engine='rasterio')
        ds_interp = ds.interp(
            coords={'x':x_grid, 'y':y_grid},
            method='nearest'
        ).expand_dims({'mid_date':[mid_date]})

        # . . If offsets, unpack bands and divide by time
        if filebase == 'dense_offsets.nc':
            # band 1 = azimuth (along-track), band 2 = range (across-track)
            ds_new = ds_interp.drop_vars(['band_data', 'band'])
            ds_new['azimuth'] = ds_interp['band_data'].sel(band=1).drop_vars('band')
            ds_new['range']   = ds_interp['band_data'].sel(band=2).drop_vars('band')

            # Replace the ISCE2 nodata sentinel (-1e4) with NaN, then normalise
            # by the temporal baseline to get pixels per day
            offsets.append(ds_new.where(ds_new != -1e4).assign(
                azimuth=lambda d: d['azimuth'] / days,
                range=lambda d: d['range'] / days,
            ))
        else:
            offsets.append(ds_interp)

    # . . Merge in batches to keep peak memory manageable
    print(f'  Merging {len(offsets)} offsets in batches of {chunk_size}...')
    batches = [
        xr.merge(offsets[i:i+chunk_size], compat='no_conflicts', join='outer')
        for i in range(0, len(offsets), chunk_size)
    ]
    offset_ds = xr.merge(batches).assign_attrs({
        'grid_resolution_meters': int(round((xmax - xmin) / (nx - 1))),
        'CRS': 'EPSG:3031',
        'temporal_units': 'pixels/day'
    })

    # . . Write to disk
    outpath = Path(output_dir) / filebase
    print(f'  Writing to {outpath}')
    offset_ds.to_netcdf(outpath)


########### CLEANING ##########################################################
def offset_cleaning(
    ds: xr.Dataset,
    window_size: int = 80,
    mad_mult: float = 15.0,
    az_bounds: tuple[float, float] = (-0.125, 0.55),
    rg_bounds: tuple[float, float] = (-0.58, 1.75),
    filter_size: tuple[int, int, int] = (3, 7, 7), 
) -> xr.Dataset:
    """
    Mask and inpaint outliers in a stacked offset dataset.

    Two masking strategies are combined with a logical OR:
      1. Hard threshold — removes offsets outside physically plausible bounds.
         Defaults correspond to ~1.5 m/day in azimuth and ~4.0 m/day in range
         for a 12-day Sentinel-1 pair.
      2. MAD outlier filter — for each timestep, computes the median absolute
         deviation of (data - local_median) over the spatial domain and masks
         pixels exceeding `mad_mult` times that value.

    Masked pixels are reconstructed using biharmonic inpainting (∇⁴u = 0),
    applied independently to each 2-D time slice.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with variables 'azimuth' and 'range', dims (mid_date, y, x).
    window_size : int
        Side length (pixels) of the median filter kernel (~4 km at 50 m res).
    mad_mult : float
        MAD threshold multiplier. Lower = more aggressive masking.
    az_bounds : (float, float)
        (min, max) hard threshold for azimuth offsets in pixels/day.
    rg_bounds : (float, float)
        (min, max) hard threshold for range offsets in pixels/day.
    filter_size : (int, int, int)
        (time, y, x) median filter size to compute after all cleaning steps

    Returns
    -------
    xr.Dataset
        Cleaned dataset with the same structure as `ds`.
    """
    # . . Hard Threshold mask
    print('  Applying treshold mask...')
    thresh_mask = xr.Dataset({
        'azimuth': (ds['azimuth'] < az_bounds[0]) | (ds['azimuth'] > az_bounds[1]),
        'range':   (ds['range']   < rg_bounds[0]) | (ds['range']   > rg_bounds[1])
    })

    # . . Local MAD mask
    print('  Computing local median (large kernel)...')
    # Braodcast over time axis
    med = xr.apply_ufunc(
        median_filter,
        ds,
        kwargs={'size': window_size, 'axes':[1, 2]}
    )
    
    print('  Applying MAD mask...')
    # MAD is computed per-timestep over the full spatial domain
    diff = np.abs(ds - med)
    mad = diff.median(['x', 'y'], skipna=True)
    mad_mask = diff > (mad_mult * mad)
    
    # . . Combine Masks
    mask = mad_mask + thresh_mask

    print('  Inpainting NaNs')
    # . . Biharmonic inpainting across time dimension
    ds_inpaint = xr.apply_ufunc(
        inpaint_slice,
        ds,
        mask,
        input_core_dims=[['y', 'x'], ['y', 'x']],
        output_core_dims=[['y', 'x']],
        vectorize=True,   # loops over mid_date, calls inpaint_biharmonic once per timestep
    )

    # . . 3D median filter
    print('  Applying 3D median filter...')
    ds_filtered = xr.apply_ufunc(
        median_filter,
        ds_inpaint,
        kwargs={'size': filter_size}
    )
    return ds_filtered


def inpaint_slice(data, mask):
    return inpaint_biharmonic(data.copy(), mask.copy())    
    

########### COORDINATE CONVERSION #############################################
def convert_coords(ds, az_res=14.1, rg_res=2.3, inc_angle=38.3, heading=0.0097):
    """
    Convert pixel offsets (pixels/day) to surface velocity (m/year) in
    a map-projected coordinate system (EPSG:3031).

    Steps
    -----
    1. Scale azimuth by azimuth pixel spacing and range by the ground-range
       pixel spacing (range spacing / sin(incidence angle)).
    2. Rotate from SAR geometry (azimuth/range) into map coordinates (x/y)
       using the satellite heading angle.
    3. Convert m/day to m/year.

    Parameters
    ----------
    ds : xr.Dataset
        Cleaned offsets with variables 'azimuth' and 'range' in pixels/day.
    az_res : float
        Azimuth pixel spacing in metres (Sentinel-1 IW: ~14.1 m).
    rg_res : float
        Slant-range pixel spacing in metres (Sentinel-1 IW: ~2.3 m).
    inc_angle : float
        Radar incidence angle in degrees (used to project slant → ground range).
    heading : float
        Satellite heading in radians relative to the EPSG:3031 y-axis.

    Returns
    -------
    xr.Dataset
        Dataset with variables 'vx', 'vy' (m/year) and 'vv' (speed, m/year).
    """
    # Convert pixel offsets to meters/day
    az_m = ds['azimuth'] * az_res
    rg_m = ds['range']   * rg_res / np.sin(np.radians(inc_angle))  # Ground-range (only horizontal)

    # Rotate from SAR geometry to map-projected x/y
    vx = (np.cos(heading) * rg_m - np.sin(heading) * az_m) * 365
    vy = (np.cos(heading) * az_m + np.sin(heading) * rg_m) * 365
    
    # Unit conversion to velocity, saved in dataset
    return xr.Dataset({
        'vx': vx.assign_attrs({
            "long_name":    "Surface velocity in x-direction",
            "standard_name":"land_ice_surface_x_velocity",
            "units":        "m yr-1",
            "grid_mapping": "spatial_ref",
        }),
        'vy': vy.assign_attrs({
            "long_name":    "Surface velocity in y-direction",
            "standard_name":"land_ice_surface_y_velocity",
            "units":        "m yr-1",
            "grid_mapping": "spatial_ref",
        }),
        'vv': np.sqrt(vx**2 + vy**2).assign_attrs({
            "long_name":    "Surface ice speed",
            "standard_name":"land_ice_surface_velocity",
            "units":        "m yr-1",
            "grid_mapping": "spatial_ref",
        }),
    })


########### CLI ###############################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description='Postprocess Sentinel-1 dense offsets: stack, clean, convert to velocity.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('config', type=str,
                        help='Path to a JSON config file.')
    return parser.parse_args()
 
    
def load_config(path):
    """Load a JSON config and return the postprocessing sub-section."""
    with open(path) as f:
        cfg = json.load(f)
    return cfg

if __name__ == '__main__':
    args = parse_args()
    cfg = load_config(args.config)
    postprocessing = cfg.get('postprocessing', {})

    # . . Stack Params
    stack_cfg     = postprocessing.get('stack', {})
    do_stack      = stack_cfg.get('enabled', False)
    stack_out_dir = stack_cfg.get('output_dir', '.')
    grid          = stack_cfg.get('grid', {})


    # . . Cleaning params
    offsets_file    = postprocessing.get('offsets', None)
    output_file     = postprocessing.get('output', 'filt_dense_offsets.nc')
    cleaning        = postprocessing.get('cleaning', {})
    window          = int(cleaning.get('window_size', 80))
    mad_mult        = float(cleaning.get('mad_mult', 15.0))
    az_bounds       = tuple(cleaning.get('az_bounds', [-0.125, 0.55]))
    rg_bounds       = tuple(cleaning.get('rg_bounds', [-0.58,  1.75]))
    temporal_filter = tuple(cleaning.get('temporal_filter', [3, 7, 7]))

    # . . Velocity params
    vel_cfg     = postprocessing.get('velocity', {})
    do_velocity = vel_cfg.get('do_velocity', True)
    vel_file    = vel_cfg.get('output', 'velocity.nc')
    az_res      = vel_cfg.get('az_res',    14.1)
    rg_res      = vel_cfg.get('rg_res',    2.3)
    inc_angle   = vel_cfg.get('inc_angle', 38.3)
    heading     = vel_cfg.get('heading',   0.0097)
        
    if do_stack:
        print('Stacking offset pairs...')
        fileroots = ['dense_offsets.nc', 'dense_offsets_snr.nc', 'dense_offsets_cov.nc']
        for fr in fileroots:
            resave_offset_nc(
                fr,
                output_dir=stack_out_dir,
                xmin=grid.get('xmin', 1.345e6+50),
                xmax=grid.get('xmax', 1.43e6-50),
                nx=grid.get('nx', 1699),
                ymin=grid.get('ymin', 1.646e6+50),
                ymax=grid.get('ymax', 1.76e6-50),
                ny=grid.get('ny', 2279),
            )
        # Derive offsets path from output_dir
        offsets_file = str(Path(stack_out_dir) / 'dense_offsets.nc')

    else:
        if offsets_file is None:
            print("Error: stack.enabled is false but no 'offsets' path was provided in config.")
            sys.exit(1)
    
    # . . Read in dataset of all offsets
    print(f'Reading offsets from {offsets_file}...')
    ds = xr.open_dataset(offsets_file)

    print('Cleaning offsets...')
    ds_clean = offset_cleaning(
        ds, 
        window_size=window, 
        mad_mult=mad_mult,
        az_bounds=az_bounds, 
        rg_bounds=rg_bounds,
        filter_size=temporal_filter
    )    
    print(f'Writing cleaned offsets to {output_file}...')
    ds_clean.to_netcdf(output_file)
    
    if do_velocity:
        print(f'Converting to velocity {vel_file}...')
        vel_ds = convert_coords(
            ds_clean, 
            az_res=az_res, 
            rg_res=rg_res,
            inc_angle=inc_angle, 
            heading=heading
        )

        vel_ds = vel_ds.assign_attrs({
            "description": "Sentinel-1 dense offset velocities processed using ISCE2 topsApp. "
                           "Units are meters per year in EPSG:3031.",
            "ampcor_window_width":           cfg["ampcor"]["window_width"],
            "ampcor_window_height":          cfg["ampcor"]["window_height"],
            "ampcor_skip_width":             cfg["ampcor"]["skip_width"],
            "ampcor_skip_height":            cfg["ampcor"]["skip_height"],
            "ampcor_search_width_per_pass":  cfg["ampcor"]["search_width_per_pass"],
            "ampcor_search_height_per_pass": cfg["ampcor"]["search_height_per_pass"],
            "filter_window_size":            window,
            "filter_mad_mult":               mad_mult,
            "filter_az_bounds":              str(az_bounds),
            "filter_rg_bounds":              str(rg_bounds),
            "filter_temporal":               str(temporal_filter),
        }).to_netcdf(vel_file)