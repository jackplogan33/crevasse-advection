import numpy as np
import xarray as xr
from scipy.ndimage import median_filter, laplace
import sys

def resave_offset_nc(filebase):
    """Read the netCDF files produced by run_isce2.py and stack multiple offsets.
    Offsets are merged along a new dimension names 'mid_date', which is the midpoint between the two dates in the folder name."""
    # . . Interpolation params
    xmin, xmax, nx = 1.345e6+50, 1.43e6-50, 1699
    ymin, ymax, ny = 1.646e6+50, 1.76e6-50, 2279

    # . . Get all offsets in local directory
    offset_dirs = glob.glob('20*-20*/merged/'+filebase)

    # . . Read each offset and append
    offsets = []
    for i, offset_dir in enumerate(sorted(offset_dirs)):
        print(f"Opening offset {i+1}/{len(offset_dirs)}")

        # . . Calculate mid date of S1 pass
        # Get start and end dates from root directory
        month_dir = offset_dir.split('/')[0]
        start_str, end_str = month_dir.split('-')
        start = np.datetime64(f'{start_str[:4]}-{start_str[4:6]}-{start_str[6:]}')
        end = np.datetime64(f'{end_str[:4]}-{end_str[4:6]}-{end_str[6:]}')
        # Calculate number of days between pass (typically 12)
        days = (end - start).astype('int64')
        
        # Calculate middle date
        mid_date = start + (end - start) / 2

        # . . Open dataset and interpolate to grid
        ds = xr.open_dataset(offset_dir, engine='rasterio')
        ds_interp = ds.interp(  #### Do I need to do this?
            coords={
            'x':np.linspace(xmin, xmax, nx),
            'y':np.linspace(ymin, ymax, ny)
            },
            method='nearest'
        )
        
        # Make time dimension
        ds_interp = ds_interp.expand_dims({'mid_date':[mid_date]})

        # . . If offsets, unpack bands and divide by time
        if filebase == 'dense_offsets.nc':
            # Unpack bands to variable, remove 'band_data'
            ds_new = ds_interp.drop_vars(['band_data', 'band'])
            ds_new['azimuth'] = ds_interp['band_data'].sel(band=1).drop_vars('band')
            ds_new['range']   = ds_interp['band_data'].sel(band=2).drop_vars('band')

            # Return offets as meters per day
            offsets.append(ds_new.where(ds_new == -1e4, ds_new / days))

        # . . Append interpolated values for others
        else:
            # Return SNR and COV unchanged
            offsets.append(ds_interp)

    # . . Merge offsets in quarters
    offsets_quarters = [xr.merge(offsets[i:i+34]) for i in range(0, 134, 34)]

    # . . Merge final 4 parts
    print('Merging Offsets')
    offset_ds = xr.merge(offsets_quarters).assign_attrs({
        'grid_resolution_meters':50,
        'CRS':'EPSG:3031'
    })

    # . . Write to disk
    print('Writing to disk')
    offset_ds.to_netcdf('./'+filebase)

def laplacian_infill(arr, mask, max_iter=2000, tol=1e-4, mode='nearest'):
    """"""
    filled = arr.copy()

    for i in range(max_iter):
        prev = filled.copy()
        lap = laplace(filled, mode=mode, axes=[-2, -1])
        filled[mask] += (0.2 * lap[mask])
        
        # Compute relative difference
        diff = np.linalg.norm((filled - prev)[mask]) / (np.linalg.norm(filled[mask]) + 1e-12)

        # Print convergence info (overwriting previous line)
        print(f"\rIteration {i+1}/{max_iter}, diff = {diff:.2e}", end="")
        sys.stdout.flush()

        if diff < tol:
            print(f"\nConverged after {i+1} iterations.")
            break
    else:
        print(f"\nMax iterations reached. Final diff = {diff:.2e}")

    return filled

def offset_cleaning(ds, window_size=80, mult=5):
    print('Masking by Treshold... ')
    # Make mask for Thresholding
    az_mask = (ds['azimuth'] > (.55)) | (ds['azimuth'] < (-.125))
    rg_mask = (ds['range'] > (1.75)) | (ds['range'] < (-.58))
    thresh_mask = xr.Dataset({
        'azimuth':az_mask,
        'range':rg_mask
    })
    
    print('Applying Large Median Filter...')
    # Compute large median filter over spatial dimensions (4 km square window)
    med = xr.apply_ufunc(
        median_filter,
        ds,
        kwargs={'size':window_size, 'axes':[-2, -1]}
    )
    
    print('Masking by MAD...')
    # Compute localized Median Absolute Deviation for each timestep
    diff = np.abs(ds - med)
    mad = diff.median(['x', 'y'], skipna=True)
    
    # Threshold at 5 times the MAD
    threshold = mult * mad
    mad_mask = (diff > threshold)
    
    # Combine MAD and threshold Masks
    mask = mad_mask + thresh_mask
    
    # Mask values, replace with median filter result
    ds_masked = ds.where(~mask, med)
    
    print('Harmonic Inpainting...')
    lap_kwargs = {'tol':1e-4, 'max_iter':2000, 'mode':'nearest'}
    # Start Infill of values to new dataset
    ds_filled = ds_masked.copy()
    for var in ds_filled.data_vars:
        ds_filled[var] = xr.apply_ufunc(
            laplacian_infill,
            ds_filled[var],
            mask[var],
            kwargs=lap_kwargs
        )

    return ds_filled

def convert_coords(ds):
    az = ds['azimuth']
    rg = ds['range']
    # Multiply offsets by resolution to get disp in meters
    az *= 14.1
    rg *= 2.3 / np.sin(np.radians(38.3))  # Ground-range (only horizontal)

    # S1 orbit heading with reference to EPSG:3031 grid
    heading = 0.0097

    # Coordinate transformation
    x_offsets = np.cos(heading) * rg - np.sin(heading) * az
    y_offsets = np.cos(heading) * az + np.sin(heading) * rg

    # Convert to meters / year
    vx = x_offsets * 365
    vy = y_offsets * 365
    
    # Unit conversion to velocity, saved in dataset
    return xr.Dataset({
        'vx':vx,
        'vy':vy,
        'vv':np.sqrt((vx ** 2) + (vy ** 2))
    })

if __name__ == '__main__':
    # . . Stack all offsets, interpolate sampling
    fileroots = ['dense_offsets.nc', 'dense_offsets_snr.nc', 'dense_offsets_cov.nc']
    for fr in fileroots:
        resave_offset_nc(fr)
    
    # . . Read in dataset of all offsets
    ds = xr.open_dataset('./dense_offsets.nc')

    # . . Send offsets for cleaning
    ds_clean = offset_cleaning(ds)

    print('Writing cleaned offsets to disk...')
    ds_clean.to_netcdf('./filt_dense_offsets.nc')

    vel = convert_coords(ds_coords)
