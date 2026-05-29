import os
from functools import partial

import numpy as np
import xarray as xr
import rioxarray as rxr
import dask.array as da
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

from shapely.geometry import Polygon
from shapely.geometry import Point

from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from scipy.signal import savgol_filter as sg
from scipy.ndimage import median_filter, convolve1d, binary_erosion
from scipy.interpolate import RegularGridInterpolator

import imageio.v2 as imageio

##########################################################################

def get_urls(*args, crs='EPSG:3031'):
    """
    Retrieve URLs from the ITS_LIVE catalog that intersect with a given GeoDataFrame or set of points.

    Parameters
    ----------
    shape : gpd.GeoDataFrame
        A GeoDataFrame specifying the region of interest in a given CRS.
    points : list, tuple, np.ndarray
        A point (x, y) or a list of points [(x1, y1), (x2, y2), ...].
        Each point must be in the dimensions of the EPSG number passed
    crs : int, optional
        EPSG code used for filtering or projecting the GeoDataFrame. Default is 'EPSG:3031'.

    Returns
    -------
    np.ndarray
        Unique ITS_LIVE data cube URLs intersecting the specified geometry or EPSG code.

    Raises
    ------
    TypeError
        If `epsg` is not an integer, or if inputs are not of expected types.
    ValueError
        If the `shape` does not have a defined CRS, or `points` are not properly formatted.

    Notes
    -----
    - If both `shape` and `points` are None, the function returns URLs filtered by the EPSG code.
    - Points are expected in (x, y) format.

    Examples
    --------
    >>> from shapely.geometry import Polygon
    >>> import geopandas as gpd
    >>> polygon = Polygon([(-60, -60), (-60, -61), (-61, -61), (-61, -60), (-60, -60)])
    >>> gdf = gpd.GeoDataFrame(geometry=[polygon], crs='EPSG:4326')
    >>> urls = get_urls(shape=gdf)
    >>> print(urls)
    """
    # preload ITS_LIVE catalog
    catalog_url = 'https://its-live-data.s3.amazonaws.com/datacubes/catalog_v02.json'
    catalog = gpd.read_file(catalog_url).to_crs(crs)  # Set to CRS of area

    # Empty list to store urls
    urls = []

    # iterate through passed arguments
    for arg in args:
        # If given shape, get URLs
        if isinstance(arg, gpd.GeoDataFrame):
            shape = arg.to_crs(crs)
            catalog_sub = catalog.sjoin(shape, how='inner')
            urls.append(catalog_sub['zarr_url'].drop_duplicates().to_numpy())

        # If given a point, get URL
        elif isinstance(arg, (list, np.ndarray, tuple)):
            geometry = [Point(x, y) for x, y in arg]
            points_gdf = gpd.GeoDataFrame(geometry=geometry, crs=crs)
            catalog_sub = catalog.sjoin(points_gdf, how='inner')
            urls.append(catalog_sub['zarr_url'].drop_duplicates().to_numpy())

        # Anything else raises a type error
        else:
            raise TypeError(f"`get_urls()` received an invalid argument type: {type(arg)}")

    return np.concatenate(urls)  # return URLs

##########################################################################

def load_itslive(
    shape,
    urls: str | list | np.ndarray = None,
    engine: str = 'zarr',
    dt_delta: float = 18,
    start_date: str | np.datetime64 = "2018-07-01",
    end_date: str | np.datetime64 = "2023-01-31",
    chunks: str | dict = 'auto',
    crs: str = 'EPSG:3031'
):
    """
    Download and process ITS_LIVE velocity data cubes, optionally clipping them to a specified geometry.

    Parameters
    ----------
    shape : gpd.GeoDataFrame
        GeoDataFrame defining the region of interest, with a defined CRS.
    urls : list, tuple, or np.ndarray
        URLs of ITS_LIVE data cubes. If not provided, `get_urls` is used with the `shape` and `epsg` parameters.
    dt_delta : int, optional
        Maximum allowed time difference (in days) for image pairs. Default is None (no filtering).
    start_date : string, optional
        First date to include data in the dataset. Default is None (beginning of record).
    end_date : string, optional
        Last date to include in dataset. Default is None (most recent image)
    engine : str, optional
        Engine used for reading Zarr files. Defaults to 'zarr'.
    crs : int, optional
        EPSG code used for filtering or projecting the GeoDataFrame. Default is 'EPSG:3031'.
    
    Returns
    -------
    xr.Dataset
        Concatenated and resampled xarray Dataset representing the velocity data.

    Raises
    ------
    TypeError
        If inputs are not of the expected types.
    ValueError
        If `shape` does not have a defined CRS.

    Notes
    -----
    - Data cubes are resampled to monthly intervals using the mean.
    - Input URLs must correspond to valid ITS_LIVE data cubes.

    Examples
    --------
    >>> from shapely.geometry import Polygon
    >>> import geopandas as gpd
    >>> polygon = Polygon([(-60, -60), (-60, -61), (-61, -61), (-61, -60), (-60, -60)])
    >>> gdf = gpd.GeoDataFrame(geometry=[polygon], crs='EPSG:4326')
    >>> dataset = get_data_cube(shape=gdf)
    >>> print(dataset)
    """
    # If no URLs passed, get from shape
    if urls is None:
        urls = get_urls(shape, crs=crs)

    # preprocess function
    preprocess = partial(
        _preprocess, 
        shape=shape,
        crs=crs,
        dt_delta=dt_delta,
        start_date=start_date,
        end_date=end_date
    )

    # Open with xr.open_mfdataset
    ds = xr.open_mfdataset(
        urls,
        engine=engine,
        preprocess=preprocess,
        chunks=chunks,
        combine='nested',
        concat_dim='mid_date',
        parallel=True,
        decode_timedelta=True  # Added to silence future warning from xarray
    )

    ds = _monthly_resample(ds)  # resample to monthly timesteps
    ds = ds.rio.write_crs(crs)  # Set CRS  
    
    return ds

##########################################################################

def calc_strain_stress(
    ds,
    rotate: bool = True,
    deriv_type: str = 'sav_gol',
    window_length: int = 11,
    polyorder: int = 2,
    dx: int = 50,
    dy: int = 50,
    n: float = 3,
    A: float = 3.5e-25
) -> xr.Dataset:
    """
    Computes strain rates and stresses from velocity fields in an xarray Dataset.

    This function calculates the strain rate tensor and derived stress fields based on
    the velocity components (`vx`, `vy`) in the input dataset. It supports two methods
    for calculating gradients: Savitzky-Golay filtering or finite differences. The user
    can optionally apply tensor rotation to align results with a local coordinate system.

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset containing velocity components `vx` and `vy`.
    sav_gol : bool, optional
        If True, use Savitzky-Golay filtering for gradient computation. 
        If False, use finite differences via xarray's `.differentiate` (default is False).
    dx : int, optional
        Spatial resolution in the x-direction (default is 120). Only relevant if `sav_gol` is True.
    dy : int, optional
        Spatial resolution in the y-direction (default is 120). Only relevant if `sav_gol` is True.
    window_length : int, optional
        Length of the Savitzky-Golay filter window (default is 11). Only relevant if `sav_gol` is True.
    polyorder : int, optional
        Polynomial order for the Savitzky-Golay filter (default is 2). Only relevant if `sav_gol` is True.
    rotate : bool, optional
        If True, rotates the strain rate tensor to align with local flow direction (default is False).
    n : float, optional
        Glen's flow law exponent (default is 3).
    A : float, optional
        Temperature-dependent flow law constant (default is 3.5e-25).

    Returns
    -------
    xr.Dataset
        The input dataset with additional variables:
        - `effective`: Effective strain rate.
        - `eps_xx`: Strain rate tensor component xx.
        - `eps_yy`: Strain rate tensor component yy.
        - `von_mises`: Von Mises stress.
        - `sigma1`: Principal stress component 1.
        - `sigma2`: Principal stress component 2.

    Notes
    -----
    - The Savitzky-Golay filter is applied along the specified axes with constant window
      length and polynomial order. Ensure that `vx` and `vy` are on the appropriate axes.
    - The strain rate tensor is symmetric, with off-diagonal components averaged appropriately.
    - Stress calculations follow Glen's Flow Law and are scaled to kilopascals (kPa).
    - Rotational calculations are based on the arctangent of velocity components to derive
      the local flow direction.

    Examples
    --------
    Using finite differences to compute strain and stress:
    >>> result = compute_strain_stress(ds, dx=100, dy=100, sav_gol=False)

    Using Savitzky-Golay filtering for gradient computation:
    >>> result = compute_strain_stress(ds, sav_gol=True, window_length=13, polyorder=3)

    Rotating strain rates:
    >>> result = compute_strain_stress(ds, rotate=True)
    """
    # Calculate derivative
    if deriv_type == 'sav_gol':
        L = _sav_gol(ds, window_length, polyorder, dx, dy)

    elif deriv_type == 'stencil':
        L = _stencil_grad(ds, dx, dy)

    else:
        L = _deriv(ds)

    E = _calc_effective(L)  # Calculate effective strain rate

    E = _rotate_strain(ds, E) if rotate else E  # rotate strain rates

    S, T = _calc_stress(E, n, A)  # Calculate cauchy stress tensor

    # Turn into datasets:
    E = xr.Dataset(data_vars=E, coords=ds.coords)
    
    # Rename dataset vars for clarity
    E = E.rename_vars({'XX':'eps_xx', 'YY':'eps_yy', 'XY':'eps_xy'})
    S = S.rename_vars({'1':'sigma1', '2':'sigma2', 'VM':'von_mises'})
    T = T.rename_vars({'XX':'tau_xx', 'YY':'tau_yy', 'XY':'tau_xy'})

    return xr.merge([ds, E, S, T], compat='no_conflicts')  # Merge three datasets, return

def _monthly_resample(ds):
    return (ds.sortby('mid_date')
            .resample(mid_date='1ME')
            .mean(dim='mid_date', skipna=True, method='cohorts', engine='flox')
           ).chunk({'x':100, 'y':100, 'mid_date':-1})

def _sav_gol(ds, window_length, polyorder, dx, dy):
    vx = ds.vx
    vy = ds.vy

    return {
        'XX':sg_ufunc(vx, window_length=window_length, polyorder=polyorder, deriv=1, delta=dx, axis=-1),
        'XY':sg_ufunc(vx, window_length=window_length, polyorder=polyorder, deriv=1, delta=dy, axis=-2),
        'YX':sg_ufunc(vy, window_length=window_length, polyorder=polyorder, deriv=1, delta=dx, axis=-1),
        'YY':sg_ufunc(vy, window_length=window_length, polyorder=polyorder, deriv=1, delta=dy, axis=-2)
    }

def _deriv(ds):
    vx = ds.vx
    vy = ds.vy
    
    return {
        'XX':vx.differentiate('x'),
        'XY':vx.differentiate('y'),
        'YX':vy.differentiate('x'),
        'YY':vy.differentiate('y')
    }

def _stencil_grad(ds, dx, dy, mode='nearest'):
    stencil = np.array([1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280])

    vx = ds.vx
    vy = ds.vy

    return {
        'XX':conv_ufunc(vx, stencil, axis=-1, delta=dx),
        'XY':conv_ufunc(vx, stencil, axis=-2, delta=dy),
        'YX':conv_ufunc(vy, stencil, axis=-1, delta=dx),
        'YY':conv_ufunc(vy, stencil, axis=-2, delta=dy)
    }

def _calc_effective(L):
    E = xr.Dataset({
        'XX':L['XX'],
        'XY':0.5 * (L['XY'] + L['YX']),  # Symmetric part for off-diagonal terms
        'YY':L['YY'],
        'effective':da.sqrt(0.5 * ((L['XX'] ** 2) + (L['YY'] ** 2) + (-L['XX'] - L['YY'])**2) + ((0.5 * (L['XY'] + L['YX']))** 2))
    })
    return E

def _rotate_strain(ds, E):
    # Calculate angle of velocity
    theta = da.arctan2(ds.vy, ds.vx)
    
    # Save trig terms as variables for ease
    cos = da.cos(theta)
    sin = da.sin(theta)
    cos2 = cos**2
    sin2 = sin**2
    cos_sin = cos * sin
    
    # rotate strain rates following Alley et. al. 2018
    E_rot = xr.Dataset({
        'XX': (E['XX'] * cos2) +( 2 * E['XY'] * cos_sin) + (E['YY'] * sin2),
        'YY': (E['XX'] * sin2) - (2 * E['XY'] * cos_sin) + (E['YY'] * cos2),
        'XY': ((E['YY'] - E['XX']) * cos_sin) + E['XY'] * (cos2 - sin2),
        'effective':E['effective']
    })
    return E_rot

    # theta = da.arctan2(ds.vy, ds.vx)
    # cos = da.cos(theta)
    # sin = da.sin(theta)

    # # Build strain tensor (..., 2, 2)
    # eps = da.stack([
    #     da.stack([E['11'], E['12']], axis=-1),
    #     da.stack([E['12'], E['22']], axis=-1)
    # ], axis=-2)

    # # Build rotation matrix (..., 2, 2)
    # R = da.stack([
    #     da.stack([ cos,  sin], axis=-1),
    #     da.stack([-sin,  cos], axis=-1)
    # ], axis=-2)

    # # Rotate: R^T eps R
    # eps_rot = da.einsum("...ki,...kl,...lj->...ij", R, eps, R)

    # # Convert back to Dataset
    # return xr.Dataset({
    #     '11': (('mid_date', 'y','x'), eps_rot[...,0,0]),
    #     '22': (('mid_date', 'y','x'), eps_rot[...,1,1]),
    #     '12': (('mid_date', 'y','x'), eps_rot[...,0,1]),
    #     'effective': E['effective']
    # }, coords=ds.coords)

def _calc_stress(E, n, A):
    n = n
    exp = (1 - n) / n
    A = A * (365*24*3600)

    # Deviatoric stress tensor
    T = {
        'XX':(A ** (-1/n)) * ((E['effective'] ** exp) * (E['XX'])),
        'YY':(A ** (-1/n)) * ((E['effective'] ** exp) * (E['YY'])),
        'XY':(A ** (-1/n)) * ((E['effective'] ** exp) * (E['XY']))
    }

    T = xr.Dataset(T)  # Convert Deviatoric stress tensor to dataset

    # Cauchy stress tensor
    S = {
        'XX':(2 * T['XX'] + T['YY']),  # Along flow component
        'YY':(T['XX'] + 2 * T['YY'])   # Across flow component
    }

    # Principal Components
    S['1'] = 0.5*(S['XX'] + S['YY']) + np.sqrt((0.5*(S['XX'] - S['YY']))**2 + T['XY']**2)
    S['2'] = 0.5*(S['XX'] + S['YY']) - np.sqrt((0.5*(S['XX'] - S['YY']))**2 + T['XY']**2)

    S['VM'] = da.sqrt((S['1']**2 + S['2']**2 - S['1'] * S['2']))  # von Mises stress

    S = xr.Dataset(S)  # Convert cauchy stress tensor to dataset
    return S, T

##########################################################################

class LagrangianTracking:
    """Track one or more particles through time and space in 
    an xarray dataset describing a velocity field (vx, vy).
    
    Supports forward and backward propogation from any index."""
    def __init__(self, ds: xr.Dataset, vx='vx', vy='vy', time='mid_date', extra_ds=None, extra_time='mid_date'):
        self.ds = ds
        self.vx_name = vx
        self.vy_name = vy
        self.time_name = time

        self.extra_ds = extra_ds
        self.extra_time = extra_time

        self.x = ds['x'].values
        self.y = ds['y'].values
        self.t = ds[time].values
        self.t_numeric = (self.t - self.t[0]) / np.timedelta64(1, 'D')  # days since start

        # Build Interpolators
        self.vx_interp = RegularGridInterpolator(
            (self.t_numeric, self.y, self.x),
            ds[vx].values, method='nearest'
        )
        
        self.vy_interp = RegularGridInterpolator(
            (self.t_numeric, self.y, self.x),
            ds[vy].values, method='nearest'
        )

    def velocity_at(self, x, y, t):
        """Return velocity components at given location:
        (time [days since start], x, y)"""
        pt = np.array([t, y, x]).T
        u = self.vx_interp(pt)
        v = self.vy_interp(pt)
        
        return np.squeeze(np.stack([u, v], axis=-1))

    def _rk4_integration(self, x, y, t, dt):
        # Get Taylor series slopes, convert from [m/yr] to [m/day]
        k1 = self.velocity_at(x, y, t) / 365.24
        k2 = self.velocity_at(x + k1[0]*dt/2.0, y + k1[1]*dt/2.0, t + dt/2.0) / 365.24
        k3 = self.velocity_at(x + k2[0]*dt/2.0, y + k2[1]*dt/2.0, t + dt/2.0) / 365.24
        k4 = self.velocity_at(x + k3[0]*dt,     y + k3[1]*dt,     t + dt)     / 365.24

        # Combine the pieces, weighting central values higher
        return (k1 + 2*k2 + 2*k3 + k4) / 6.0

    def _advect_parcel(self, x0, y0, t0, direction='forward'):
        """
        Advect a single particle through time.

        Parameters
        ----------
        x0, y0 : float
            Starting coordinates.
        start_index : int
            Index along time dimension to start from.
        direction : str
            'forward' or 'backward'.

        Returns
        -------
        dict
            Trajectory for one particle with keys 'x', 'y', 'time'.
        """
        ## . . Time axis setup
        t = self.t_numeric.copy()  # Time in days since start of record

        t0 = np.datetime64(t0)                          # enforce t0 as numpy datetime
        t0 = (t0 - self.t[0]) / np.timedelta64(1, 'D')  # Get t0 as days since start of array
        t_start = np.argmin(np.abs(t - t0))
        
        # Get only times "ahead" of starting index for direction
        if direction == 'backward':
            t = t[t_start::-1]
        else:
            t = t[t_start:]
        
        dt = np.diff(t)  # Timestep size
        nt = len(t)      # Number of timsteps in this direction

        # . . Initialize empty vectors for coordinate values
        X = np.full(nt, np.nan)
        Y = np.full(nt, np.nan)
        X[0], Y[0] = x0, y0  # Put current location in first entry

        # . . Iterate through time
        for i in range(nt-1):
            ti = t[i]            # Current time [days since start]
            xi, yi = X[i], Y[i]  # Current x, y location
    
            # . . Advect parcel using RK4 integration
            u, v = self._rk4_integration(xi, yi, ti, dt[i])

            # Break loop if a NaN is encountered (left domain)
            if np.isnan(u) or np.isnan(v):
                break

            # . . Move point and save to next timestep
            X[i+1] = xi + u * dt[i]
            Y[i+1] = yi + v * dt[i]
        
        # . . Convert times back to datetimes
        times = self.t[0] + np.array(t[:len(X)], dtype='timedelta64[D]')
        
        # . . Collect dataset variables using .sel(method='nearest')
        ds_sampled = self.ds.sel(
            x=xr.DataArray(X, dims='t'),
            y=xr.DataArray(Y, dims='t'),
            **{self.time_name:xr.DataArray(times, dims='t')},
            method='nearest'
        )
        
        # . . Assign parcel's coordinates explicitly
        ds_sampled = ds_sampled.assign_coords(x=('t', X), y=('t', Y))

        if self.extra_ds is not None:
            extra_sampled = self.extra_ds.sel(
                x=xr.DataArray(X, dims='t'),
                y=xr.DataArray(Y, dims='t'),
                **{self.extra_time: xr.DataArray(times, dims='t')},
                method='backfill'
            )
            ds_sampled = xr.merge([ds_sampled, extra_sampled], compat='override')

        return ds_sampled
    
    def parcel_advection(self, x0, y0, t0):
        """
        Advect one or multiple points, each with its own start index.

        Parameters
        ----------
        x0, y0 : float or array-like
            Starting coordinates (same length as t0).
        t0 : array-like | np.datetime64 | str
            Time(s) to start from for each particle.

        Returns
        -------
        dict
            Combined trajectories across all particles
            The 'x' and 'y' entries have the shape (num parcels, num timesteps)
        """
        # # Detect if scalar input
        # scalar_input = np.isscalar(x0) or (np.ndim(x0) == 0)

        # Ensure input is 2D
        x0 = np.atleast_1d(x0)
        y0 = np.atleast_1d(y0)
        t0 = np.atleast_1d(t0)
        n_parcels = len(x0)  # Number of parcels

        # . . Iterate through parcels
        all_trajs = []
        for i in range(n_parcels):
            # Get forward and backward trajectories
            fwd = self._advect_parcel(x0[i], y0[i], t0[i], direction='forward')
            bwd = self._advect_parcel(x0[i], y0[i], t0[i], direction='backward')

            # Combine trajectories
            combined = xr.concat([bwd.isel(t=slice(1, None)).isel(t=slice(None, None, -1)), fwd], dim='t')
            combined = combined.assign_coords(
                parcel=i,
                # x0=x0[i],
                # y0=y0[i],
                # t0=t0[i]
            )  # Create parcel dimension for merge
            all_trajs.append(combined)

        # . . Concat on parcel dimension
        ds_out = xr.concat(all_trajs, dim='parcel', coords='different', compat='equals')
        return ds_out

    def ensemble_advection(
        self,
        clip_shapes,
        mask_date,
        mask_var='fracture_conf',
        mask_time_dim='mid_date',
        N=800,
        erode_iters=20,
        seed=7,
        filt_size=(1, 10),
        t0=None,
    ):
        """
        Sample N random parcels from a masked region and run ensemble parcel advection.
 
        Clips a copy of the dataset supplied at initialisation (``extra_ds`` if
        provided, otherwise ``ds``) to define a valid sampling region, erodes the
        boundary inward to avoid edge effects, randomly picks N start points, runs
        ``parcel_advection``, and returns the median-filtered result.
 
        Parameters
        ----------
        clip_shapes : list of dict
            Ordered sequence of clipping operations. Each dict
            must contain:
            - 'gdf' (gpd.GeoDataFrame): geometry used for clipping.
            - 'invert' (bool): if True, pixels *outside* the geometry are kept.
            - 'all_touched' (bool, optional): passed to rioxarray.clip. Default False.
 
            Example::
 
                [
                    {'gdf': gl_poly,       'invert': True,  'all_touched': True},
                    {'gdf': shirase_shape, 'invert': False},
                ]
 
        mask_date : str
            Date string for extracting a single time slice for masking
            (e.g., '2022-03-25'). Pixels where `mask_var` is null at this date
            are excluded from sampling.
        mask_var : str
            Variable tested with isnull() to identify valid pixels.
            Default: 'fracture_conf'.
        mask_time_dim : str
            Name of the time dimension in the dataset used for masking.
            Default: 'mid_date'.
        N : int
            Number of parcels to sample. Default: 800.
        erode_iters : int
            Binary erosion iterations applied to the valid mask before sampling,
            to avoid selecting pixels at region edges. Default: 20.
        seed : int
            Random seed for reproducible parcel selection. Default: 7.
        filt_size : tuple
            size argument forwarded to apply_med_filt for smoothing raw
            trajectories along the time axis. Default: (1, 10).
        t0 : str or None
            Start time for all parcels. Defaults to `mask_date` if not provided.
 
        Returns
        -------
        xr.Dataset
            Median-filtered parcel dataset with dimensions (parcel, t), as
            returned by parcel_advection.
 
        Raises
        ------
        ValueError
            If no valid pixels remain after masking and erosion.
 
        Examples
        --------
        >>> clip_shapes = [
        ...     {'gdf': gl_poly,       'invert': True,  'all_touched': True},
        ...     {'gdf': shirase_shape, 'invert': False},
        ... ]
        >>> parcels = tracker.ensemble_advection(
        ...     clip_shapes=clip_shapes,
        ...     mask_date='2022-03-25',
        ...     N=800,
        ... )
        """
        # . . Copy the initialisation dataset to avoid mutating instance state.
        # Use extra_ds when available (e.g. fracture data), otherwise fall back to ds.
        source = self.extra_ds if self.extra_ds is not None else self.ds
        shelf = source.copy()
 
        # . . Build valid region by sequential clipping
        for clip in clip_shapes:
            gdf = clip['gdf']
            shelf = shelf.rio.clip(
                gdf.geometry, gdf.crs,
                all_touched=clip.get('all_touched', False),
                invert=clip.get('invert', False),
            )
 
        # . . Extract validity mask at mask_date: True where data is absent
        mask = (
            shelf[mask_var]
            .sel({mask_time_dim: mask_date}, method='nearest')
            .isnull()
            .values
        )
 
        # . . Erode boundary inward to avoid edge selection
        eroded = binary_erosion(~mask, iterations=erode_iters)
 
        # . . Sample N random (x, y) points within the valid region
        valid = np.argwhere(eroded)
        if len(valid) == 0:
            raise ValueError(
                "No valid pixels remain after masking and erosion. "
                "Try reducing erode_iters, checking clip geometries, or "
                "verifying that mask_var is populated at mask_date."
            )
        rng = np.random.default_rng(seed)
        sample = valid[rng.choice(len(valid), size=N, replace=False)]
        ys = shelf.y.values[sample[:, 0]]
        xs = shelf.x.values[sample[:, 1]]
        t0_arr = [t0 or mask_date] * N
 
        # . . Advect all parcels and apply median filter along time axis
        parcels_raw = self.parcel_advection(xs, ys, t0_arr)
        return apply_med_filt(parcels_raw, size=filt_size)


##########################################################################

def _preprocess(
    ds: xr.Dataset,
    shape: gpd.GeoDataFrame = None, 
    crs: str = 'EPSG:3031',
    dt_delta: int = None,
    start_date = None,
    end_date = None
) -> xr.Dataset:
    """
    Preprocess ITS_LIVE velocity data cubes before concatenation.

    Steps:
    - Filter dataset to include specific date range if specified.
    - Filter dataset for time delta (if `dt_delta` is set).
    - Select relevant velocity variables.
    - Clip data to a GeoDataFrame region (if provided).
    - Sort by 'mid_date'.

    Parameters
    ----------
    ds : xr.Dataset
        Input xarray dataset to preprocess.
    shape : gpd.GeoDataFrame, optional
        GeoDataFrame for spatial clipping.
    crs : int, optional
        EPSG code for the CRS. Defaults to 'EPSG:3031'.
    dt_delta : int, optional
        Maximum allowed time difference (in days) for image pairs. If None, no filtering is applied.
    start_date : string, optional
        First date to include data in the dataset. Default is None (beginning of record).
    end_date : string, optional
        Last date to include in dataset. Default is None (most recent image)

    Returns
    -------
    xr.Dataset
        Preprocessed xarray Dataset.
    """
    ## TODO: 
    ## figure out single satellite selection

    # Filter by start and end date if specified
    if start_date is not None:
        start_date = np.datetime64(start_date)
        ds = ds.where(ds.mid_date >= start_date, drop=True)
    
    if end_date is not None:
        end_date = np.datetime64(end_date)
        ds = ds.where(ds.mid_date <= end_date, drop=True)

    # Filter by time delta if dt_delta is specified
    if dt_delta is not None:
        time_threshold_ns = np.timedelta64(dt_delta, 'D')
        ds = ds.where(ds.date_dt <= time_threshold_ns)

    # Select relevant velocity variables
    ds = ds[['vx', 'vy']]
    
    # Clip to GeoDataFrame if provided
    if shape is not None:
        ds = ds.rio.write_crs(crs)
        ds = ds.rio.clip(shape.geometry, shape.crs)
        
    return ds

##########################################################################

def sg_ufunc(arr, **kwargs):
    """
    Applies Savitzky-Golay filter via xarray's apply_ufunc.
    """
    return xr.apply_ufunc(
        sg, 
        arr,
        kwargs=kwargs,
        dask='parallelized',
        output_dtypes=['float32']
    )

##########################################################################

def conv_ufunc(arr, stencil, axis=-1, mode='nearest', delta=50):
    """Computes the derivative of a function using an eighth-order stencil"""    
    return xr.apply_ufunc(
        convolve1d, arr, stencil,
        # input_core_dims=[['y', 'x']],
        kwargs={'axis':axis, 'mode':mode},
        dask='allowed', output_dtypes=[float]
    ) / delta

##########################################################################

def apply_med_filt(arr, size=3):
    return xr.apply_ufunc(
        median_filter,
        arr,
        kwargs={'size':size}
    )

############ Plotting Utilities  ##############################################################

def map_patch(ax, txt, va='center', fs=36, alpha=0.9):
    ax.text(
        0.03, 0.95, txt,
        transform=ax.transAxes,
        va=va,
        fontsize=fs,
        bbox=dict(facecolor='white', alpha=alpha, edgecolor='none', pad=3)
    )

def plot_patch(ax, txt, va='center', fs=36, alpha=0.9, x=0.01, y=0.93):
    ax.text(
        x, y, txt,
        transform=ax.transAxes,
        va=va,
        fontsize=fs,
        bbox=dict(facecolor='white', alpha=alpha, edgecolor='none', pad=3),
        zorder=1000
    )

def meters_to_km(x, pos):
    return f"{x/1000:.0f}"

def format_map(ax, fs=14, xrot=0, yrot=90):
    ax.xaxis.set_major_formatter(FuncFormatter(meters_to_km))
    ax.yaxis.set_major_formatter(FuncFormatter(meters_to_km))
    ax.tick_params(axis='x', rotation=xrot)
    ax.tick_params(axis='y', rotation=yrot)
    ax.set_xlabel('Southern Polar Stereographic X [km]', fontsize=fs)
    ax.set_ylabel('Southern Polar Stereographic Y [km]', fontsize=fs)
    ax.set_title(None)
    ax.set_aspect('equal')
    ax.grid(ls='--', alpha=0.5)

def inset_colorbar(fig, ax, cs, label, inset=0.035, height_scale=0.4, pad=0.01, cbar_width=0.005):
    pos = ax.get_position()

    cax = fig.add_axes([
        pos.x1 - cbar_width - inset,
        pos.y0 + inset,
        cbar_width,
        pos.height * height_scale
    ])
    cax.set_zorder(10)

    # Make the colorbar
    cbar = fig.colorbar(cs, cax=cax, label=label)

    # Force a draw so that text bounding boxes exist
    fig.canvas.draw()

    # . . Tight bounding box of the whole colorbar (including label!)
    bbox = cbar.ax.get_tightbbox(fig.canvas.get_renderer())

    # Convert from display coords to figure fraction
    inv = fig.transFigure.inverted()
    bbox_fig = inv.transform(bbox)

    # Expand the background by scalar mults of pad
    bg_rect = [
        bbox_fig[0][0] - pad/2,
        bbox_fig[0][1] - 1.5*pad,
        (bbox_fig[1][0] - bbox_fig[0][0]) + pad,
        (bbox_fig[1][1] - bbox_fig[0][1]) + 3*pad
    ]

    # . . Calculate shift needed to align right edge w/ inset number
    bg_right_edge = bg_rect[0] + bg_rect[2]  # left + width
    desired_right_edge = pos.x1 - inset
    xshift = desired_right_edge - bg_right_edge

    bg_bottom_edge = bg_rect[1]
    desired_bottom_edge = pos.y0 + inset
    yshift = desired_bottom_edge - bg_bottom_edge

    cax_pos = cax.get_position()
    cax.set_position([cax_pos.x0+xshift, cax_pos.y0+yshift, cax_pos.width, cax_pos.height])

    bg_rect[0] += xshift  # shift background left position
    bg_rect[1] += yshift  # shift background up
    
    # Background axes
    bg = fig.add_axes(bg_rect, zorder=9)
    bg.set_facecolor((1,1,1,0.7))
    bg.set_xticks([])
    bg.set_yticks([])

    return cbar

def plot_arrows(xs, ys, ax):
    '''
    Plotting function that plots the direction arrows for the parcel location
    produced by bit.parcel_strain_stress().
    '''
    xs = np.array(xs)
    ys = np.array(ys)
    
    # Sample points for arrows (e.g., every 10th point)
    arrow_indices = np.arange(2, len(xs)-1, 5)  # Avoid including the last index
    x_arrows = xs[arrow_indices]
    y_arrows = ys[arrow_indices]
    # Compute direction vectors for arrows (using np.diff to calculate directional vectors)
    dx = np.diff(xs)  # Differences in x
    dy = np.diff(ys)  # Differences in y
    directions = np.sqrt(dx**2 + dy**2)  # Magnitudes of direction vectors
    
    # Normalize the direction vectors
    dx = dx / directions
    dy = dy / directions
    
    # Align direction vectors with sampled arrow positions
    dx_arrows = dx[arrow_indices]  # Aligning with arrow positions
    dy_arrows = dy[arrow_indices]  # Aligning with arrow positions
    
    # Plot the line
    ax.plot(xs, ys, color='white', label='Parcel path')
    
    # Add arrows using quiver
    ax.quiver(x_arrows, y_arrows, dx_arrows, dy_arrows, 
               angles='uv', scale_units='xy', width=.05, color='white')