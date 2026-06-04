# Part 2: Stress Derivation and Crevasse Advection

This directory contains the analysis and figure-generation code for [CITATION].
It does not require running Part 1 &mdash; all necessary input data can be downloaded from the Zenodo archive below.

## Input Data

Processed velocity fields and supporting datasets are archived on Zenodo.

> **Zenodo record:** [DOI: 10.5281/zenodo.20535986]

Download and place the data in `data/` before running the notebooks.
Velocity NetCDF files should go in `data/velocity`.

---

## Notebooks

| Notebook | Figures generated | Key outputs |
|----------|------------------|-------------|
| `01-stress_derivation.ipynb` | Figures 1, 3 | Strain rate tensor, deviatoric stress, von Mises stress, effective stress |
| `02-parcel_tracking.ipynb` | Figures 4&ndash;9 | Lagrangian parcel trajectories, crevasse advection timeseries |
| `03-supplemental.ipynb` | Supplemental figures and movies | S1&ndash;S7, Movies 1&ndash;5 |

Notebooks 01 and 02 can be run independently of each other. Notebook 03 depends on outputs from both and should be run last.

> [!NOTE]
> Movie generation in notebook 03 can be time-intensive;
> individual cells are labeled so you can generate a subset if needed.

## Dependencies

All dependencies are covered by the shared `environment.yml` in the root of this repository. 
See the root-level [README](../README.md) for installation instructions.

When running notebooks using the manual installation, make sure the **ISCE2** kernel is selected.

## `blue_ice_tools.py`

All core analysis functions are in `blue_ice_tools.py`. The main entry points are:

### Data loading

`get_urls(shape_or_points, crs='EPSG:3031')`: 
queries the ITS_LIVE catalog and returns datacube URLs that intersect a given GeoDataFrame or list of points.
 
`load_itslive(shape, urls, dt_delta, start_date, end_date)`: 
downloads and preprocesses ITS_LIVE velocity data cubes, clips them to a region of interest, and resamples to monthly timesteps.

> [!NOTE]
> These are only needed if you are not generating your own velocities or downloading from the Zenodo.
> They are provided for exploratory data analysis if you are looking into a new region.

### Strain rate and stress
 
`calc_strain_stress(ds, rotate, deriv_type, ...)`: 
computes the strain rate tensor and derived stress fields from velocity components `vx` and `vy`. 
Returns the input dataset extended with:
 
| Variable | Description |
|----------|-------------|
| `eps_xx`, `eps_yy`, `eps_xy` | Strain rate tensor components |
| `effective` | Effective strain rate |
| `tau_xx`, `tau_yy`, `tau_xy` | Deviatoric stress tensor |
| `sigma1`, `sigma2` | Principal stresses |
| `von_mises` | Von Mises stress |
 
Three gradient methods are available via `deriv_type`:
 
| Value | Method |
|-------|--------|
| `'sav_gol'` | Savitzky-Golay filter *(default; recommended)* |
| `'stencil'` | 8th-order finite difference stencil |
| `'deriv'` | xarray's built-in `.differentiate()` |
 
Set `rotate=True` to rotate the tensor into along- and across-flow components
following Alley et al. (2018).

### Lagrangian parcel tracking
 
`LagrangianTracking(ds, vx, vy, ...)`: builds interpolators over the velocity field and exposes three methods:

* `parcel_advection(x0, y0, t0)`: 
  advects one or more parcels both forward and backward in time from the given start locations using 4th-order Runge-Kutta integration. 
  Returns a combined `xr.Dataset` with dimensions `(parcel, t)`.
 
* `ensemble_advection(clip_shapes, mask_date, N, ...)`: 
  samples N random starting points from a masked region and runs `parcel_advection` over all of them, returning the median-filtered result. 
  Useful for uncertainty quantification.
 
* `velocity_at(x, y, t)`: 
  returns interpolated velocity components at a given location and time (days since the start of the record).
 
`02-parcel_tracking.ipynb` is the best entry point if you want to adapt the advection approach for your own work.

### Plotting utilities
 
Several helpers are included for consistent figure styling:
 
| Function | Purpose |
|----------|---------|
| `format_map(ax)` | Formats axes for Southern Polar Stereographic maps |
| `inset_colorbar(fig, ax, cs, label)` | Places a compact colorbar inside a map panel |
| `plot_arrows(xs, ys, ax)` | Overlays direction arrows on a parcel trajectory |
| `map_patch` / `plot_patch` | Adds a labeled text box to a map or plot panel |
 
---

## Notes

- Stress calculations follow Glen's flow law with default parameters n = 3 and A = 3.5 × 10$^{-25}$ Pa$^{-3}$ s$^{-1}$.
  These can be adjusted via the `n` and `A` arguments to `calc_strain_stress()`.
- Velocities are expected in m/yr and spatial coordinates in meters (EPSG:3031).
  The Lagrangian tracker converts internally to m/day for integration.
- The `LagrangianTracking` class accepts an optional `extra_ds` dataset (e.g. fracture confidence maps) that is sampled along each trajectory alongside the velocity data.
