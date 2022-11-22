# ExospecoCalibration

`ExospecoCalibration` is a [Julia](https://julialang.org/) package used for calibrating long-slit spectroscopy data (Vigan et al., [2008](https://doi.org/10.1051/0004-6361:200810090)) from the SPHERE/IRDIS instrument (Dohlen et al., [2008](https://doi.org/10.1117/12.789786)).

## Installation

To install the package, write down the following code after hitting the key `]` to activate `Pkg` mode:

```julia
add https://github.com/SJJThe/ExospecoCalibration
```

## Usage

### Generating calibration of spatial and spectral dispersion maps

Given a `LAMP` calibration `d_cal` given by the IRDIS instrument and formed by illuminating the slit by 6 monochromatic sources, the spatial and spectral dispersion maps of the detector can be generated following simple commands. 

Choose the parameters of the polynomial dispersion laws in the form of tuples containing the number of dimensions of the polynomial (`1` for 1D polynomial, `2` for 2D) and the degree of the polynomial. Then call the `calibrate_geometry` method which yields the dispersion maps in the form of a `GeoCalib` structure.

A bad pixels map `bpm` can be given to the method in order to take into account bad pixels in the detection of the calibrated spectral lines in `d_cal`. A `study` keyword can be set to `Val(:log)` which makes the method plots the results of the calibration so the user can validate the dispersion laws.

```julia
using ExospecoCalibration

spatial_pol_carac = (2,1)
spectral_pol_carac = (2,5)
Geo = calibrate_geometry(d_cal, bpm; spatial_law_carac= spatial_pol_carac,
                         spectral_law_carac=spectral_pol_carac, study=Val(:log))
```


### Selecting the region of interest

It is possible to select the region of interest that will be taken into account in further post-processing methods. Given a `GeoCalib` structure, it will create a new structure with masked regions outside the `rho_bnds` and `lambda_bnds` boundaries.

```julia
# Boundaries of region of interest
lambda_bnds = (920.0, 1870.0) .* ExospecoCalibration.nm
rho_bnds = (-15.0, 15.0) .* ExospecoCalibration.rho_pixel

Geo_masked = select_region_of_interest(Geo; rho_bnds=rho_bnds, lambda_bnds=lambda_bnds)
```


### Handling Geometric Calibration data

The package also allows writing and reading fits files based on `GeoCalib` structures.

```julia
# Reading Geometric calibration
path_to_struct = "..."
Geo = readfits(path_to_struct)

# Writing GeoCalib structure
writefits(path_to_struct, Geo; overwrite=true)
```


## Dependencies

`ExospecoCalibration` depends on the following packages:
 - [EasyFITS](https://github.com/emmt/EasyFITS.jl)
 - [OptimPackNextGen](https://github.com/emmt/OptimPackNextGen.jl)
 - [TwoDimensional](https://github.com/emmt/TwoDimensional.jl)

For more details, see the `Project.toml` file.
