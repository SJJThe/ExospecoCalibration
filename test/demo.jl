
using EasyFITS
using ExospecoCalibration
using PyPlot
const plt = PyPlot


# Load data
path_to_data = joinpath(@__DIR__, "..", "data")

d_cal = read(FitsArray, joinpath(path_to_data, "preprocessed_lamp.fits.gz"))
bpm = read(FitsArray, joinpath(path_to_data, "bad_pixels_map.fits.gz"))


# Parameters of the dispersion laws (dimensions, degree)
spectral_law_carac = (2,5)
spatial_law_carac = (2,1)

# Generated calibration maps of the detector geometry
Geo = calibrate_geometry(d_cal, bpm; spatial_law_carac=spatial_law_carac,
                         spectral_law_carac=spectral_law_carac,
                         wordy=true, plot=true)


# Boundaries of region of interest
rho_map, lambda_map = Geo.rho, Geo.lambda
lambda_bnds = (920.0, 1870.0) .* ExospecoCalibration.nm
rho_bnds = (-15.0, 15.0) .* ExospecoCalibration.rho_pixel

# Select the region of interest taken into account
Geo_masked = select_region_of_interest(Geo; rho_bnds=rho_bnds,
                                       lambda_bnds=lambda_bnds)


# Save calibration
writefits(joinpath(path_to_data, "geometric_calibration.fits"), 
          Geo_masked; overwrite=true)


