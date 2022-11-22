
using EasyFITS
using ExospecoCalibration
using PyPlot
const plt = PyPlot


# Load data
path_to_data = "/home/samuel/Documents/Data/IRDIS_LSS/HIP43620"

d_cal = read(FitsArray, joinpath(path_to_data, "preprocessed_lamp.fits"))
bpm = read(FitsArray, joinpath(path_to_data, "bad_pixels_map.fits"))


# Parameters of the dispersion laws (dimensions, degree)
spectral_law_carac = (2,5)
spatial_law_carac = (2,1)

# Generated calibration maps of the detector geometry
Geo = calibrate_geometry(d_cal, bpm; spatial_law_carac=spatial_law_carac,
                         spectral_law_carac=spectral_law_carac,
                         wordy=true, study=Val(:log))


# Boundaries of region of interest
rho_map, lambda_map = get_spatial_map(Geo), get_spectral_map(Geo)
lambda_bnds = (920.0, 1870.0) .* ExospecoCalibration.nm
rho_bnds = (-15.0, 15.0) .* ExospecoCalibration.rho_pixel

d = read(FitsArray, joinpath(path_to_data, "preprocessed_data.fits"))
figure()
plt.imshow(d[:,:,1], origin="lower", interpolation="none", aspect="auto", 
           cmap="gnuplot", vmin=0.0, vmax=1.0)
plt.colorbar()
contour(lambda_map; levels=lambda_bnds, colors="tab:blue")
contour(rho_map; levels=rho_bnds, colors="green")

# Select the region of interest taken into account
Geo_masked = select_region_of_interest(Geo; rho_bnds=rho_bnds,
                                       lambda_bnds=lambda_bnds)


# Save the geometric calibration
write(FitsFile, joinpath(path_to_data, "geometric_calibration.fits"), 
      Geo_masked; overwrite=true)

