
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
ρ_map, λ_map = calibrate_geometry(d_cal, bpm;
                                  spectral_law_carac=spectral_law_carac,
                                  spatial_law_carac=spatial_law_carac,
                                  wordy=true, study=Val(:log))
                                  

# Boundaries of region of interest
λ_bnds = (920.0, 1870.0) .* ExospecoCalibration.nm
ρ_bnds = (-15.0, 15.0) .* ExospecoCalibration.rho_pixel

d = read(FitsArray, joinpath(path_to_data, "preprocessed_data.fits"))
figure()
plt.imshow(d[:,:,1], origin="lower", interpolation="none", aspect="auto", 
           cmap="gnuplot", vmin=0.0, vmax=1.0)
plt.colorbar()
contour(λ_map; levels=λ_bnds, colors="tab:blue")
contour(ρ_map; levels=ρ_bnds, colors="green")

# Select the region of interest taken into account
mask = select_region_of_interest(ρ_map, λ_map, bpm; 
                                 lambda_bnds=λ_bnds, rho_bnds=ρ_bnds)


# Save the geometric calibration
geo_calib = Array{Float64,3}(undef, size(d_cal,1), size(d_cal,2), 3)
geo_calib[:,:,1] = ρ_map
geo_calib[:,:,2] = λ_map
geo_calib[:,:,3] = mask
write(FitsFile, joinpath(path_to_data, "geometric_calibration.fits"), geo_calib; 
      overwrite=true)
