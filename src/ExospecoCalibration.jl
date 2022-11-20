#
# ExospecoCalibration.jl --
#
# ExospecoCalibration is a module which allows calibrating the IRDIS/LSS data 
# of the SPHERE instrument.
#
#-------------------------------------------------------------------------------
#

#TODO: add export figure in a log file

module ExospecoCalibration


export build_map,
       calibrate_geometry,
       detect_calibrated_coordinates,
       model_dispersion_laws,
       select_region_of_interest


using Base: axes1
using LinearAlgebra
using OptimPackNextGen: BraDi
using PyPlot
const plt = PyPlot
using Statistics
using TwoDimensional


""" Wavelengh units (all wavelengths in nanometers) """
const nm  = 1.0    # one nanometer
const µm  = 1000nm # one micrometer

""" Angular distance units (all angular distances in milliarcseconds) """
const mas = 1.0    # one milliarcsecond

""" Angular distance by pixel """
const rho_pixel = 12.27mas

""" Wavelength of lamps in geometric calibration data (in μm) """
const IRDIS_calibrated_lambda = [0.9877, 1.1237, 1.3094, 1.5451, 1.730, 2.015] .*μm

""" converion constant from degrees to radian """
const deg = π/180


include("types.jl")
include("calibration.jl")
include("tools.jl")

end # module