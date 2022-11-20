

"""
    calibrate_geometry(d_cal [, bpm]) -> rho_map, lambda_map

Yields the geometric calibration in the form of an angular separation map, 
`rho_map` and a wavelength map, `lambda_map`. These two maps gives the spatial and 
spectral coordinates of each pixel of the detector.

Takes as input the pre-processed calibration data `d_cal` produced by mono-
chromatic sources illuminating the slit. A bad pixels map `bpm` can also be 
given as the form of an array to mark the pixels that do not follow a linear 
behavior. By default, `bpm` is an `AbstractArray` filled with ones.

"""
function calibrate_geometry(d_cal::AbstractArray{T,N},
    bpm::AbstractArray{T,N} = ones(size(data,1), size(data,2));
    spectral_law_carac::NTuple{2,Int} = (2, 5),
    spatial_law_carac::NTuple{2,Int} = (2, 1),
    study::Val = Val(false),
    wordy::Bool = false) where {T,N}

    @assert axes(d_cal) == axes(bpm)
    
    wordy && println("+ Geometric calibration")
    # Detect the calibrated coordinates
    C_paths, Edges = detect_calibrated_coordinates(d_cal, bpm; study=study, 
                                                   wordy=wordy)
    # Model the dispersion laws
    spatial_law, spectral_law = model_dispersion_laws(C_paths, Edges; 
                                        spectral_law_carac=spectral_law_carac,
                                        spatial_law_carac=spatial_law_carac, 
                                        wordy=wordy)
                                        
    wordy && println("|-- Building spatial and spectral maps")
    rho_map = build_map(spatial_law, size(d_cal))
    lambda_map = build_map(spectral_law, size(d_cal))

    if study === Val(:log)
        #FIXME: add export figure in a log file
        figure()
        plt.imshow(d_cal .* bpm, origin="lower", interpolation="none", 
                   aspect="auto", cmap="gnuplot")
        plt.colorbar()
        contour(lambda_map; levels=IRDIS_calibrated_lambda, colors="tab:blue")
        contour(rho_map; levels=[0.0], colors="green")        
        plt.tight_layout()
    end

    if study == Val(true)
        return rho_map, lambda_map, C_paths, Edges, spatial_law, spectral_law
    else
        return rho_map, lambda_map
    end
end



"""
"""
function detect_calibrated_coordinates(d_cal::AbstractArray{T,N},
    bpm::AbstractArray{T,N};
    study::Val = Val(false),
    wordy::Bool = false) where {T,N}

    wordy && println("|-- Detection of calibrated coordinates")
    # projected positions of spectral lines
    ϕ_lambda, L_peak = detect_spectral_lines(d_cal, bpm; study=study)
    # spectral line paths
    C_paths = find_spectral_line_paths(d_cal, bpm, L_peak, ϕ_lambda*deg)
    # edges of cronagraphic mask
    Edges = find_edges_coro(ϕ_lambda*deg, d_cal, bpm, L_peak; study=study)
    
    return C_paths, Edges
end



"""
"""
function model_dispersion_laws(C_paths::Vector{Vector{Point{Float64}}},
    Edges::Vector{Vector{Point{Float64}}};
    calibrated_wavelengths::AbstractVector{T} = IRDIS_calibrated_lambda,
    spectral_law_carac::NTuple{2,Int} = (2, 5),
    spatial_law_carac::NTuple{2,Int} = (2, 1),
    wordy::Bool = false) where {T}

    wordy && println("|-- Modeling dispersion laws")
    # fit spectral law parameters
    spectral_prms = fit_spectral_law(C_paths, calibrated_wavelengths, 
                                     spectral_law_carac)
    # fit spatial law parameters
    spatial_prms = fit_spatial_law(Edges, spatial_law_carac)
    
    # build laws 
    spatial_law = PolynLaw{spatial_law_carac[1],T,spatial_law_carac[2]}(spatial_prms)
    spectral_law = PolynLaw{spectral_law_carac[1],T,spectral_law_carac[2]}(spectral_prms)
    
    return spatial_law, spectral_law
end



"""
"""
lambda_ref(lambda_map::AbstractMatrix{T}) where {T,N} = maximum(lambda_map)



"""
"""
function select_region_of_interest(rho_map::AbstractMatrix{T},
    lambda_map::AbstractMatrix{T},
    bpm::AbstractMatrix{T} = ones(size(rho_map));
    lambda_bnds::NTuple{2,T} = (minimum(lambda_map), maximum(lambda_map)),
    rho_bnds::NTuple{2,T} = (-25rho_pixel, 25rho_pixel),
    wordy::Bool = false) where {T}
    
    wordy && println("|-- Select region of interest")

    # select all pixels between lambda_bnds[1] and lambda_bnds[2] as valid
    mask_lambda = lambda_bnds[1] .<= lambda_map .<= lambda_bnds[2]

    # select all pixels between the maxima angular distances and rho_bnds[1]
    # and rho_bnds[2]
    map_gamma_rho = (lambda_ref(lambda_map) ./ lambda_map) .* rho_map
    mask_rho = (minimum(rho_map) .<= map_gamma_rho .<= rho_bnds[1]) + 
               (rho_bnds[2] .<= map_gamma_rho .<= maximum(rho_map))

    # mask of invaled pixels
    mask = bpm .* mask_lambda .* mask_rho

    return mask
end
