

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
    spatial_law_carac::NTuple{2,Int} = (2, 1),
    spectral_law_carac::NTuple{2,Int} = (2, 5),
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

    # Create geometric calibration structure
    G = GeoCalib(rho_map, lambda_map, bpm)

    if study === Val(:log)
        figure()
        plt.imshow(d_cal .* bpm, origin="lower", interpolation="none", 
                   aspect="auto", cmap="gnuplot")
        plt.colorbar()
        contour(lambda_map; levels=IRDIS_calibrated_lambda, colors="tab:blue")
        contour(rho_map; levels=[0.0], colors="green")        
        plt.tight_layout()
    end

    if study === Val(:laws)
        return rho_map, lambda_map, C_paths, Edges, spatial_law, spectral_law
    else
        return G
    end
end



"""
    detect_calibrated_coordinates(d_cal, bpm) -> C_paths, Edges

yields the coordinates of each peaks along the path of the spectral bands
present in `d_cal` and stored in `C_paths` as `Point` structures (see
`TwoDimensional` package). The edges of the coronagraphic mask are also detected
and stored as `Point`structures in `Edges`. All of these coordinates are
detected while taking into account the potential bad-pixels in the data, thanks
to a map `bpm` where they are flagged as `0` values.

It is possible to use the keywords `study` set to `Val(:log)` for the user to
check the detection of these different coordinates on the data `d_cal`.

"""
function detect_calibrated_coordinates(d_cal::AbstractArray{T,N},
    bpm::AbstractArray{T,N};
    study::Val = Val(false),
    wordy::Bool = false) where {T,N}

    wordy && println("|-- Detection of calibrated coordinates")
    # projected positions of spectral lines
    ϕ_lambda, L_peak = detect_spectral_lines(d_cal, bpm; study=study)
    # spectral line paths
    C_paths = find_spectral_line_paths(d_cal, bpm, L_peak, ϕ_lambda)
    # edges of cronagraphic mask
    Edges = find_edges_coro(ϕ_lambda, d_cal, bpm, L_peak; study=study)
    
    return C_paths, Edges
end



"""
    model_dispersion_laws(C_paths, Edges) -> spatial_law, spectral_law

yields the spatial and spectral dispersion laws in the form of `PolynLaw`
structures. These laws are fitted on the coordinates contained in `C_paths` and
`Edges` obtained with the function `detect_calibrated_coordinates`.

To characterize the spatial and spectral laws, the user can specify the
dimension (`1` for 1D polynomial, `2` for 2D) and degree of the polynomials as a
`Tuple` via the `spatial_law_carac` and `spectral_law_carac` keywords. By
default, the dispersion laws have the following characteristics:

```
spatial_pol_carac = (2,1)
spectral_pol_carac = (2,5)

```

It is possible to give the wavelengths of the laser sources used in the
calibration process via the keyword `calibrated_wavelengths`. By default the
values are given by the `const` `IRDIS_calibrated_lambda`.

"""
function model_dispersion_laws(C_paths::Vector{Vector{Point{Float64}}},
    Edges::Vector{Vector{Point{Float64}}};
    spatial_law_carac::NTuple{2,Int} = (2, 1),
    spectral_law_carac::NTuple{2,Int} = (2, 5),
    calibrated_wavelengths::AbstractVector{T} = IRDIS_calibrated_lambda,
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
    lambda_ref(lambda_map) -> ref

yields the reference wavelength used by the model EXOSPECO. By default, this
value is the maximum of the `lambda_map`.

"""
lambda_ref(lambda_map::AbstractMatrix{T}) where {T} = maximum(lambda_map)



"""
    select_region_of_interest(rho_map, lambda_map [, bpm]) -> geocalib

Taking a spatial (`rho_map`) and spectral (`lambda_map`) coordinates maps, this
function is used to identify the region that contains the data modeled by the
EXOSPECO direct model. A bad-pixels map can be given as an argument of the
function, but is not necessary (although recommended).

Spatial boundaries corresponding to the surroundings of the coronagraphic mask
can be given in mas using the keyword `rho_bnds` (by default,
`rho_bnds=(-25rho_pixel,25rho_pixel)`, with `rho_pixel` the mas per pixel
value). In the same way, it is possible to define wavelength boundaries via the
keyword `lambda_bnds`.

To define the region of interest, the function needs a reference wavelength. The
user can specify a function to extract this value from a wavelength coordinates
map via the keyword `lambda_ref_func` (`= lambda_ref` by default).

"""
function select_region_of_interest(rho_map::AbstractMatrix{T},
    lambda_map::AbstractMatrix{T},
    bpm::AbstractMatrix{T} = ones(T, size(rho_map));
    rho_bnds::NTuple{2,T} = (-25rho_pixel, 25rho_pixel),
    lambda_bnds::NTuple{2,T} = (minimum(lambda_map), maximum(lambda_map)),
    lambda_ref_func::Function = lambda_ref,
    wordy::Bool = false) where {T}
    
    wordy && println("|-- Select region of interest")

    # select all pixels between lambda_bnds[1] and lambda_bnds[2] as valid
    mask_lambda = lambda_bnds[1] .<= lambda_map .<= lambda_bnds[2]

    # select all pixels between the maxima angular distances and rho_bnds[1]
    # and rho_bnds[2]
    gamma = lambda_ref_func(lambda_map) ./ lambda_map
    map_gamma_rho = gamma .* rho_map
    mask_rho = (minimum(rho_map) .<= map_gamma_rho .<= rho_bnds[1]) + 
               (rho_bnds[2] .<= map_gamma_rho .<= maximum(rho_map))

    # mask of invaled pixels
    mask = bpm .* mask_lambda .* mask_rho

    return GeoCalib(rho_map, lambda_map, mask)
end

function select_region_of_interest(G::GeoCalib{T};
    rho_bnds::NTuple{2,T} = (0.0, 0.0),
    lambda_bnds::NTuple{2,T} = (minimum(G.lambda),maximum(G.lambda)),
    kwds...) where {T}

    return select_region_of_interest(G.rho, G.lambda, G.mask; 
                                     rho_bnds=rho_bnds,
                                     lambda_bnds=lambda_bnds, kwds...)
end
