


"""
    tolerance(atol, rtol, val) -> max(atol, rtol*abs(val))

yields a nonnegative tolerance given an absolute toterance `atol`, a relative
tolerance `rtol` and a typical value `val`.  An exception is thrown if `atol`
or `rtol` are not nonnegative.

"""
tolerance(atol::Real, rtol::Real, val::Real) = begin
    atol ≥ 0 || error("invalid absolute tolerance (atol=$atol)")
    rtol ≥ 0 || error("invalid relative tolerance (rtol=$rtol)")
    T = float(promote_type(typeof(atol), typeof(rtol), typeof(val)))
    return max(T(atol), T(rtol)*T(abs(val)))
end



########## Detection of calibrated coordinates ##########

"""
    detect_spectral_lines(d_cal, bpm) -> ϕ_λ, L_peak

Yields the projected coordinates of the calibrated spectral lines in 
calibration data `d_cal`, sorted from the highest to the lowest wavelength in 
`L_peak`. The BraDi method is used to find the best projection angle `ϕ_λ` that
will ensure the highest heights of the spectral lines peaks.

The keyword `study` can be given the value `Val(:func)` to return only the 
function `f_obj` that gives the sum of peaks heights, given an angle `p` in 
radian. If the value `Val(:log)` is given, it will plot the resulting profile 
`q_perp` as well as the estimated location of the most significant peaks, whose
coordinates are stored in `L_peak`.

"""
function detect_spectral_lines(d_cal::AbstractMatrix{T},
    bpm::AbstractMatrix{T};
    plot::Bool = false,
    func::Bool = false) where {T<:AbstractFloat}
    
    function f_obj(p)
        wgts, ind_prof, prof = project(d_cal, p; weights=bpm)
        inds = find_peaks(prof; weights=wgts, dist=10, nmax=6)
        return sum(prof[inds])
    end
    if func
        return f_obj
    end
    
    ϕ_λ = BraDi.maximize(p -> f_obj(p*deg), [80.0, 100.0]; atol=0.0, rtol=1e-4)[1]
    
    w_perp, ind_q_perp, q_perp = project(d_cal, ϕ_λ*deg; weights=bpm)
    L_peak = find_peaks(q_perp; weights=w_perp, dist=10, nmax=6)
    
    if plot
        figure()
        plt.plot(ind_q_perp, q_perp)
        plt.scatter(ind_q_perp[L_peak], q_perp[L_peak], color="r", marker="x")
        plt.tight_layout()
    end
    
    return ϕ_λ, sort(L_peak, rev=false)
end



"""
    project(A, θ; weights=undef) -> w, x, y

yields the projection of the 2-dimensional array `A` along a direction of angle
`θ` (in radian counter-clockwise) in a Cartesian frame.

Keyword `weights` may be used to specify a 2-dimensional array of nonnegative
weights of same indices as `A`.

The result is a 3-tuple of vectors `w`, `x` and `y` such that:

    w[k] = sum_{i,j} f(x[k] - i*cos(θ) - j*sin(θ))*W[i,j]
    x[k] = k-th projected index
    y[k] = (sum_{i,j} f(x[k] - i*cos(θ) - j*sin(θ))*A[i,j]*W[i,j])/w[k]

where `f(t) = max(1 - abs(t), 0)` is the linear interpolation B-spline and
`W[i,j] = 1` if `weights = undef` (the default) or `W[i,j] = weights[i,j]`
otherwise.

In other words, `w[k]` and `y[k]` are the sum of weights and weighted average
of values at projected coordinate `x[k]`.

If weights are specified, the sums `sum_{i,j}` are for all indices `(i,j)` such
that `W[i,j] > 0`.

"""
function project(A::AbstractMatrix, θ::Real;
                 weights::Optional{AbstractMatrix}=undef)
    
    if weights === undef
        T = (eltype(A) === Float32 ? Float32 : Float64)
    else
        axes(weights) == axes(A) ||
            throw_dimension_mismatch("weights and data must have the same indices")
        T = (promote_type(eltype(weights), eltype(A)) === Float32 ? Float32 : Float64)
    end
    sinθ, cosθ = sincos(T(θ))
    I, J = axes(A)
    zmin, zmax = projection_limits(cosθ, I, sinθ, J)
    xmin = floor(Int, zmin)
    xmax = floor(Int, zmax) + 1
    x = xmin:xmax
    n = length(x)
    w = zeros(T, n)
    y = zeros(T, n)
    off = 1 - xmin
    if weights === undef
        @inbounds for j in J
            for i in I
                z = cosθ*i + sinθ*j
                r = z - floor(z)
                k = floor(Int, z) + off
                if (1 ≤ k < n)&(0 ≤ r ≤ 1)
                    w0 = one(r) - r
                    w1 = r
                    a_ij = T(A[i,j])
                    y[k]   += w0*a_ij
                    y[k+1] += w1*a_ij
                    w[k]   += w0
                    w[k+1] += w1
                end
            end
        end
    else
        @inbounds for j in J
            for i in I
                w_ij = T(weights[i,j])
                if w_ij > 0
                    z = cosθ*i + sinθ*j
                    r = z - floor(z)
                    k = floor(Int, z) + off
                    if (1 ≤ k < n)&(0 ≤ r ≤ 1)
                        w0 = (one(r) - r)*w_ij
                        w1 = r*w_ij
                        a_ij = T(A[i,j])
                        y[k]   += w0*a_ij
                        y[k+1] += w1*a_ij
                        w[k]   += w0
                        w[k+1] += w1
                    end
                end
            end
        end
    end
    @inbounds for k in 1:n
        if w[k] > 0
            y[k] /= w[k]
        else
            w[k] = 0
            y[k] = 0
        end
    end
    return w, x, y
end

function projection_limits(α::T,
                           I::AbstractRange{<:Integer},
                           β::T,
                           J::AbstractRange{<:Integer}) where {T<:AbstractFloat}
    
    zmin, zmax = typemax(T), typemin(T)
    for i in (first(I), last(I)), j in (first(J), last(J))
        z = α*i + β*j
        zmin = min(zmin, z)
        zmax = max(zmax, z)
    end
    return (zmin, zmax)
end




# Used to specify an optional argument of type T.
const Optional{T} = Union{T,UndefInitializer}

@noinline throw_bad_argument(args...) = throw_bad_argument(string(args...))
@noinline throw_bad_argument(mesg::AbstractString) = throw(ArgumentError(mesg))

@noinline throw_dimension_mismatch(args...) = throw_dimension_mismatch(string(args...))
@noinline throw_dimension_mismatch(mesg::AbstractString) = throw(DimensionMismatch(mesg))


"""
    find_peaks(vect; kwds...) -> inds

yields the indices of the most significant local maxima found in vector `vect`.
The indices are returned in decreasing order of peak heights.

The following keywords can be used to tune the algorithm:

* `weights` specifies a vector of nonnegative weights of same indices as
  `vect`.  Only the values of `vect` where the weights are strictly positive
  are taken into account.

* `dist` (5 by default) specifies the minimal distance (in indices units)
  between two peaks.

* `atol` (0.0 by default) and `rtol` (0.1 by default) specify the absolute and
  relative tolerances for the detection threshold.  If `atol` is `NaN`, the
  detection threshold is `rtol*maximum(vect)`; otherwise, the detection
  threshold is `max(atol, rtol*maximum(vect))`.  All selected peaks have values
  greater or equal the detection threshold.

* `nmax` (no limits by default) specifies the maximum number of peaks to
  detect.

"""
find_peaks(vect::AbstractVector; kwds...) = find_peaks!(copy(vect); kwds...)

"""
    find_peaks!(vect; kwds...) -> inds

yields the indices of the most significant local maxima found in vector `vect`
destroying the contents of `vect` in the process.  The indices are returned in
decreasing order of peak heights.

Call `find_peaks` to avoid overwriting `vect`.

See `find_peaks` for a list of accepted keywords.

"""
function find_peaks!(vect::AbstractVector{<:Real};
                     weights::Optional{AbstractVector}=undef,
                     dist::Integer = 5,
                     rtol::Real = 0.1,
                     atol::Real = 0.0,
                     nmax::Integer = typemax(Int))
    
    dist ≥ 1 || throw_bad_argument("minimal distance must be ≥ 1")
    0 ≤ rtol ≤ 1 || throw_bad_argument("invalid relative threshold")
    dst = Int(dist) - 1
    vmin = typemin(eltype(vect))
    vtol = float(atol)
    I = Base.axes1(vect)
    I_first, I_last = first(I), last(I)
    inds = Int[]
    if weights !== undef
        axes(weights) == axes(vect) ||
            throw_dimension_mismatch("weights and data must have the same indices")
        @inbounds for i in eachindex(vect, weights)
            if !(weights[i] > 0)
                vect[i] = vmin
            end
        end
    end
    while length(inds) < nmax
        vmax, imax = findmax(vect)
        if length(inds) < 1
            # Compute selection threshold.
            v = oftype(vtol, rtol*vmax)
            if isnan(vtol)
                vtol = v
            else
                vtol = max(vtol, v)
            end
        end
        vmax ≥ vtol || break
        push!(inds, imax)
        @inbounds for i in max(I_first, imax - dst):min(I_last, imax + dst)
            vect[i] = vmin
        end
    end
    return inds
end



"""
    LineStat(cog_int, cog_pos, max_int, max_pos)

Structure that contains the level and position of the center of gravity and maximum
of a closed area around a peak of the calibration wavelength data. The area is 
defined in the function `find_projected_positions`.

    #FIXME: update doc
"""
struct LineStat
    cog_int::Float64 # intensity of center of gravity
    cog_pos::Point{Float64} # projected position of center of gravity
    max_int::Float64 # intensity of maximum
    max_pos::Point{Float64} # projected position of maximum
end



"""
    find_spectral_line_paths(d_cal, bpm, L_peak, ϕ_λ) -> C_paths

find the positions of each wavelength rays in the calibration 'd_cal' along its first
index coordinates. The results are stored in 'C_paths', a vector containing for each 
calibrated wavelength, a vector of their projected positions. The projection is 
defined thanks to the angle `ϕ_λ` which must be in radian.

Keywords can be given as to change parameters in the 'find_projected_positions' 
function.
#FIXME: update doc
"""
function find_spectral_line_paths(d_cal::AbstractArray{Float64,2}, 
    bpm::AbstractArray{Float64,2},
    L_peak::AbstractVector{Int},
    ϕ_λ::Float64;
    kwds...)
    
    C_paths = Vector{Point{Float64}}[]
    for l in L_peak
        push!(C_paths, find_projected_positions(d_cal, bpm, l, ϕ_λ; kwds...))
    end
    
    return C_paths
end



"""


#FIXME: update doc

"""
function find_projected_positions(d_cal::AbstractArray{Float64,2}, 
    bpm::AbstractArray{Float64,2},
    l::Int, 
    ϕ_λ::Float64;
    half_width_box::Int = 10,
    ltol::Float64 = 0.1,
    sel_max::Bool = false)

    # select all boxes located around the mean projected position l
    sinphi, cosphi = sincos(ϕ_λ)
    buff = LineStat[] # stat of lines container
    imin, imax = firstindex(axes(d_cal,1)), lastindex(axes(d_cal,1))
    jmin, jmax = firstindex(axes(d_cal,2)), lastindex(axes(d_cal,2))
    lsummax = typemin(Float64)
    for i in imin:imax
        j_center = round(Int, (l - cosphi*i)/sinphi)
        J = max(jmin, j_center-half_width_box):min(jmax, j_center+half_width_box) # boxing of peak
        stat = find_peak_in_box(d_cal, bpm, i, J)
        if stat.cog_int > 0
            push!(buff, stat)
            lsummax = max(lsummax, stat.cog_int)
        end
    end

    # select which boxes contains a peak
    thresh = ltol*lsummax
    if sel_max
        Z = [buff[i].max_pos for i in eachindex(buff) if buff[i].max_int >= thresh]
    else
        Z = [buff[i].cog_pos for i in eachindex(buff) if buff[i].cog_int >= thresh]
    end
    
    return Z
end



"""
    find_peak_in_box(d_cal, bpm, i, J_box) -> LineStat

yields a LineStat structure containing the peak information at the first dimension
index `i` and in the `J_box`.

#FIXME: update doc
"""
function find_peak_in_box(d_cal::AbstractArray{Float64,2},
    bpm::AbstractArray{Float64,2},
    i::Int,
    J_box::AbstractVector{Int})

    lmax, indjmax = typemin(Float64), Float64(0)
    lsum = Float64(0)
    jlsum = Float64(0)
    for j in J_box
        if bpm[i,j] > 0
            val = Float64(d_cal[i,j])
            if val > lmax
                lmax, indjmax = val, j
            end
            if val > 0
                lsum += val
                jlsum += val*j
            end
        end
    end
    if lsum > 0
        return LineStat(lsum, Point(i, jlsum/lsum), lmax, Point(i, indjmax))
    else
        return LineStat(0.0, Point(1.0,1.0), 0.0, Point(1.0,1.0))
    end
end



"""
#FIXME: update doc
"""
function find_edges_coro(ϕ_λ::Float64, 
    d_cal::AbstractArray{Float64,2},
    bpm::AbstractArray{Float64,2},
    L_peak::AbstractVector;
    half_width_box::Int = 10,
    plot::Bool = false,
    study::Bool = false)

    sinphi, cosphi = sincos(ϕ_λ)
    imin, imax = firstindex(axes(d_cal,1)), lastindex(axes(d_cal,1))
    jmin, jmax = firstindex(axes(d_cal,2)), lastindex(axes(d_cal,2))
    edges_coro = Vector{Point{Float64}}[]
    if study
        prof_stud = []
        ind_prof_stud = []
        ind_edges_stud = []
        tol_stud = 0
    end
    for l in 1:length(L_peak)
        # account for only suroundings of ray l
        mask = zeros(size(d_cal))
        for i in imin:imax
            if i > 50 && i < 920 # select field of view
                j_center = round(Int, (L_peak[l] - cosphi*i)/sinphi)
                J = max(jmin, j_center-half_width_box):min(jmax, j_center+half_width_box)
                mask[i,J] .= bpm[i,J]
            end
        end

        # clipping of saturated values
        # done to avoid any remaining invalid pixels compromising the detection of 
        # edges
        med = median(d_cal[mask .== 1.0])
        vari = 1/(length(d_cal[mask .== 1.0]) - 1) * sum((d_cal[mask .== 1.0] .- med).^2)
        #mad = median(abs.(d_cal[mask .== 1.0] .- med))
        lamp_clip = copy(d_cal) .* mask
        lamp_clip[(lamp_clip .- med)./sqrt(vari) .> 4] .= med

        # get angular separation profil of ray (angle must be ϕ_λ + 90)
        wgts, ind_prof, prof = project(lamp_clip, ϕ_λ + 90deg; weights=mask)

        # get ends of edges of segment
        if study
            edges, tol = find_edges(prof; weights=wgts, study=study)
        else
            edges = find_edges(prof; weights=wgts)
        end

        i_l_d = L_peak[l]*cosphi - ind_prof[edges[1]]*sinphi
        j_l_d = L_peak[l]*sinphi + ind_prof[edges[1]]*cosphi
        i_l_u = L_peak[l]*cosphi - ind_prof[edges[2]]*sinphi
        j_l_u = L_peak[l]*sinphi + ind_prof[edges[2]]*cosphi
        
        push!(edges_coro, [Point(i_l_d, j_l_d), Point(i_l_u, j_l_u)])
        if study && (l == 4)
            prof_stud = prof
            ind_prof_stud = ind_prof
            ind_edges_stud = [i_l_d, i_l_u]
            tol_stud = tol
        end
    end

    if plot
        figure()
        plt.imshow(d_cal .* bpm, origin="lower", interpolation="none", 
                   aspect="auto", cmap="gnuplot")
        plt.colorbar()
        for E_l in edges_coro
            plt.scatter(E_l[1][2], E_l[1][1], color="r", marker="x")
            plt.scatter(E_l[2][2], E_l[2][1], color="r", marker="x")
        end
        plt.tight_layout()
    end
    
    if study
        return edges_coro, prof_stud, ind_prof_stud, ind_edges_stud, tol_stud
    else
        return edges_coro
    end
end



"""
    find_edges(vect; atol=0.0, rtol=0.2, weights=undef) -> I

yields the indices of the edges found in vector `vect`.  The method searches
the following motif:

```
    +--------+    +--------+
    |        |    |        |
----+        +----+        +--- ...
   I[1]     I[2] I[3]     I[4]  ...
```

in `vect` and edges are detected by comparing values to a threshold computed
as:

    tol = vmin + max(atol, rtol*(mean(vect) - vmin))

with `vmin = minimum(vect)` and where absolute and relative tolerances can be
specified by keywords `atol` and `rtol`.
#FIXME: update doc

"""
function find_edges(vect::AbstractVector;
                    weights::Optional{AbstractVector}=undef,
                    atol::Real = 0.0,
                    rtol::Real = 0.2,
                    study::Val = Val(false))
    vmin = typemax(eltype(vect))
    vsum = zero(eltype(vect))
    vcnt = 0
    I = axes1(vect)
    if weights === undef
        @inbounds for i in I
            val = vect[i]
            if !isnan(val)
                vmin = ifelse(val < vmin, val, vmin)
                vsum += val
                vcnt += 1
            end
        end
    else
        axes1(weights) == I ||
            throw_dimension_mismatch("weights and data must have the same indices")
        @inbounds @simd for i in I
            val = vect[i]
            if (!isnan(val))&(weights[i] > 0) 
                vmin = ifelse(val < vmin, val, vmin)
                vsum += val
                vcnt += 1
            end
        end
    end
    vcnt > 0 || error("no valid data")
    tol = vmin + tolerance(atol, rtol, vsum/vcnt - vmin)
    #println("vmin=$vmin, vavg=$(vsum/vcnt), tol=$tol")

    # Set initial flag with first valid value, then register the index of each
    # edge.
    inds = Int[]
    flag = false
    unset = true
    if weights === undef
        @inbounds for i in I
            val = vect[i]
            if !isnan(val)
                if unset
                    flag = !(val > tol)
                    unset = false
                elseif (val > tol) == flag
                    push!(inds, i)
                    flag = !flag
                end
            end
        end
    else
        @inbounds for i in I
            val = vect[i]
            if (!isnan(val))&(weights[i] > 0)
                if unset
                    flag = !(val > tol)
                    unset = false
                elseif (val > tol) == flag
                    push!(inds, i)
                    flag = !flag
                end
            end
        end
    end

    if study === Val(true)
        return inds, tol
    else
        return inds
    end
end




########## Fit of dispersion laws ##########

"""
#FIXME: update doc
"""
function fit_spectral_law(C_paths::Vector{Vector{Point{Float64}}},
    calibrated_wavelengths::AbstractVector{Float64},
    polynom_carac::NTuple{2,Int})

    H, d = build_model_fit(polynom_carac, C_paths, calibrated_wavelengths)

    A = H'*H
    b = H'*d

    return ldiv!(cholesky!(Symmetric(A)), b)    
end


"""
#FIXME: update doc
"""
function fit_spatial_law(Edges::Vector{Vector{Point{Float64}}},
    polynom_carac::NTuple{2,Int})
    
    Pts = Vector{Point{Float64}}[]
    for p in 1:2
        pts = Point{Float64}[]
        for l in 1:length(Edges)
            push!(pts, Edges[l][p])
        end
        push!(Pts, pts)
    end
    sum_j = 0.0
    for l in 1:length(Edges)
        p_l_d, p_l_u = Edges[l][1], Edges[l][2]
        sum_j += p_l_u[1] - p_l_d[1]
    end
    delta_rho = sum_j/length(Edges) .* rho_pixel
    data = [-delta_rho/2, delta_rho/2]
    H, d = build_model_fit(polynom_carac, Pts, data)

    A = H'*H
    b = H'*d

    return ldiv!(cholesky!(Symmetric(A)), b)
end


"""
#FIXME: update doc
"""
function build_model_fit(polynom_carac::NTuple{2,Int},
    Pts::Vector{Vector{Point{Float64}}},
    data::AbstractVector{Float64})
    
    d = Float64[]
    H = Matrix{Float64}[]
    for l in 1:length(Pts)
        Pt = Pts[l]
        append!(d, fill(data[l], length(Pt)))
        for p in Pt
            Pmdl = get_mdl(PolynMdl{polynom_carac[1],Float64,polynom_carac[2]}(p))
            v = []
            for i in 1:length(Pmdl)
                v = append!(v, Pmdl[i])
            end
            push!(H, v')
        end
    end

    return vcat(H...), d
end


"""
"""
function build_map(P::PolynLaw{N,T,D},
    s::Tuple{Int,Int}) where {N,T,D}
    
    map = Array{T,N}(undef, s)
    for i in 1:size(map, 1)
        for j in 1:size(map, 2)
            pt = Point(i,j)
            map[i,j] = P(pt)
        end
    end
    
    return map
end

