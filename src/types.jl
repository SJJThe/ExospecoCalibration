

"""
    GeoCalib(rho, lambda, mask)

Yields a `GeoCalib` structure which contains the angular separation map `rho`, the 
wavelength map `lambda` and a mask where invalid pixels are flagged by a 0 value.

The `Base` functions `axes`, `size` and `eltype` can be applied to this structure.

It is possible to save a `GeoCalib` structure `G` as a .fits file using the `EasyFITS`
package:

    EasyFITS.writefits(PATH, G; overwrite=true)

and to load the structure:

    EasyFITS.readfits(GeoCalib, PATH)

"""
struct GeoCalib{T<:AbstractFloat,R<:AbstractMatrix{T},L<:AbstractMatrix{T},
                M<:AbstractMatrix{T}}
    rho::R
    lambda::L
    mask::M

    function GeoCalib(rho::R, 
                      lambda::L, 
                      mask::M = ones(T, size(rho))) where {T<:AbstractFloat,
                                                           R<:AbstractMatrix{T},
                                                           L<:AbstractMatrix{T},
                                                           M<:AbstractMatrix{T}}
        @assert axes(rho) == axes(lambda) == axes(mask)
        return new{T,R,L,M}(rho, lambda, mask)
    end
end

Base.axes(G::GeoCalib) = axes(G.rho)
Base.size(G::GeoCalib) = size(G.rho)
Base.eltype(G::GeoCalib) = eltype(typeof(G))
Base.eltype(::Type{<:GeoCalib{T}}) where {T} = T
Base.show(io::IO, G::GeoCalib{T}) where {T} = begin
    print(io,"GeoCalib{$T}:")
    print(io,"\n - spatial dispersion map `rho` : ",typeof(G.rho))
    print(io,"\n - spectral dispersion map `lambda` : ",typeof(G.lambda))
    print(io,"\n - mask of valid data `mask` : ",typeof(G.mask))
end


# Extend EasyFITS method to provide HDU name and revision number.
EasyFITS.hduname(::Type{<:GeoCalib}) = ("GEOMETRIC-CALIBRATION", 1)

function EasyFITS.writefits(path::AbstractString,
    G::GeoCalib{T};
    kwds...) where {T}

    rho, lambda, mask = G.rho, G.lambda, G.mask
    arr = Array{T}(undef, size(rho)..., 3)
    arr[:,:,1] = rho
    arr[:,:,2] = lambda
    arr[:,:,3] = mask
    
    name, vers = EasyFITS.hduname(typeof(G))
    hdr = FitsHeader("HDUNAME" => name, 
                     "HDUVERS" => vers,
                     "SLICE_1" => "Angular separation map [mas]",
                     "SLICE_2" => "Wavelength map [nm]",
                     "SLICE_3" => "Mask of valid data")
    
    EasyFITS.writefits(path, hdr, arr; kwds...)
    nothing
end

function EasyFITS.readfits(::Type{GeoCalib}, path::AbstractString)
    
    arr = readfits(path)

    return GeoCalib(arr[:,:,1], arr[:,:,2], arr[:,:,3])
end




"""
    PolynLaw{N,T,D}(a)

Yields a structure containing the `N`-dimensional polynomial law of degree `D` defined by
the coefficients in `a` of type `T`.

Write `P(x)` to compute the value of the polynomial `P` at `x` (either a `Real` for 1D 
polynomials or a `Point` for 2D polynomials).

"""
struct PolynLaw{N,T,D} <: Function
    coefs::AbstractVector{T}

    function PolynLaw{N,T,D}(a::AbstractVector{T}) where {N,T,D}
        D >= 0 || error("got a negatif polynomial degree")
        if N == 1
            n = D + 1
        elseif N == 2
            n = ((D + 1)*(D + 2)) >> 1
        else
            error("$(N) dimensional polynomial are not implemented")
        end
        length(a) == n || error("Expected $(n) coefficients, got $(length(a))")
        return new{N,T,D}(a)
    end
end
get_coefs(P::PolynLaw{N,T,D}) where {N,T,D} =P.coefs

function (P::PolynLaw{N,T,D})(pt) where {N,T,D}
    
    mdl = PolynMdl{N,T,D}(pt)
    
    return mdl(get_coefs(P))
end




"""
    PolynMdl{N,T,D}(x)

Yieds a list of the parameter `x` up to the powers of the decomposition of the 
`N`-dimensional polynomial of degree `D`.

This structure is used in the computation of the polynomial law structure `PolynLaw`.

"""
struct PolynMdl{N,T,D} <: Function
    mdl::AbstractVector{AbstractVector{T}}

    function PolynMdl{1,T,D}(x::U) where {T,D,U<:Real}
        D >= 0 || error("got a negatif polynomial degree")
        mdl = AbstractVector[]
        for d in 0:D
            push!(mdl, [x^d])
        end
        return new{1,T,D}(mdl)
    end
    function PolynMdl{2,T,D}(coord::Point{U}) where {T,D,U<:Real}
        D >= 0 || error("got a negatif polynomial degree")
        x, y = coord
        mdl = AbstractVector[]
        for d in 0:D
            vars = []
            degx = d
            degy = 0
            while degx >= 0 && degy <= d
                push!(vars, x^degx*y^degy)
                degx -= 1
                degy += 1
            end
            push!(mdl, vars)
        end
        return new{2,T,D}(mdl)
    end
    function PolynMdl{N,T,D}(coord::Point{U}) where {N,T,D,U<:Real}
        error("$(N) dimensional polynomial are not implemented")
    end
end
get_mdl(P::PolynMdl) = P.mdl

function (P::PolynMdl{N,T,D})(a::AbstractVector{T}) where {N,T,D}
    
    if N == 1
        n = D + 1
    elseif N == 2
        n = ((D + 1)*(D + 2)) >> 1
    end
    length(a) == n || error("Expected $(n) coefficients, got $(length(a))")
    result = T(0)
    mdl = P.mdl
    c_min = 1
    for d in 1:length(mdl)
        coord = mdl[d]
        c_max = c_min + length(coord) - 1
        coefs = a[c_min:c_max]
        result += coord'*coefs
        c_min = c_max + 1
    end

    return result
end

