

"""
#FIXME: update doc
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

# function EasyFITS.write(::Type{FitsFile}, 
#     path::AbstractString, 
#     P::PolynLaw{N,T,D}; 
#     kwds...) where {N,T,D}
    
#     coefs = get_coefs(P)
#     hdr = FitsHeader("Dimension" => N, "Degree" => D)
#     EasyFITS.write(FitsFile, path, hdr, coefs; kwds...)
#     nothing
# end

# function EasyFITS.read(::Type{<:PolynLaw{N,T,D}}, 
#     path::AbstractString) where {N,T,D}
    
#     coefs = read(FitsArray, path)
#     hdr = read(FitsHeader, path)

#     return PolynLaw{hdr["DIMENSION"],eltype(coefs),hdr["DEGREE"]}(coefs)
# end

function (P::PolynLaw{N,T,D})(pt) where {N,T,D}
    
    mdl = PolynMdl{N,T,D}(pt)
    
    return mdl(get_coefs(P))
end

# function (P::PolynLaw{N,T,D})(coord::AbstractArray{U,N}) where {N,T,D,U}
    
#     map = Array{T,N}(undef, size(coord))
#     for k in eachindex(coord)
#         map[k] = P(coord[k])
#     end

#     return map
# end

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



"""
    Polyn(coefs)

Structure which contains the `L` coefficients of a `N` dimensional polynomial of degree `D`.

Example of stored coefficients for a 2D polynomials:
            P(x,y) = a_1 + a_2*x + a_3*y + a_4*x^2 + a_5*x*y + a_6*y^2 + ...

            #FIXME: update doc
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

