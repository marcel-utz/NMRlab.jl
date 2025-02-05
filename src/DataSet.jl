## (c)2025 Marcel Utz

# abstract type CoordMap{T} <: AbstractVector{T} end

@doc raw"""
    struct SpectData{T,N} <: AbstractArray{T,N}

Basic  structure for spectral data.
"""
struct SpectData{T,N} <: AbstractArray{T,N}
    dat::AbstractArray{T,N}
    coord::NTuple{N,AbstractVector}
end

import Base.size
import Base.getindex
import Base.setindex!
import Base.IndexStyle
import Base.showarg

size(S::SpectData) = size(S.dat)
getindex(S::SpectData, k::Integer) = getindex(S.dat,k)
setindex!(S::SpectData, v, k::Integer) = setindex!(S.dat,v,k)
IndexStyle(S::SpectData) = IndexStyle(S.dat)

Base.showarg(io::IO, A::SpectData, toplevel) = print(io, typeof(A), " with coords:", A.coord)

## support broadcasting ------------------------------------------------------
Base.BroadcastStyle(::Type{<:SpectData}) = Broadcast.ArrayStyle{SpectData}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SpectData}}, ::Type{ElType}) where ElType
    # Scan the inputs for the ArrayAndChar:
    A = find_spdta(bc)
    # Use the char field of A to create the output
    SpectData(similar(Array{ElType}, axes(bc)), A.coord)
end

"`A = find_aac(As)` returns the first ArrayAndChar among the arguments."
find_spdta(bc::Base.Broadcast.Broadcasted) = find_spdta(bc.args)
find_spdta(args::Tuple) = find_spdta(find_spdta(args[1]), Base.tail(args))
find_spdta(x) = x
find_spdta(::Tuple{}) = nothing
find_spdta(a::SpectData, rest) = a
find_spdta(::Any, rest) = find_spdta(rest)

## ----------------------------------------------------------------------------

## support slicing ------------------------------------------------------------
getindex(S::SpectData{T,N},I::Vararg{<:Integer,N}) where {T,N} = getindex(S.dat,I...)

function getindex(S::SpectData,I...)
    adat=S.dat[I...]
    if all(x->x isa Integer,I)
        return adat
    else
        newcoords=[]
        for (k,i) in enumerate(I)
            if !(i isa Integer) # integer indices are dropped
                push!(newcoords,S.coord[k][i])
            end
        end
    
    return SpectData(adat,(newcoords...,))
    end
end 


@doc raw"""
    function coords(S::SpectData)
    function coords(S::SpectData,k::Integer)

returns a tuple with the coordinates of `S`, analogous to `axes(S)`.
The second form returns the coordinate of the `k`-th dimension.
"""
coords(S::SpectData) = S.coord
coords(S::SpectData,k::Integer) = S.coord[k]

function SpectData(A::AbstractArray)
    sz=size(A)
    coord=map(x->1:x,sz)
    return SpectData(A,coord)
end

import Base: convert,similar
convert(::Type{SpectData},A::AbstractArray) = SpectData(A)


@doc raw"""
function import(path::String,vendor::Symbol)

Import a data set located at `path`, and return a `DataSet` object. `vendor` designates
the origin data format. Currently implemented are

- `:Bruker`: the path points to a directory with a Bruker NMR data set.
- `:JEOL`: the path points to a JEOL `.jdf` file
"""
function load(f::String,vendor::Symbol)
    if vendor == :Bruker 
        params = FileIO.readBrukerParameterFile(f*"/acqus")
        rawdata = FileIO.readBrukerFID(f*"/fid")
        rawdata = rawdata[params["GRPDLY"]:end]
        tcoord = range(0.0,step=1.0/params["SW_h"],length=length(rawdata))
        return params,SpectData(rawdata,(tcoord,))
    elseif vendor == :JEOL
        error("Not yet implemented")
    else 
        error("Unsupported data format")
    end
end