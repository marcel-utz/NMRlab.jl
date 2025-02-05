
@doc raw"""
    abstract type NMRProcessor <: Function end

Abstract data type for NMR processing routines. 
New processing tools should be declared as subtypes of
`NMRProcessor`.
"""
abstract type NMRProcessor <: Function end

@doc raw"""
    (m::NMRProcessor)(A::AbstractArray)

here is some fallback behaviour. NMR processors are defined
# for SpectData objects. However, they can also be applied to
# any AbstractArray, by promoting it to a SpectData, and then
# extracting the data part.
"""
(m::NMRProcessor)(A::AbstractArray) = m(SpectData(A)).dat


@doc raw"""
    Chain(fs::Vararg{NMRProcessor}) 

create a chain of processing tools, which will be
applied in order (the first in the argument list is 
applied first)
"""
Chain(fs::Vararg{NMRProcessor}) = reduce(∘, reverse(fs))


import FFTW

struct FourierTransform <: NMRProcessor
    dims::Vector{Integer}
    SI::Vector{Integer}
    fftshift::Bool
    plan
end

function FourierTransform(SI::Vector,dims::Vector; fftshift=true)
    dummy=zeros(ComplexF64,SI...)
    plan=FFTW.plan_fft(dummy,dims)
    return( FourierTransform(dims,SI,fftshift,plan))
end

function (ft::FourierTransform)(S::SpectData)
    ftdat = ft.plan*S.dat
    if ft.fftshift
        FFTW.fftshift!(ftdat,dims)
    end
    return SpectData(ftdat,S.coord)
end


struct ZeroFill <: NMRProcessor
    SI::Vector{Union{Integer,Colon}}
end

function (zf::ZeroFill)(A::SpectData{T,N}) where {T,N}
    oldsize = size(A)
    newsize = zf.SI
    for (k,d) in enumerate(newsize)
        if d isa Colon
            newsize[k]=oldsize[k]
        end
    end
    newA = zeros(T, newsize...)
    oldRange = map(x->1:x,oldsize)
    newA[oldRange...] = A
    return SpectData(newA,A.coord) 
end

struct Apodize <: NMRProcessor
    R::Vector{Union{Real,Colon}}
end

function (ap::Apodize)(A::SpectData)
    apo = A.dat
    for (k,r) in enumerate(ap.R)
        if !(r isa Colon)
            ix=ones(Int64,ndims(apo))
            ix[k] = length(A.coord[k])
            f = exp.(-ap.R[k] .* A.coord[k])
            apo .*= reshape(f,ix...)
        end
    end

    return(SpectData(apo,A.coord))
end

