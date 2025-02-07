
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
for SpectData objects. However, they can also be applied to
any AbstractArray, by promoting it to a SpectData, and then
extracting the data part.
"""
(m::NMRProcessor)(A::AbstractArray) = m(SpectData(A)).dat


@doc raw"""
    Chain(fs::Vararg{NMRProcessor}) 

returns a chain of processing tools, which will be
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

@doc raw"""
    function FourierTransform(SI::Vector,dims::Vector; fftshift=true)

Fourier transform processor for data sets of size `SI`. `dims` is a Vector
of the dimensions along which a Fourier transform will be computed.
The corresponding coordinates are automatically replaced by frequencies,
based on the Nyqvist theorem. The zero frequency appears in the centre of the
spectrum. 
"""
function FourierTransform(SI::Vector,dims::Vector; fftshift=true)
    dummy=zeros(ComplexF64,SI...)
    plan=FFTW.plan_fft(dummy,dims)
    return( FourierTransform(dims,SI,fftshift,plan))
end

function (ft::FourierTransform)(S::SpectData)
    ftdat = ft.plan*S.dat
    newcoord=[]
    for (k,d) in enumerate(S.coord)
        if k in ft.dims && d isa AbstractRange
            Δf = 1.0/step(d)
            push!(newcoord, range(-Δf,Δf,length=length(d)))
        else
            push!(newcoord,d)
        end
    end
    if ft.fftshift
        ftdat=FFTW.fftshift(ftdat,ft.dims)
    end
    return SpectData(ftdat,(newcoord...,))
end


struct ZeroFill <: NMRProcessor
    SI::Vector{Union{Integer,Colon}}
end

function (zf::ZeroFill)(A::SpectData{T,N}) where {T,N}
    oldsize = size(A)
    newsize = zf.SI
    newcoord=[]
    for (k,d) in enumerate(newsize)
        if d isa Colon
            newsize[k]=oldsize[k]
            push!(newcoord,A.coord[k])
        elseif A.coord[k] isa AbstractRange
            extrange = range(first(A.coord[k]),step=step(A.coord[k]),length=d)
            push!(newcoord,extrange)
        else
            push!(newcoord,A.coord[k])
        end
    end
    newA = zeros(T, newsize...)
    oldRange = map(x->1:x,oldsize)
    newA[oldRange...] = A
    return SpectData(newA,(newcoord...,)) 
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

