
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
            push!(newcoord, range(-Δf/2,Δf/2,length=length(d)))
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


struct MedianBaselineCorrect <: NMRProcessor
    dim::Int64
    wdw::Int64
    gauss::Vector{Float64}
end

function MedianBaselineCorrect(dim::Integer;wdw=1024)
    g=exp.(-((-wdw:wdw)./wdw).^2)
    g=g/sum(g)
    return MedianBaselineCorrect(dim,wdw,g)
end

@doc """
    function extrema(X::AbstractArray{T,N}, dim::Integer) 

returns an array of all extremal values of `X` along the dimension `dim`.
"""
function extrema(X::AbstractArray{T,N}, dim::Integer) where {T,N} 
    shifter=zeros(Int64,ndims(X))
    shifter[dim]=1
    left=circshift(X,-shifter)
    right=circshift(X,shifter)

    minmaxima = (X .> left .&& X .> right) .|| (X .< left .&& X .< right)
    return X[minmaxima]
end



@doc"""
    function conv(X::AbstractArray{T1,N}, y::AbstractVector{T2},dim::Integer) where {N>1,T1,T2}

computes the convolution of the array `X` with the vector `y` along the dimension `dim`.
The ends of `X` are zero-padded such that the result is guaranteed to have the same size as `X`.
`NMRlab.conv()` uses a direct algorithm for the convolution, not fft. It is therefore efficient
when the length of `y` is much less than the corresponding dimension of `X`. If this is not
the case and performance is critical, a different algorithm should be used.
"""
function conv(X::AbstractArray{T1,N}, y::AbstractVector{T2},dim::Integer) where {N,T1,T2}
    return mapslices(a->conv(a,y), X, dims=dim)
end

function conv(x::AbstractVector{T1}, y::AbstractVector{T2} ) where {T1,T2}
    N=length(x)
    M=length(y)>>1
    lrange = isodd(length(y)) ? (-M:M) : (-M:(M-1))
    b = [ sum( ((k+l>0 && k+l<=N) ? x[k+l] : zero(T1) ) * y[l+M+1] for l=lrange)   for k=1:N]
end


import Statistics

@doc raw"""
    function (mb::MedianBaselineCorrect)(s::SpectData)

subtract baseline for the real part of `s` by the algorithm of M. S. Friedrichs,
*Journal of Biomolecular NMR*,  **5** (1995) 147  153.
"""
function (mb::MedianBaselineCorrect)(s::SpectData)
    r=real(s.dat)
    ax=axes(s,mb.dim)

    indices=Any[axes(s)...]
    b = [ begin
           indices[mb.dim]= max(k-mb.wdw,1):min(k+mb.wdw,length(ax))
           Statistics.median(extrema(view(r,indices...),mb.dim))
         end

        for k in ax
    ]
    e2kernel=exp.( ((-mb.wdw:mb.wdw)./mb.wdw).^2)
    e2kernel ./= sum(e2kernel)
    baseline = conv(cat(b,dims=mb.dim) , e2kernel)
    return s-baseline
end

@doc raw"""
    abstract type NMRProcessor1D <: NMRProcessor

Data type for processing tools that apply to a single dimension, i.e., that are inherently
1D. They need to be defined as functors that act on a `AbstractVector{T}`. They must
contain a field `.dim` that indicates which dimension in a multidimensional array they
should be applied to.
"""
abstract type NMRProcessor1D <: NMRProcessor end

function (np1d::NMRProcessor1D)(A::SpectData{T,N}) where {T,N}
    return mapslices(np1d,A,dims=np1d.dim)
end


struct PhaseCorrect <: NMRProcessor1D
    ph0::Float64
    ph1::Float64
    dim::Int32
end

function (pc::PhaseCorrect)(x::SpectData{T,1}) where {T<:Number}
    c = exp(im*pc.ph0).* exp.(im*pc.ph1.*coords(x,1))
    return x.*c 
end
