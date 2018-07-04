


abstract type AbstractRecursiveVector{T,L} <: StaticArrays.StaticArray{Tuple{L}, T, 1} end
abstract type MutableRecursiveVector{T,L} <: AbstractRecursiveVector{T,L} end

struct StaticRecursiveVector{T,L} <: AbstractRecursiveVector{T,L}
    data::NTuple{L,T}
    StaticRecursiveVector{T,L}(data::NTuple{L,T}) where {T,L} = new(data)
end
mutable struct RecursiveVector{T,L} <: MutableRecursiveVector{T,L}
    data::NTuple{L,T}
    RecursiveVector{T,L}(data::NTuple{L,T}) where {T,L} = new(data)
    RecursiveVector{T,L}() where {T,L} = new()
end
struct PointerRecursiveVector{T,M,N,L} <: MutableRecursiveVector{T,L}
    data::Ptr{StaticRecursiveVector{T,cutoff}}
    # PointerRecursiveVector(ptr::Ptr{StaticRecursiveVector{T,cutoff,cutoff,cutoff2}}, ::Val{M}, ::Val{N}, ::Val{L}) where {T,M,N,L} = new{T,M,N,L}(ptr)
end


@inline float_point(A::RecursiveVector{T}) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
@inline float_point(A::PointerRecursiveVector{T}) where T = Base.unsafe_convert(Ptr{T}, A.data)

Base.similar(::RecursiveVector{T,L}) where {T,L} = RecursiveVector{T,L}()
@inline Base.size(::RecursiveVector{T,L}) where {T,L} = (L,)
@inline Base.length(::RecursiveVector{T,L}) where {T,L} = L
@inline function Base.getindex(x::RecursiveVector{T,L}, i::Int) where {T,L}
    @boundscheck i > L && throw(BoundsError())
    x.data[i]
end
@inline function Base.setindex!(x::RecursiveVector{T,L}, val, i::Int) where {T,L}
    @boundscheck i > L && throw(BoundsError())    
    unsafe_store!(float_point(x), convert(T, val), i )
    val
end

function Base.copyto!(x::RecursiveVector{T,L}, y::RecursiveVector{T,L}) where {T,L}
    @inbounds for l ∈ 1:L
        x[l] = y[l]
    end
    x
end
Base.copy(x::RecursiveVector{T,L}) where {T,L} = copyto!(RecursiveVector{T,L}(), x)

function Base.scale!(x::RecursiveVector{T,L}, α) where {T,L}
    @inbounds @simd for l ∈ 1:L
        x[l] *= α
    end
    x
end

function LinearAlgebra.vecdot(x::RecursiveVector{T,L}, y::RecursiveVector{T,L}) where {T,L}


end


