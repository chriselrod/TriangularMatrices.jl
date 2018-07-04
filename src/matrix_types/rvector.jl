


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
function Base.getindex(x::RecursiveVector{N,T}, i, ::Val{R}) where {N,R,T}
    @boundscheck i > N && throw(BoundsError())
    # ntuple(j -> x.data[i+j-1].value, Val(R))
    RT = Ptr{StaticSIMD{R,T}}
    Base.unsafe_load(Base.unsafe_convert(RT, pointer_from_objref(x)), i)
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

function LinearAlgebra.scale!(x::RecursiveVector{T,L}, α) where {T,L}
    @inbounds @simd for l ∈ 1:L
        x[l] *= α
    end
    x
end

Base.vec(x::RecursiveVector) = x

struct StaticSIMD{N,T} <: AbstractVector{T}
    data::NTuple{N,Core.VecElement{T}}
end
Base.size(::StaticSIMD{N}) where N = (N,)
Base.length(::StaticSIMD{N}) where N = N
function Base.getindex(x::StaticSIMD{N}, i::Int) where N
    @boundscheck i > N && throw(BoundsError())
    x.data[i].value
end
@generated function Base.:*(x::StaticSIMD{N,T}, y::StaticSIMD{N,T}) where {N,T}
    quote
        Base.@_inline_meta
        StaticSIMD{$N,$T}( $(Expr(:tuple, [:(x.data[$i].value * y.data[$i].value) for i ∈ 1:N]...)) )
    end
end
@generated function Base.:+(x::StaticSIMD{N,T}, y::StaticSIMD{N,T}) where {N,T}
    quote
        Base.@_inline_meta
        StaticSIMD{$N,$T}( $(Expr(:tuple, [:(x.data[$i].value + y.data[$i].value) for i ∈ 1:N]...)) )
    end
end
@generated function Base.fma(x::StaticSIMD{N,T}, y::StaticSIMD{N,T}, z::StaticSIMD{N,T}) where {N,T}
    quote
        Base.@_inline_meta
        StaticSIMD{$N,$T}( $(Expr(:tuple, [:(x.data[$i].value * y.data[$i].value + z.data[$i].value) for i ∈ 1:N]...)) )
    end
end
@generated function Base.sum(x::StaticSIMD{N,T}) where {N,T}
    quote
        Base.@_inline_meta
        $(Expr(:call, :+, [:(x.data[$i].value) for i ∈ 1:N]...))
    end
end
@generated function LinearAlgebra.dot(x::RecursiveVector{T,L}, y::RecursiveVector{T,L}) where {T,L}
    vec_length = 8
    VVal = Val(vec_length)
    Lo4, r = divrem(L, vec_length)
    quote
        # @inbounds out_tup = x[1:4] .* y[1:4]
        @inbounds out_tup = x[1, $VVal] * y[1, $VVal]
        @fastmath for i ∈ 2:$Lo4
            @inbounds out_tup = fma(x[i, $VVal], y[i, $VVal], out_tup)
        end
        out = sum(out_tup)
        @inbounds for i ∈ $(L-r+1):$L
            out += x[i] * y[i]
        end
        out
    end
end


@generated function abs2norm(x::RecursiveVector{T,L}) where {T,L}
    vec_length = 8
    VVal = Val(vec_length)
    Lo4, r = divrem(L, vec_length)
    quote
        # @inbounds out_tup = x[1:4] .* y[1:4]
        simdv = x[1, $VVal]
        @inbounds out_tup = simdv * simdv
        @fastmath for i ∈ 2:$Lo4
            simdv = x[i, $VVal]
            @inbounds out_tup = fma(simdv, simdv, out_tup)
        end
        out = sum(out_tup)
        @inbounds for i ∈ $(L-r+1):$L
            out += abs2(x[i])
        end
        out
    end
end


