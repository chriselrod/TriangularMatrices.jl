abstract type AbstractSymmetricMatrix{T,N,L} <: AbstractMatrix{T} end
abstract type MutableSymmetricMatrix{T,N,L} <: AbstractSymmetricMatrix{T,N,L} end

"""
Actually Hermitian, because adjoint returns itself.
"""
struct StaticSymmetricMatrix{T,N,L} <: AbstractSymmetricMatrix{T,N,L}
    data::NTuple{L,T}
end
mutable struct SymmetricMatrix{T,N,L} <: MutableSymmetricMatrix{T,N,L}
    data::NTuple{L,T}
    @generated SymmetricMatrix{T,N}() where {T,N} = :(SymmetricMatrix{$T,$N,$(N*N)}())
    function SymmetricMatrix{T,N,L}() where {T,N,L}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,L}()
    end
    function SymmetricMatrix{T, N, L}(d::NTuple{L,T}) where {T,N,L}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,L}(d)
    end
end
struct PointerSymmetricMatrix{T,N,L} <: MutableSymmetricMatrix{T,N,L}
    data::Ptr{T}
end
Base.size(::AbstractSymmetricMatrix{T,N}) where {T,N} = (N,N)


function StaticSymmetricMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    StaticSymmetricMatrix(A, Val(n))
end
@generated function StaticSymmetricMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, L = upper_triangle_quote(N, T)
    push!(q.args, :(StaticSymmetricMatrix{$T,$N,$L}( @ntuple $L A )))
    q
end
function SymmetricMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    SymmetricMatrix(A, Val(n))
end
@generated function SymmetricMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, L = upper_triangle_quote(N, T)
    push!(q.args, :(SymmetricMatrix{$T,$N,$L}( data )))
    q
end
function StaticSymmetricMatrix(data::AbstractVector)
    L = length(data)
    N = inv_triangle(L)
    StaticSymmetricMatrix(data, Val(N), Val(L))
end
function StaticSymmetricMatrix(data::AbstractVector{T}, ::Val{N}, ::Val{L}) where {T,N,L}
    StaticSymmetricMatrix{T,N,L}(ntuple(i -> data[i], L))
end
StaticSymmetricMatrix(data::NTuple{L,T}) where {T,L} = StaticSymmetricMatrix(data,ValI(Val{L}()))
SymmetricMatrix(data::NTuple{L,T}) where {T,L} = SymmetricMatrix(data,ValI(Val{L}()))
StaticSymmetricMatrix(data::SVector{L,T}) where {T,L} = StaticSymmetricMatrix(data.data,ValI(Val{L}()))
SymmetricMatrix(data::SVector{L,T}) where {T,L} = SymmetricMatrix(data.data,ValI(Val{L}()))
StaticSymmetricMatrix(data::NTuple{L,T}, ::Val{N}) where {T,N,L} = StaticSymmetricMatrix{T,N,L}(data)
SymmetricMatrix(data::NTuple{L,T}, ::Val{N}) where {T,N,L} = SymmetricMatrix{T,N,L}(data)
StaticSymmetricMatrix(data::SymmetricMatrix{T,N,L}) where {T,N,L} = StaticSymmetricMatrix{T,N,L}(data.data)
SymmetricMatrix(data::StaticSymmetricMatrix{T,N,L}) where {T,N,L} = SymmetricMatrix{T,N,L}(data.data)

point(A::MutableSymmetricMatrix{T}) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
point(A::PointerSymmetricMatrix) = A.data


@inline Base.getindex(A::StaticSymmetricMatrix, i::Integer) = A.data[i]
@inline Base.getindex(A::SymmetricMatrix, i::Integer) = A.data[i]
@inline Base.getindex(A::PointerSymmetricMatrix, i::Integer) = unsafe_load(A.data, i)
@inline function Base.getindex(A::StaticSymmetricMatrix{T,N,L}, i::Integer, j::Integer) where {T,N,L}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    @inbounds out = A.data[small_triangle(j)+i]
    out
end

# @inline function Base.getindex(A::MutableSymmetricMatrix{T,N,L}, i::Integer, j::Integer) where {T,N,L}
#     @boundscheck begin
#         if (max(i,j) > 2) || (min(i,j) < 1)
#             throw(BoundsError("When recursing only indicies of 1 and 2 are supported. Received $i,$j.")) 
#         end
#     end
#     i, j = minmax(i, j)
#     unsafe_load(point(A.data), small_triangle(j)+i)
# end


# for N ∈ 1:cutoff
#     L = big_triangle(N)

#     @eval Base.size(A::MutableSymmetricMatrix{T,$N,$L}) = ($N,$N)

#     @eval @inline function Base.getindex(A::MutableSymmetricMatrix{T,$N,$L}, i::Integer, j::Integer) where {T,$N,$L}
#         @boundscheck begin
#             max(i,j) > $N && throw("BoundsError: index ($i, $j) is out of bounds.")
#         end
#         i, j = minmax(i, j)
#         unsafe_load(point(A.data), small_triangle(j)+i)
#     end

# end


"""
Note that Cartesian indexing is relatively slow, and results in branches even with `@inbounds`.
This means that `@simd`, for example, will not work while using Cartesian indices.
While linear indexing is fast, note that the underlying data is neither column nor row-major, so
be careful that your computations are actually doing what you intend.

You can use triangle_sub2ind (which incurs branching) to explore the behavior of
Cartesian -> Linear indices for Symmetric and Triangular matrices.
"""
@inline function Base.getindex(A::MutableSymmetricMatrix{T,N}, i::Integer, j::Integer) where {T,N}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j) # j >= i
    A[triangle_sub2ind(Val{N}(), i, j)]
end
@generated function Base.getindex(A::MutableSymmetricMatrix{T,N}, ::Val{1}, ::Val{1}) where {T, N}
    Nhalf = cld(N,2)
    triangle_size = big_triangle(Nhalf)
    :(PointerSymmetricMatrix{$T,$Nhalf,$triangle_size}(point(A)))
end
@generated function Base.getindex(A::MutableSymmetricMatrix{T,N}, ::Val{2}, ::Val{1}) where {T, N}
    Nremain, r = divrem(N,2)
    Nhalf = Nremain + r
    L = Nhalf * Nremain
    pointer_offset = big_triangle(Nhalf) * sizeof(T)
    :(LinearAlgebra.adjoint(PointerMatrix{$T,$Nhalf,$Nremain,$L}(point(A) + $pointer_offset )))
end
@generated function Base.getindex(A::MutableSymmetricMatrix{T,N}, ::Val{1}, ::Val{2}) where {T, N}
    Nremain, r = divrem(N,2)
    Nhalf = Nremain + r
    L = Nhalf * Nremain
    pointer_offset = big_triangle(Nhalf) * sizeof(T)
    :(PointerMatrix{$T,$Nhalf,$Nremain,$L}(point(A) + $pointer_offset ))
end
@generated function Base.getindex(A::MutableSymmetricMatrix{T,N,L}, ::Val{2}, ::Val{2}) where {T, N, L}
    Nremain = N ÷ 2
    triangle_size = big_triangle(Nremain)
    pointer_offset = (L - triangle_size) * sizeof(T)
    :(PointerSymmetricMatrix{$T,$Nremain,$triangle_size}(point(A) + $pointer_offset ))
end


@inline function Base.setindex!(A::MutableSymmetricMatrix{T,N,L}, val, i::Integer) where  {T,N,L}
    @boundscheck begin
        i > L && throw(BoundsError())
    end
    unsafe_store!(point(A), convert(T, val), i)
    return val
end
@inline function Base.setindex!(A::SymmetricMatrix{T,N,L}, val, i::Integer, j::Integer) where  {T,N,L}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
    end
    i, j = minmax(i, j)
    unsafe_store!(point(A), convert(T, val), triangle_sub2ind(Val{N}(), i, j))
    return val
end

LinearAlgebra.transpose(x::SymmetricMatrix) = x
@static if VERSION >= v"0.7.0-DEV"
    LinearAlgebra.adjoint(x::SymmetricMatrix) = x
end
