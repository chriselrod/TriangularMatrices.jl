abstract type AbstractSymmetricMatrix{T,N,N2} <: AbstractMatrix{T} end
abstract type MutableSymmetricMatrix{T,N,N2} <: AbstractSymmetricMatrix{T,N,N2} end

"""
Actually Hermitian, because adjoint returns itself.
"""
struct StaticSymmetricMatrix{T,N,N2} <: AbstractSymmetricMatrix{T,N,N2}
    data::NTuple{N2,T}
end
mutable struct SymmetricMatrix{T,N,N2} <: MutableSymmetricMatrix{T,N,N2}
    data::NTuple{N2,T}
    function SymmetricMatrix{T, N, N2}(d::NTuple{N2,T}) where {T,N,N2}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,N2}(d)
    end
end
struct PointerSymmetricMatrix{T,N,N2} <: MutableSymmetricMatrix{T,N,N2}
    data::Ptr{T}
end
Base.size(::AbstractSymmetricMatrix{T,N}) where {T,N} = (N,N)


function StaticSymmetricMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    StaticSymmetricMatrix(A, Val(n))
end
@generated function StaticSymmetricMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = upper_triangle_quote(N, T)
    push!(q.args, :(StaticSymmetricMatrix{$T,$N,$N2}( @ntuple $N2 A )))
    q
end
function SymmetricMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    SymmetricMatrix(A, Val(n))
end
@generated function SymmetricMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = upper_triangle_quote(N, T)
    push!(q.args, :(SymmetricMatrix{$T,$N,$N2}( data )))
    q
end
function StaticSymmetricMatrix(data::AbstractVector)
    N2 = length(data)
    N = inv_triangle(N2)
    StaticSymmetricMatrix(data, Val(N), Val(N2))
end
function StaticSymmetricMatrix(data::AbstractVector{T}, ::Val{N}, ::Val{N2}) where {T,N,N2}
    StaticSymmetricMatrix{T,N,N2}(ntuple(i -> data[i], N2))
end
StaticSymmetricMatrix(data::NTuple{N2,T}) where {T,N2} = StaticSymmetricMatrix(data,ValI(Val{N2}()))
SymmetricMatrix(data::NTuple{N2,T}) where {T,N2} = SymmetricMatrix(data,ValI(Val{N2}()))
StaticSymmetricMatrix(data::SVector{N2,T}) where {T,N2} = StaticSymmetricMatrix(data.data,ValI(Val{N2}()))
SymmetricMatrix(data::SVector{N2,T}) where {T,N2} = SymmetricMatrix(data.data,ValI(Val{N2}()))
StaticSymmetricMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = StaticSymmetricMatrix{T,N,N2}(data)
SymmetricMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = SymmetricMatrix{T,N,N2}(data)
StaticSymmetricMatrix(data::SymmetricMatrix{T,N,N2}) where {T,N,N2} = StaticSymmetricMatrix{T,N,N2}(data.data)
SymmetricMatrix(data::StaticSymmetricMatrix{T,N,N2}) where {T,N,N2} = SymmetricMatrix{T,N,N2}(data.data)

point(A::MutableSymmetricMatrix{T}) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
point(A::PointerSymmetricMatrix) = A.data


@inline Base.getindex(A::StaticSymmetricMatrix, i::Integer) = A.data[i]
@inline Base.getindex(A::SymmetricMatrix, i::Integer) = A.data[i]
@inline Base.getindex(A::PointerSymmetricMatrix, i::Integer) = unsafe_load(A.data, i)
@inline function Base.getindex(A::StaticSymmetricMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    @inbounds out = A.data[small_triangle(j)+i]
    out
end

# @inline function Base.getindex(A::MutableSymmetricMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
#     @boundscheck begin
#         if (max(i,j) > 2) || (min(i,j) < 1)
#             throw(BoundsError("When recursing only indicies of 1 and 2 are supported. Received $i,$j.")) 
#         end
#     end
#     i, j = minmax(i, j)
#     unsafe_load(point(A.data), small_triangle(j)+i)
# end


# for N ∈ 1:cutoff
#     N2 = big_triangle(N)

#     @eval Base.size(A::MutableSymmetricMatrix{T,$N,$N2}) = ($N,$N)

#     @eval @inline function Base.getindex(A::MutableSymmetricMatrix{T,$N,$N2}, i::Integer, j::Integer) where {T,$N,$N2}
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
@generated function Base.getindex(A::MutableSymmetricMatrix{T,N,N2}, ::Val{2}, ::Val{2}) where {T, N, N2}
    Nremain = N ÷ 2
    triangle_size = big_triangle(Nremain)
    pointer_offset = (N2 - triangle_size) * sizeof(T)
    :(PointerSymmetricMatrix{$T,$Nremain,$triangle_size}(point(A) + $pointer_offset ))
end


@inline function Base.setindex!(A::MutableSymmetricMatrix{T,N,N2}, val, i::Integer) where  {T,N,N2}
    @boundscheck begin
        i > N && throw(BoundsError())
    end
    unsafe_store!(point(A), convert(T, val), i)
    return val
end
@inline function Base.setindex!(A::SymmetricMatrix{T,N,N2}, val, i::Integer, j::Integer) where  {T,N,N2}
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
