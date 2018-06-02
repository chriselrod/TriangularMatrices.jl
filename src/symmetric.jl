abstract type AbstractSymmetricMatrix{T,N,N2} <: AbstractMatrix{T} end
abstract type MutableSymetricMatrix{T,N,N2} <: AbstractSymmetricMatrix{T,N,N2} end

"""
Actually Hermitian, because adjoint returns itself.
"""
struct StaticSymmetricMatrix{T,N,N2} <: AbstractSymmetricMatrix{T,N,N2}
    data::NTuple{N2,T}
end
struct SymmetricMatrix{T,N,N2} <: MutableSymetricMatrix{T,N,N2}
    data::Base.RefValue{NTuple{N2,T}}
    function SymmetricMatrix{T, N, N2}(d::Base.RefValue{NTuple{N2,T}}) where {T,N,N2}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,N2}(d)
    end
end
struct PointerSymmetricMatrix{T,N,N2} <: MutableSymetricMatrix{T,N,N2}
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
    push!(q.args, :(SymmetricMatrix{$T,$N,$N2}( Ref( @ntuple $N2 A ))))
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
SymmetricMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = SymmetricMatrix{T,N,N2}(Ref(data))
StaticSymmetricMatrix(data::SymmetricMMatrix{T,N,N2}) where {T,N,N2} = StaticSymmetricMatrix{T,N,N2}(data.data[])
SymmetricMatrix(data::SymmetricMatrix{T,N,N2}) where {T,N,N2} = SymmetricMatrix{T,N,N2}(Ref(data.data))


@inline Base.getindex(A::StaticSymmetricMatrix, i::Integer) = A.data[i]
@inline function Base.getindex(A::SymmetricMatrix{T}, i::Integer) where T
    unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), i)
end
@inline Base.getindex(A::PointerSymmetricMatrix, i::Integer) where T = unsafe_load(A.data, i)
@inline function Base.getindex(A::StaticSymmetricMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    @inbounds out = A.data[small_triangle(j)+i]
    out
end
@inline function Base.getindex(A::SymmetricMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    # @inbounds out = A.data[][small_triangle(j)+i]
    # out
    unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), small_triangle(j)+i)
end
@inline function Base.getindex(A::PointerSymmetricMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    # @inbounds out = A.data[][small_triangle(j)+i]
    # out
    unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), small_triangle(j)+i)
end


@generated function (A::SymmetricMatrix{T,N,N2})(i::Integer, j::Integer) where {T,N,N2}
    num_blocks = blocks_per_dim(N)
    quote
        $(Expr(:meta, :inline))
        @boundscheck begin
            max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
        end
        i, j = minmax(i, j) # j >= i
        row_block_id = $(1+num_blocks) - cld($num_blocks*($(N+1)-i),$N)
        col_block_id = $(1+num_blocks) - cld($num_blocks*($(N+1)-j),$N)

        row_block_lbound = cld((row_block_id-1)*$N,$num_blocks)
        # row_block_ubound = cld( row_block_id   *$N,$num_blocks)
        col_block_lbound = cld((col_block_id-1)*$N,$num_blocks)
        col_block_ubound = cld( col_block_id   *$N,$num_blocks)

        ind = big_triangle(col_block_lbound) + (col_block_ubound - col_block_lbound) * row_block_lbound

        block_i = i - row_block_lbound
        block_j = j - col_block_lbound

        if row_block_id == col_block_id # then we're in a symmetric portion
            ind += ltriangle(block_j) + block_i
        else #then we are in square block
            nrows_in_block = cld( row_block_id   *$N,$num_blocks)  - row_block_lbound
            ind += nrows_in_block * (block_j-1) + block_i
        end
        unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), ind)
    end
end
@generated function (A::PointerSymmetricMatrix{T,N,N2})(i::Integer, j::Integer) where {T,N,N2}
    num_blocks = blocks_per_dim(N)
    quote
        $(Expr(:meta, :inline))
        @boundscheck begin
            max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
        end
        i, j = minmax(i, j) # j >= i
        row_block_id = $(1+num_blocks) - cld($num_blocks*($(N+1)-i),$N)
        col_block_id = $(1+num_blocks) - cld($num_blocks*($(N+1)-j),$N)

        row_block_lbound = cld((row_block_id-1)*$N,$num_blocks)
        # row_block_ubound = cld( row_block_id   *$N,$num_blocks)
        col_block_lbound = cld((col_block_id-1)*$N,$num_blocks)
        col_block_ubound = cld( col_block_id   *$N,$num_blocks)

        ind = big_triangle(col_block_lbound) + (col_block_ubound - col_block_lbound) * row_block_lbound

        block_i = i - row_block_lbound
        block_j = j - col_block_lbound

        if row_block_id == col_block_id # then we're in a symmetric portion
            ind += ltriangle(block_j) + block_i
        else #then we are in square block
            nrows_in_block = cld( row_block_id   *$N,$num_blocks)  - row_block_lbound
            ind += nrows_in_block * (block_j-1) + block_i
        end
        unsafe_load(A.data, ind)
    end
end


@inline function Base.setindex!(A::SymmetricMatrix{T,N,N2}, val, i::Integer) where  {T,N,N2}
    @boundscheck begin
        i > N && throw(BoundsError())
    end
    unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), convert(T, val), i)
    return val
end
@inline function Base.setindex!(A::SymmetricMatrix{T,N,N2}, val, i::Integer, j::Integer) where  {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
    end
    i, j = minmax(i, j)
    unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), convert(T, val), small_triangle(j)+i)
    return val
end

LinearAlgebra.transpose(x::SymmetricMatrix) = x
@static if VERSION >= v"0.7.0-DEV"
    LinearAlgebra.adjoint(x::SymmetricMatrix) = x
end
