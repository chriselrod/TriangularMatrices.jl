"""
Upper triangular and lower triangular matrices are the same views of the underlying data.
Upper triangular matrices are column major, and lower triangular row major.

As an implementation detail, indexing with @inbounds into the 0-triangle (ie, lower triangle of an Upper Triangular matrix) is undefined. This is because that check creates a branch.
"""
abstract type AbstractTriangularMatrix{T, N, N2} <: AbstractMatrix{T} end

abstract type AbstractUpperTriangular{T, N, N2} <: AbstractTriangularMatrix{T, N, N2} end
abstract type AbstractLowerTriangular{T, N, N2} <: AbstractTriangularMatrix{T, N, N2} end
abstract type MutableUpperTriangular{T, N, N2} <: AbstractUpperTriangular{T, N, N2} end
abstract type MutableLowerTriangular{T, N, N2} <: AbstractLowerTriangular{T, N, N2} end

const MutableTriangularMatrix{T,N,N2} = Union{MutableUpperTriangular{T, N, N2}, MutableLowerTriangular{T, N, N2}}

struct StaticUpperTriangularMatrix{T, N, N2} <: AbstractUpperTriangular{T, N, N2}
    data::NTuple{N2,T}
end
struct StaticLowerTriangularMatrix{T, N, N2} <: AbstractLowerTriangular{T, N, N2}
    data::NTuple{N2,T}
end
mutable struct UpperTriangularMatrix{T, N, N2} <: MutableUpperTriangular{T, N, N2}
    data::NTuple{N2,T}
    function UpperTriangularMatrix{T, N, N2}(d::NTuple{N2,T}) where {T,N,N2}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,N2}(d)
    end
end
mutable struct LowerTriangularMatrix{T, N, N2} <: MutableLowerTriangular{T, N, N2}
    data::NTuple{N2,T}
    function LowerTriangularMatrix{T, N, N2}(d::NTuple{N2,T}) where {T,N,N2}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,N2}(d)
    end
end
struct PointerUpperTriangularMatrix{T, N, N2} <: MutableUpperTriangular{T, N, N2}
    data::Ptr{T}
end
struct PointerLowerTriangularMatrix{T, N, N2} <: MutableLowerTriangular{T, N, N2}
    data::Ptr{T}
end
Base.size(::AbstractTriangularMatrix{T,N}) where {T,N} = (N,N)

const PointerTriangularMatrix{T,N,N2} = Union{
    PointerUpperTriangularMatrix{T, N, N2},
    PointerLowerTriangularMatrix{T, N, N2}
}

point(A::MutableTriangularMatrix{T}) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
point(A::PointerTriangularMatrix) = A.data

function upper_triangle_quote(N, S = :A)
    N2 = big_triangle(N)
    q = quote
        @inbounds begin
        end
    end
    @static if VERSION > v"0.7-"
        qa = q.args[2].args[3].args
    else
        qa = q.args[2].args[2].args
    end
    ind = 0
    for i ∈ 1:N, j ∈ 1:i
        ind += 1
        A_i = Symbol(S, :_, ind)
        push!(qa, :( $A_i = $(S)[$j,$i] ) )
    end
    q, N2
end
function lower_triangle_quote(N, S = :A)
    N2 = big_triangle(N)
    q = quote
        @inbounds begin
        end
    end
    @static if VERSION > v"0.7-"
        qa = q.args[2].args[3].args
    else
        qa = q.args[2].args[2].args
    end
    for i ∈ 1:N, j ∈ i:N
        A_i = Symbol(S, :_, small_triangle(j)+i )
        push!(qa, :( $A_i = $(S)[$j,$i] ) )
    end
    push!(q.args, :($M{$T,$N,$N2}( @ntuple $N2 A )))
    q, N2
end

function recursive_lower_triangle_quote!(qa, N, ::Type{T}) where T

end


function UpperTriangularMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    UpperTriangularMatrix(A, Val(n))
end
@generated function UpperTriangularMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = upper_triangle_quote(N, T)
    push!(q.args, :(UpperTriangularMatrix{$T,$N,$N2}( @ntuple $N2 A )))
    q
end
function StaticUpperTriangularMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    StaticUpperTriangularMatrix(A, Val(n))
end
@generated function StaticUpperTriangularMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = upper_triangle_quote(N, T)
    push!(q.args, :(StaticUpperTriangularMatrix{$T,$N,$N2}( @ntuple $N2 A )))
    q
end
function LowerTriangularMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    LowerTriangularMatrix(A, Val(n))
end
@generated function LowerTriangularMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = lower_triangle_quote(N, T)
    push!(q.args, :(LowerTriangularMatrix{$T,$N,$N2}( @ntuple $N2 A )))
    q
end
function StaticLowerTriangularMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    StaticLowerTriangularMatrix(A, Val(n))
end
@generated function StaticLowerTriangularMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = lower_triangle_quote(N, T)
    push!(q.args, :(StaticLowerTriangularMatrix{$T,$N,$N2}( @ntuple $N2 A )))
    q
end
UpperTriangularMatrix(data::NTuple{N2,T}) where {T,N2} = UpperTriangularMatrix(data,ValI(Val{N2}()))
StaticUpperTriangularMatrix(data::NTuple{N2,T}) where {T,N2} = StaticUpperTriangularMatrix(data,ValI(Val{N2}()))
UpperTriangularMatrix(data::SVector{N2,T}) where {T,N2} = UpperTriangularMatrix(data.data,ValI(Val{N2}()))
StaticUpperTriangularMatrix(data::SVector{N2,T}) where {T,N2} = StaticUpperTriangularMatrix(data.data,ValI(Val{N2}()))
UpperTriangularMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = UpperTriangularMatrix{T,N,N2}(data)
StaticUpperTriangularMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = StaticUpperTriangularMatrix{T,N,N2}(data)
UpperTriangularMatrix(data::StaticUpperTriangularMatrix{T,N,N2}) where {T,N,N2} = UpperTriangularMatrix{T,N,N2}(data.data)
StaticUpperTriangularMatrix(data::UpperTriangularMatrix{T,N,N2}) where {T,N,N2} = StaticUpperTriangularMatrix{T,N,N2}(data.data)

LowerTriangularMatrix(data::NTuple{N2,T}) where {T,N2} = LowerTriangularMatrix(data,ValI(Val{N2}()))
StaticLowerTriangularMatrix(data::NTuple{N2,T}) where {T,N2} = StaticLowerTriangularMatrix(data,ValI(Val{N2}()))
LowerTriangularMatrix(data::SVector{N2,T}) where {T,N2} = LowerTriangularMatrix(data.data,ValI(Val{N2}()))
StaticLowerTriangularMatrix(data::SVector{N2,T}) where {T,N2} = StaticLowerTriangularMatrix(data.data,ValI(Val{N2}()))
LowerTriangularMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = LowerTriangularMatrix{T,N,N2}(data)
StaticLowerTriangularMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = StaticLowerTriangularMatrix{T,N,N2}(data)
LowerTriangularMatrix(data::StaticLowerTriangularMatrix{T,N,N2}) where {T,N,N2} = LowerTriangularMatrix{T,N,N2}(data.data)
StaticLowerTriangularMatrix(data::LowerTriangularMatrix{T,N,N2}) where {T,N,N2} = StaticLowerTriangularMatrix{T,N,N2}(data.data)

@inline Base.getindex(A::AbstractTriangularMatrix, i::Integer) = A.data[i]
@inline Base.getindex(A::PointerTriangularMatrix, i::Integer) = unsafe_load(A.data, i)



@inline function Base.getindex(A::StaticUpperTriangularMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i > j
            return zero(T)
        end
    end
    @inbounds out = A.data[small_triangle(j)+i]
    out
end
@inline function Base.getindex(A::StaticLowerTriangularMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i < j
            return zero(T)
        end
    end
    @inbounds out = A.data[small_triangle(i)+j]
    out
end
@inline function Base.getindex(A::UpperTriangularMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i > j
            return zero(T)
        end
    end
    # @inbounds out = A.data[][small_triangle(j)+i]
    # out
    A.data[triangle_sub2ind(Val{N}(), i, j)]
end
@inline function Base.getindex(A::LowerTriangularMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i < j
            return zero(T)
        end
    end
    # @inbounds out = A.data[][small_triangle(i)+j]
    # out
    # unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), small_triangle(i)+j)
    A.data[triangle_sub2ind(Val{N}(), j, i)]
end
@inline function Base.setindex!(A::MutableTriangularMatrix{T,N,N2}, val, i::Integer) where  {T,N,N2}
    @boundscheck i > N2 && throw(BoundsError())
    unsafe_store!(point(A), convert(T, val), i)
    return val
end
@inline function Base.setindex!(A::MutableUpperTriangular{T,N,N2}, val, i::Integer, j::Integer) where  {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i > j
            throw("Cannot set index within lower triangular of UpperTriangularMatrix.")
        end
    end
    unsafe_store!(point(A), convert(T, val), triangle_sub2ind(Val{N}(), i, j))
    return val
end
@inline function Base.setindex!(A::MutableLowerTriangular{T,N,N2}, val, i::Integer, j::Integer) where  {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i < j
            throw("Cannot set index within upper triangular of LowerTriangularMatrix.")
        end
    end
    unsafe_store!(point(A), convert(T, val), triangle_sub2ind(Val{N}(), j, i))
    return val
end

"""
This always returns the static form so that it is non-allocating.
"""
LinearAlgebra.transpose(x::AbstractUpperTriangular{T,N,N2}) where {T,N,N2} = LowerTriangularMatrix{T,N,N2}(x.data)
LinearAlgebra.transpose(x::AbstractLowerTriangular{T,N,N2}) where {T,N,N2} = UpperTriangularMatrix{T,N,N2}(x.data)
@static if VERSION > v"0.6.9"
    """
    Taking the adjoint returns a static matrix so that it is non-allocating.
    """
    LinearAlgebra.adjoint(x::AbstractUpperTriangular{T,N,N2}) where {T<:Real,N,N2} = LowerTriangularMatrix{T,N,N2}(x.data)
    LinearAlgebra.adjoint(x::AbstractLowerTriangular{T,N,N2}) where {T<:Real,N,N2} = UpperTriangularMatrix{T,N,N2}(x.data)
end
