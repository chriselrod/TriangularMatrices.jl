"""
Upper triangular and lower triangular matrices are the same views of the underlying data.
Upper triangular matrices are column major, and lower triangular row major.

As an implementation detail, indexing with @inbounds into the 0-triangle (ie, lower triangle of an Upper Triangular matrix) is undefined. This is because that check creates a branch.
"""
abstract type TriangularMatrix{T, N, N2} <: AbstractMatrix{T} end

abstract type AbstractUpperTriangular{T, N, N2} <: TriangularMatrix{T, N, N2} end
abstract type AbstractLowerTriangular{T, N, N2} <: TriangularMatrix{T, N, N2} end

struct UpperTriangularMatrix{T, N, N2} <: AbstractUpperTriangular{T, N, N2}
    data::NTuple{N2,T}
end
struct LowerTriangularMatrix{T, N, N2} <: AbstractLowerTriangular{T, N, N2}
    data::NTuple{N2,T}
end
struct UpperTriangularMMatrix{T, N, N2} <: AbstractUpperTriangular{T, N, N2}
    data::Base.RefValue{NTuple{N2,T}}
    function UpperTriangularMMatrix{T, N, N2}(d::Base.RefValue{NTuple{N2,T}}) where {T,N,N2}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,N2}(d)
    end
end
struct LowerTriangularMMatrix{T, N, N2} <: AbstractLowerTriangular{T, N, N2}
    data::Base.RefValue{NTuple{N2,T}}
    function LowerTriangularMMatrix{T, N, N2}(d::Base.RefValue{NTuple{N2,T}}) where {T,N,N2}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,N2}(d)
    end
end
Base.size(::TriangularMatrix{T,N}) where {T,N} = (N,N)

function upper_triangle_quote(N, ::Type{T}) where T
    N2 = big_triangle(N)
    q = quote
        @inbounds begin
        end
    end
    @static if VERSION > v"0.6.9"
        qa = q.args[2].args[3].args
    else
        qa = q.args[2].args[2].args
    end
    ind = 0
    for i ∈ 1:N, j ∈ 1:i
        ind += 1
        A_i = Symbol(:A_, ind)
        push!(qa, :( $A_i = A[$j,$i] ) )
    end
    q, N2
end
function lower_triangle_quote(N, ::Type{T}) where T
    N2 = big_triangle(N)
    q = quote
        @inbounds begin
        end
    end
    @static if VERSION > v"0.6.9"
        qa = q.args[2].args[3].args
    else
        qa = q.args[2].args[2].args
    end
    for i ∈ 1:N, j ∈ i:N
        A_i = Symbol(:A_, small_triangle(j)+i )
        push!(qa, :( $A_i = A[$j,$i] ) )
    end
    push!(q.args, :($M{$T,$N,$N2}( @ntuple $N2 A )))
    q, N2
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
function UpperTriangularMMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    UpperTriangularMMatrix(A, Val(n))
end
@generated function UpperTriangularMMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = upper_triangle_quote(N, T)
    push!(q.args, :(UpperTriangularMMatrix{$T,$N,$N2}( Ref( @ntuple $N2 A ))))
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
function LowerTriangularMMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    LowerTriangularMMatrix(A, Val(n))
end
@generated function LowerTriangularMMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = lower_triangle_quote(N, T)
    push!(q.args, :(LowerTriangularMMatrix{$T,$N,$N2}( Ref( @ntuple $N2 A ))))
    q
end
UpperTriangularMatrix(data::NTuple{N2,T}) where {T,N2} = UpperTriangularMatrix(data,ValI(Val{N2}()))
UpperTriangularMMatrix(data::NTuple{N2,T}) where {T,N2} = UpperTriangularMMatrix(data,ValI(Val{N2}()))
UpperTriangularMatrix(data::SVector{N2,T}) where {T,N2} = UpperTriangularMatrix(data.data,ValI(Val{N2}()))
UpperTriangularMMatrix(data::SVector{N2,T}) where {T,N2} = UpperTriangularMMatrix(data.data,ValI(Val{N2}()))
UpperTriangularMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = UpperTriangularMatrix{T,N,N2}(data)
UpperTriangularMMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = UpperTriangularMMatrix{T,N,N2}(Ref(data))
UpperTriangularMatrix(data::UpperTriangularMMatrix{T,N,N2}) where {T,N,N2} = UpperTriangularMatrix{T,N,N2}(data.data[])
UpperTriangularMMatrix(data::UpperTriangularMatrix{T,N,N2}) where {T,N,N2} = UpperTriangularMMatrix{T,N,N2}(Ref(data.data))

LowerTriangularMatrix(data::NTuple{N2,T}) where {T,N2} = LowerTriangularMatrix(data,ValI(Val{N2}()))
LowerTriangularMMatrix(data::NTuple{N2,T}) where {T,N2} = LowerTriangularMMatrix(data,ValI(Val{N2}()))
LowerTriangularMatrix(data::SVector{N2,T}) where {T,N2} = LowerTriangularMatrix(data.data,ValI(Val{N2}()))
LowerTriangularMMatrix(data::SVector{N2,T}) where {T,N2} = LowerTriangularMMatrix(data.data,ValI(Val{N2}()))
LowerTriangularMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = LowerTriangularMatrix{T,N,N2}(data)
LowerTriangularMMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = LowerTriangularMMatrix{T,N,N2}(Ref(data))
LowerTriangularMatrix(data::LowerTriangularMMatrix{T,N,N2}) where {T,N,N2} = LowerTriangularMatrix{T,N,N2}(data.data[])
LowerTriangularMMatrix(data::LowerTriangularMatrix{T,N,N2}) where {T,N,N2} = LowerTriangularMMatrix{T,N,N2}(Ref(data.data))

@inline Base.getindex(A::UpperTriangularMatrix, i::Integer) = A.data[i]
@inline Base.getindex(A::LowerTriangularMatrix, i::Integer) = A.data[i]
@inline function Base.getindex(A::UpperTriangularMMatrix{T}, i::Integer) where T
    unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), i)
end
@inline function Base.getindex(A::LowerTriangularMMatrix{T}, i::Integer) where T
    unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), i)
end
@inline function Base.getindex(A::UpperTriangularMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i > j
            return zero(T)
        end
    end
    @inbounds out = A.data[small_triangle(j)+i]
    out
end
@inline function Base.getindex(A::LowerTriangularMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i < j
            return zero(T)
        end
    end
    @inbounds out = A.data[small_triangle(i)+j]
    out
end
@inline function Base.getindex(A::UpperTriangularMMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i > j
            return zero(T)
        end
    end
    # @inbounds out = A.data[][small_triangle(j)+i]
    # out
    unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), small_triangle(j)+i)
end
@inline function Base.getindex(A::LowerTriangularMMatrix{T,N,N2}, i::Integer, j::Integer) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i < j
            return zero(T)
        end
    end
    # @inbounds out = A.data[][small_triangle(i)+j]
    # out
    unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), small_triangle(i)+j)
end
@inline function Base.setindex!(A::TriangularMatrix{T,N,N2}, val, i::Integer) where  {T,N,N2}
    @boundscheck i > N2 && throw(BoundsError())
    unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), convert(T, val), i)
    return val
end
@inline function Base.setindex!(A::UpperTriangularMMatrix{T,N,N2}, val, i::Integer, j::Integer) where  {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i > j
            throw("Cannot set index within lower triangular of UpperTriangularMatrix.")
        end
    end
    unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), convert(T, val), small_triangle(j)+i)
    return val
end
@inline function Base.setindex!(A::LowerTriangularMMatrix{T,N,N2}, val, i::Integer, j::Integer) where  {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
        if i < j
            throw("Cannot set index within upper triangular of LowerTriangularMatrix.")
        end
    end
    unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), convert(T, val), small_triangle(i)+j)
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
