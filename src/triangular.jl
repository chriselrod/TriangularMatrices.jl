"""
Upper triangular and lower triangular matrices are the same views of the underlying data.
Upper triangular matrices are column major, and lower triangular row major.

As an implementation detail, indexing with @inbounds into the 0-triangle (ie, lower triangle of an Upper Triangular matrix) is undefined. This is because that check creates a branch.
"""
abstract type TriangularMatrix{T, N, N2} <: AbstractMatrix{T} end

struct UpperTriangularMatrix{T, N, N2} <: TriangularMatrix{T, N, N2}
    data::NTuple{N2,T}
end
struct LowerTriangularMatrix{T, N, N2} <: TriangularMatrix{T, N, N2}
    data::NTuple{N2,T}
end
Base.size(::TriangularMatrix{T,N}) where {T,N} = (N,N)


function UpperTriangularMatrix(x::AbstractMatrix{T}) where T
    N = size(x,1)
    @assert N == size(x,2)
    N2 = btriangle(N)
    y = Vector{T}(undef, N2)
    ind = 0
    for i ∈ 1:N, j ∈ 1:i
        ind += 1
        y[ind] = x[j,i]
    end
    UpperTriangularMatrix{T,N,N2}(ntuple(i -> y[i], N2))
end

function Base.getindex(x::UpperTriangularMatrix{T,N,N2}, i::Int, j::Int) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
        if i > j
            return zero(T)
        end
    end
    @inbounds out = x.data[ltriangle(j)+i]
    out
end
function Base.getindex(x::LowerTriangularMatrix{T,N,N2}, i::Int, j::Int) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
        if i < j
            return zero(T)
        end
    end
    @inbounds out = x.data[ltriangle(i)+j]
    out
end
LinearAlgebra.transpose(x::UpperTriangularMatrix{T,N,N2}) where {T,N,N2} = LowerTriangularMatrix{T,N,N2}(x.data)
LinearAlgebra.transpose(x::LowerTriangularMatrix{T,N,N2}) where {T,N,N2} = UpperTriangularMatrix{T,N,N2}(x.data)
@static if VERSION > v"0.6.9"
    LinearAlgebra.adjoint(x::UpperTriangularMatrix{T,N,N2}) where {T<:Real,N,N2} = LowerTriangularMatrix{T,N,N2}(x.data)
    LinearAlgebra.adjoint(x::LowerTriangularMatrix{T,N,N2}) where {T<:Real,N,N2} = UpperTriangularMatrix{T,N,N2}(x.data)
end
