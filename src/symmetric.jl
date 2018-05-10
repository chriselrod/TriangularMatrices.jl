

struct SymmetricMatrix{T,N,N2} <: AbstractMatrix{T}
    data::NTuple{N,T}
end
Base.size(::SymmetricMatrix{T,N}) where {T,N} = (N,N)


function SymmetricMatrix(x::AbstractMatrix{T}) where T
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

function Base.getindex(x::SymmetricMatrix{T,N,N2}, i::Int, j::Int) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    @inbounds out = x.data[ltriangle(j)+i]
    out
end

LinearAlgebra.transpose(x::SymmetricMatrix) = x
LinearAlgebra.adjoint(x::SymmetricMatrix{T}) where {T<:Real} = x

ltriangle(x) = (x-1)*x ÷ 2
btriangle(x) = (x+1)*x ÷ 2