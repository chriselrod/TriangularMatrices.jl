
"""
Actually Hermitian, because adjoint returns itself.
"""
struct SymmetricMatrix{T,N,N2} <: AbstractMatrix{T}
    data::NTuple{N2,T}
end
Base.size(::SymmetricMatrix{T,N}) where {T,N} = (N,N)


function SymmetricMatrix(x::AbstractMatrix{T}) where T
    N = size(x,1)
    @assert N == size(x,2)
    N2 = btriangle(N)
    data = Vector{T}(undef, N2)
    ind = 0
    for i ∈ 1:N, j ∈ 1:i
        ind += 1
        data[ind] = x[j,i]
    end
    SymmetricMatrix{T,N,N2}(ntuple(i -> data[i], N2))
end
function SymmetricMatrix(data::AbstractVector{T}) where T
    N2 = length(data)
    N = itriangle(N2)
    SymmetricMatrix{T,N,N2}(ntuple(i -> data[i], N2))
end
function SymmetricMatrix(data::Sized{T,N2}) where {T,N2}
    SymmetricMatrix(data, ValI(Val{N2}()))
end
function SymmetricMatrix(data::SVector{T,N2}, ::Val{N}) where {T,N,N2}
    SymmetricMatrix{T,N,N2}(ntuple(i -> data[i], Val{N2}()))
end
function SymmetricMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2}
    SymmetricMatrix{T,N,N2}(ntuple(i -> data[i], Val{N2}()))
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
@static if VERSION > v"0.6.9"
    LinearAlgebra.adjoint(x::SymmetricMatrix) = x
end
