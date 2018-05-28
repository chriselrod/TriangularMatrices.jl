abstract type AbstractSymmetricMatrix{T,N,N2} <: AbstractMatrix{T} end

"""
Actually Hermitian, because adjoint returns itself.
"""
struct SymmetricMatrix{T,N,N2} <: AbstractSymmetricMatrix{T,N,N2}
    data::NTuple{N2,T}
end
struct SymmetricMMatrix{T,N,N2} <: AbstractSymmetricMatrix{T,N,N2}
    data::Base.RefValue{NTuple{N2,T}}
    function SymmetricMMatrix{T, N, N2}(d::Base.RefValue{NTuple{N2,T}}) where {T,N,N2}
        isbits(T) || error("Can only construct mutable isbits matrices.")
        new{T,N,N2}(d)
    end
end
Base.size(::AbstractSymmetricMatrix{T,N}) where {T,N} = (N,N)


function SymmetricMatrix(A::AbstractMatrix)
    m, n = size(A)
    @assert m == n
    SymmetricMatrix(A, Val(n))
end
@generated function SymmetricMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = upper_triangle_quote(N, T)
    push!(q.args, :(SymmetricMatrix{$T,$N,$N2}( @ntuple $N2 A )))
    q
end
@generated function SymmetricMMatrix(A::AbstractMatrix{T}, ::Val{N}) where {T,N}
    q, N2 = upper_triangle_quote(N, T)
    push!(q.args, :(SymmetricMMatrix{$T,$N,$N2}( Ref( @ntuple $N2 A ))))
    q
end
function SymmetricMatrix(data::AbstractVector)
    N2 = length(data)
    N = itriangle(N2)
    SymmetricMatrix(data, Val(N), Val(N2))
end
function SymmetricMatrix(data::AbstractVector{T}, ::Val{N}, ::Val{N2}) where {T,N,N2}
    SymmetricMatrix{T,N,N2}(ntuple(i -> data[i], N2))
end
SymmetricMatrix(data::NTuple{N2,T}) where {T,N2} = SymmetricMatrix(data,ValI(Val{N2}()))
SymmetricMMatrix(data::NTuple{N2,T}) where {T,N2} = SymmetricMMatrix(data,ValI(Val{N2}()))
SymmetricMatrix(data::SVector{N2,T}) where {T,N2} = SymmetricMatrix(data.data,ValI(Val{N2}()))
SymmetricMMatrix(data::SVector{N2,T}) where {T,N2} = SymmetricMMatrix(data.data,ValI(Val{N2}()))
SymmetricMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = SymmetricMatrix{T,N,N2}(data)
SymmetricMMatrix(data::NTuple{N2,T}, ::Val{N}) where {T,N,N2} = SymmetricMMatrix{T,N,N2}(Ref(data))
SymmetricMatrix(data::SymmetricMMatrix{T,N,N2}) where {T,N,N2} = SymmetricMatrix{T,N,N2}(data.data[])
SymmetricMMatrix(data::SymmetricMatrix{T,N,N2}) where {T,N,N2} = SymmetricMMatrix{T,N,N2}(Ref(data.data))



function Base.getindex(x::SymmetricMatrix{T,N,N2}, i::Int, j::Int) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    @inbounds out = x.data[ltriangle(j)+i]
    out
end
function Base.getindex(x::SymmetricMMatrix{T,N,N2}, i::Int, j::Int) where {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw("BoundsError: index ($i, $j) is out of bounds.")
    end
    i, j = minmax(i, j)
    @inbounds out = x.data[][ltriangle(j)+i]
    out
end
@inline function Base.setindex!(A::SymmetricMMatrix{T,N,N2}, val, i::Integer, j::Integer) where  {T,N,N2}
    @boundscheck begin
        max(i,j) > N && throw(BoundsError())
    end
    i, j = minmax(i, j)
    unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(A.data)), convert(T, val), ltriangle(j)+i)
    return val
end

LinearAlgebra.transpose(x::SymmetricMatrix) = x
@static if VERSION >= v"0.7.0-DEV"
    LinearAlgebra.adjoint(x::SymmetricMatrix) = x
end
