module TriangularMatrices

using Compat, Compat.LinearAlgebra
const LinearAlgebra = Compat.LinearAlgebra
using StaticArrays, Base.Cartesian

const Sized{T,N} = Union{SVector{T,N}, NTuple{N,T}}

export  SymmetricMatrix,
        LowerTriangularMatrix,
        UpperTriangularMatrix,
        revchol

ltriangle(x)::Int = (x-1)*x ÷ 2
btriangle(x)::Int = (x+1)*x ÷ 2
itriangle(x)::Int = (Int(sqrt(1+8x))-1) ÷ 2
@generated ValL(::Val{x}) where x = Val{ltriangle(x)}()
@generated ValB(::Val{x}) where x = Val{btriangle(x)}()
@generated ValI(::Val{x}) where x = Val{itriangle(x)}()

include("triangular.jl")
include("symmetric.jl")
include("meta.jl")
include("cholesky.jl")
include("inverse.jl")
include("arithmetic.jl")

end # module
