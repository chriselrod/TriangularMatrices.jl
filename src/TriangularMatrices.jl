module TriangularMatrices

using Compat, Compat.LinearAlgebra
const LinearAlgebra = Compat.LinearAlgebra
using StaticArrays, Base.Cartesian

const Sized{N,T} = Union{SVector{N,T}, NTuple{N,T}}

export  SymmetricMatrix,
        SymmetricMMatrix,
        LowerTriangularMatrix,
        LowerTriangularMMatrix,
        UpperTriangularMatrix,
        UpperTriangularMMatrix,
        cholesky,
        cholesky!,
        reverse_cholesky,
        reverse_cholesky!,
        xxt,
        xtx,
        inv!



# The recursions turminate below the cutoff.
# Halfcut is the smallest possible value we need to have a kernel for.
# If cutoff is even, cutoff + 1 would dispatch to kernels for halfcut+1 and halfcut.
# If cutoff is odd, cutoff + 1 would dispatch to kernels for halfcut.
const cutoff = 24
const halfcut = cld(cutoff, 2)

small_triangle(x)::Int = (x-1)*x ÷ 2
big_triangle(x)::Int = (x+1)*x ÷ 2
inv_triangle(x)::Int = (Int(sqrt(1+8x))-1) ÷ 2
@generated ValSmallT(::Val{x}) where x = Val{small_triangle(x)}()
@generated ValBigT(::Val{x}) where x = Val{big_triangle(x)}()
@generated ValInvT(::Val{x}) where x = Val{inv_triangle(x)}()

include("triangular.jl")
include("symmetric.jl")
include("meta.jl")
include("cholesky.jl")
include("inverse.jl")
include("arithmetic.jl")

end # module
