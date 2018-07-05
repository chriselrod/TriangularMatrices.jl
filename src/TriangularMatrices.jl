module TriangularMatrices

using   LinearAlgebra,
        Base.Cartesian,
        MacroTools,
        StaticArrays # Drop, eventually?

const Sized{N,T} = Union{SVector{N,T}, NTuple{N,T}}

import LinearAlgebra: Adjoint

export  SymmetricMatrix,
        StaticSymmetricMatrix,
        LowerTriangularMatrix,
        StaticLowerTriangularMatrix,
        UpperTriangularMatrix,
        StaticUpperTriangularMatrix,
        RecursiveMatrix,
        StaticRecursiveMatrix,
        # cholesky,
        cholesky!,
        reverse_cholesky,
        reverse_cholesky!,
        xxt,
        xtx,
        # inv!,
        randmat,
        srandmat,
        RecursiveVector,
        choldet!,
        invdet!





# cutoff must be even
# The indexing recursion terminates at the cutoff
const halfcutoff = 4
const cutoff = 2halfcutoff
const cutoff2 = abs2(cutoff)
# const halfcut = cld(cutoff, 2)

small_triangle(x)::Int = (x-1)*x ÷ 2
big_triangle(x)::Int = (x+1)*x ÷ 2
inv_triangle(x)::Int = (Int(sqrt(1+8x))-1) ÷ 2
@generated ValSmallT(::Val{x}) where x = Val{small_triangle(x)}()
@generated ValBigT(::Val{x}) where x = Val{big_triangle(x)}()
@generated ValInvT(::Val{x}) where x = Val{inv_triangle(x)}()

include("meta.jl")
include("recursive_indexing.jl")
include("matrix_types/triangular.jl")
include("matrix_types/symmetric.jl")
include("matrix_types/recursive_matrix.jl")
include("matrix_types/rvector.jl")
include("decomp_and_inversions/cholesky.jl")
include("decomp_and_inversions/inverse.jl")
include("arithmetic/arithmetic.jl")
include("arithmetic/gemm.jl")
include("arithmetic/symv.jl")

end # module
