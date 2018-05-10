module TriangularMatrices

using Compat, Compat.LinearAlgebra
const LinearAlgebra = Compat.LinearAlgebra
using StaticArrays, Base.Cartesian

include("triangular.jl")
include("symmetric.jl")
include("cholesky.jl")
include("inverse.jl")

end # module
