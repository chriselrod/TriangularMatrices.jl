using TriangularMatrices
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
using StaticArrays

x = @SMatrix randn(5,3)
