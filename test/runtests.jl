using TriangularMatrices
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
using StaticArrays

x = @SMatrix randn(5,3)
y = @SMatrix randn(3,5)
SAS = x' * x
SAS2 = y * y'

TMS = xtx(x)
TMS2 = xxt(y)

@test all(TMS .== SAS)
@test all(TMS2 .== SAS2)