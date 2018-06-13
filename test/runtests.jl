using TriangularMatrices
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end


# write your own tests here
using StaticArrays
using Compat, Compat.LinearAlgebra


x = @SMatrix randn(5,3)
y = @SMatrix randn(3,5)
SAS = x' * x
SAS2 = y * y'

TMS = xtx(x)
TMS2 = xxt(y)

# @show TMS
# @show SAS
# @show TMS2
# @show SAS2

@test all(TMS .== SAS)
@test all(TMS2 .== SAS2)

function check_approx_equality(A, R)
    @assert size(A) == size(R)
    m, n = size(A)
    all_true = true
    for i = 1:n, j = 1:m
        all_true = A[j,i] ≈ R[j,i]
        all_true || ( @show((j,i)); break )
    end
    all_true
end



# @time @testset for M = 1:128, N = round(Int, M/√2, RoundUp):round(Int, M*√2, RoundDown), P = max(round(Int, M/√2, RoundUp), round(Int, N/√2, RoundUp)):min(round(Int, M*√2, RoundDown),round(Int, N*√2, RoundDown))

@time @testset for M = 20:128, N = round(Int, M/√2, RoundUp):round(Int, M*√2, RoundDown), P = max(round(Int, M/√2, RoundUp), round(Int, N/√2, RoundUp)):min(round(Int, M*√2, RoundDown),round(Int, N*√2, RoundDown))
        d_r = randmat(M,P)
        a_r = randmat(M,N)
        x_r = randmat(N,P)
        TriangularMatrices.mul!(d_r, a_r, x_r)
        d_m = Matrix(d_r)
        if VERSION < v"0.7-"
            A_mul_B!(d_m, Matrix(a_r), Matrix(x_r))
        else
            mul!(d_m, Matrix(a_r), Matrix(x_r))
        end
        @test check_approx_equality(d_m, d_r)
    end
end


