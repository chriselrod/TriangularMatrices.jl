# TriangularMatrices

[![Build Status](https://travis-ci.org/chriselrod/TriangularMatrices.jl.svg?branch=master)](https://travis-ci.org/chriselrod/TriangularMatrices.jl)

[![Coverage Status](https://coveralls.io/repos/chriselrod/TriangularMatrices.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chriselrod/TriangularMatrices.jl?branch=master)

[![codecov.io](http://codecov.io/github/chriselrod/TriangularMatrices.jl/coverage.svg?branch=master)](http://codecov.io/github/chriselrod/TriangularMatrices.jl?branch=master)




---


This package is probably totally broken right now.
It has been ages since I tried `using` it.

But here is minimal test code you can run. I am however adding comments on where that code came from.
You may be especially interested in looking at src/arithmetic/gemm.jl, and looking at all the commented our methods I've tried for generating kernels.
```julia

const cutoff = 16

###
### from meta.jl
###

extract_symbol(A, i) = :($A[$i])
id_symbol(A, i) = Symbol(A, :_, i)


function create_quote(fast::Bool = true)
    if fast
        q = quote @fastmath @inbounds begin end end
        @static if VERSION > v"0.7-"
            qa = q.args[2].args[3].args[3].args
        else
            qa = q.args[2].args[2].args[2].args
        end
    else
        q = quote @inbounds begin end end
        @static if VERSION > v"0.7-"
            qa = q.args[2].args[3].args
        else
            qa = q.args[2].args[2].args
        end
    end
    q, qa
end

function create_quote(A_i, A)
    q, qa = create_quote()
    push!(qa, :($A_i = $A[1]))
    q, qa
end

function extract_linear!(qa, N, prefix = :B)
    prefix_ = Symbol(prefix, :_)
    for i ∈ 1:N
        push!(qa, :($(Symbol(prefix_, i)) = $(prefix)[$i]))
    end
    qa
end
function insert_linear!(qa, N, prefix = :B)
    prefix_ = Symbol(prefix, :_)
    for i ∈ 1:N
        push!(qa, :( $(prefix)[$i] = $(Symbol(prefix_, i))) )
    end
    qa
end

@generated ValProd(::Val{M}, ::Val{N}) where {M,N} = Val{M*N}()
@generated function ValDiv(::Val{L}, ::Val{M}) where {L,M}
    N, r = divrem(L, M)
    @assert r == 0
    Val{N}()
end

###
### Recursive indices, not actually needed for the example because we're only looking at the kernels.
### From src/recursive_indexing.jl
###


function dense_sub2ind_quote(M, N, offset = 0, joffset = 1)
    minind, maxind = minmax(M, N)
    maxind < cutoff && return :( (j-$joffset ) * $M + i + $offset )
    
    if √2*minind > maxind && minind > cutoff # we divide the matrix into four blocks.
        Mh, Mr = divrem(M, 2)
        Nh, Nr = divrem(N, 2)
        Mh_p_Mr = Mh + Mr
        Nh_p_Nr = Nh + Nr
        return :(if i <= $Mh_p_Mr && j <= $Nh_p_Nr # enter block 1
                $(dense_sub2ind_quote(Mh_p_Mr, Nh_p_Nr, offset, joffset))
            elseif j <= $Nh_p_Nr # enter block 2
                $(dense_sub2ind_quote(Mh, Nh_p_Nr, offset + Mh_p_Mr * (Nh_p_Nr-1), joffset))
            elseif i <= $Mh_p_Mr # enter block 3
                $(dense_sub2ind_quote(Mh_p_Mr, Nh, offset + M * Nh_p_Nr, joffset + Nh_p_Nr))
            else # enter block 4
                $(dense_sub2ind_quote(Mh, Nh, offset + M * Nh_p_Nr + Mh_p_Mr * (Nh-1), joffset + Nh_p_Nr))
            end)
    elseif M > N # we are splitting the matrix into two blocks stacked on top of one another.
        Mh, Mr = divrem(M, 2)
        Mh_p_Mr = Mh + Mr

        return :(if i < $Mh_p_Mr
                $(dense_sub2ind_quote(Mh_p_Mr, N, offset, joffset))
            else
                $(dense_sub2ind_quote(Mh, N, offset + Mh_p_Mr * (N-1), joffset))
            end)
    else # we are splitting the matrix into two blocks, side by side.
        Nh, Nr = divrem(N, 2)
        Nh_p_Nr = Nh + Nr

        return :(if j < $Nh_p_Nr
                $(dense_sub2ind_quote(M, Nh_p_Nr, offset, joffset))
            else
                $(dense_sub2ind_quote(M, Nh, offset + M * Nh_p_Nr, joffset + Nh_p_Nr))
            end)
    end
end

"""
Cartesian indexing isn't recomended, but it is convenient for printing.
The approach for sub2ind here incurs branches, which would disable SIMD.
Therefore, Cartesian indexing is strongly discouraged for hot loops.
Still, it is reasonably fast -- only a handful of ns.
"""
@generated dense_sub2ind(::Val{M}, ::Val{N}, i, j) where {M,N} = dense_sub2ind_quote(M, N)



###
### From src/recursive_matrix.jl
###
@static if VERSION < v"0.7-"
    struct Adjoint{T,A} <: AbstractMatrix{T}
        parent::A
    end
else
    import LinearAlgebra: Adjoint
end



abstract type AbstractRecursiveMatrix{T,M,N,L} <: AbstractArray{T,2} end# StaticArrays.StaticArray{Tuple{M,N}, T, 2} end
abstract type MutableRecursiveMatrix{T,M,N,L} <: AbstractRecursiveMatrix{T,M,N,L} end
mutable struct RecursiveMatrix{T,M,N,L} <: MutableRecursiveMatrix{T,M,N,L}
    data::NTuple{L,T}
    RecursiveMatrix{T,M,N,L}(data::NTuple{L,T}) where {T,M,N,L} = new(data)
    RecursiveMatrix{T,M,N,L}() where {T,M,N,L} = new()
end

function RecursiveMatrix(data::NTuple{L,T}, ::Val{M}, ::Val{N}) where {M,N,L,T}
    RecursiveMatrix{T,M,N,L}(data)
end
function RecursiveMatrix(::Type{T}, ::Val{M}, ::Val{N}, ::Val{L}) where {M,N,L,T}
    RecursiveMatrix{T,M,N,L}()
end
function RecursiveMatrix(::Type{T}, ::Val{M}, ::Val{N}) where {T,M,N}
    RecursiveMatrix(T,Val{M}(),Val{N}(),ValProd(Val{M}(),Val{L}()))
end
function RecursiveMatrix(data::NTuple{L,T}, ::Val{M}) where {M,L,T}
    RecursiveMatrix(data, Val{M}(), ValDiv(Val{L}(),Val{M}()))
end

"""
Creating a recursive matrix without compile time size information is not type stable.
"""
function RecursiveMatrix(data::AbstractMatrix)
    M, N = size(data)
    RecursiveMatrix(data, Val{M}(), Val{N}(), Val{M*N}())
end
function RecursiveMatrix(data::AbstractMatrix, ::Val{M}, ::Val{N}) where {M,N}
    RecursiveMatrix(data, Val{M}(), Val{N}(), ValProd(Val{M}(), Val{N}()))
end
function RecursiveMatrix(data::AbstractMatrix{T}, ::Val{M}, ::Val{N}, ::Val{L}) where {T,M,N,L}
    out = RecursiveMatrix(T, Val{M}(), Val{N}() ,Val{L}())
    # Not an efficient order of indexing into the recursive matrix.
    # Optimize later by doing it in branch-less chunks.
    for i ∈ 1:N, j ∈ 1:M
        out[j,i] = data[j,i] 
    end
    out
end

const RecursiveMatrixOrTranpose{T,M,N,L} = Union{
    RecursiveMatrix{T,M,N,L},
    Adjoint{T,RecursiveMatrix{T,N,M,L}}
}

istransposed(::AbstractRecursiveMatrix) = false #n()
istransposed(::Type{<:AbstractRecursiveMatrix}) = false #n()
istransposed(::Adjoint{T,<:AbstractRecursiveMatrix{T}}) where T = true #t()
istransposed(::Type{<:Adjoint{T,<:AbstractRecursiveMatrix{T}}}) where T = true #t()

sub2ind(t, dims, i, j) = t ? j + (i-1)*dims[2] : i + (j-1)*dims[1]

# Base.size(::RecurseOrTranpose{T,M,N}) where {T,M,N} = (M,N)
Base.size(::RecursiveMatrixOrTranpose{T,M,N}) where {T,M,N} = (M,N)


@inline point(A::RecursiveMatrix{T}) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
@inline point(A::Adjoint{T,<:RecursiveMatrix{T}}) where T = point(A.parent)

@generated function Base.getindex(A::RecursiveMatrixOrTranpose{T,M,N,L}, i::Int) where {T,M,N,L}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."

    quote #Add recursion.
        Base.@_inline_meta
        @boundscheck begin
            $L < i && throw(BoundsError($bounds_error_string))
        end
        A.data[i]
        # A.data[][i]
    end
end


@generated function Base.getindex(A::RecursiveMatrixOrTranpose{T,M,N,L}, i::Int, j::Int) where {T,M,N,L}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    if istransposed(A)
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), j, i))
    else
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), i, j))
    end
    quote #Add recursion.
        Base.@_inline_meta
        @boundscheck begin
            ($M < i || $N < j) && throw(BoundsError($bounds_error_string))
        end
        A.data[$linear_ind_expr]
    end
end


@generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val, i::Int) where {T,M,N,L}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    quote #Add recursion.
        Base.@_inline_meta
        @boundscheck begin
            $L < i && throw(BoundsError($bounds_error_string))
        end
        unsafe_store!(point(A), convert(T, val), i )
    end
end
# @generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val::K, i::Int, ::Type{K}) where {T,M,N,L,K}
#     bounds_error_string = "($M, $N) array at index (\$i,\$j)."
#     quote #Add recursion.
#         Base.@_inline_meta
#         # @boundscheck begin
#         #     $L < i && throw(BoundsError($bounds_error_string))
#         # end
#         unsafe_store!(Base.unsafe_convert(Ptr{$K}, pointer_from_data(A.data) + i*$(sizeof(T)) ), 1 )
#     end
# end

@generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val, i::Int, j::Int) where {T,M,N,L}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    if istransposed(A)
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), j, i))
    else
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), i, j))
    end
    quote #Add recursion.
        Base.@_inline_meta
        @boundscheck begin
            ($M < i || $N < j) && throw(BoundsError($bounds_error_string))
        end
        # unsafe_store!(point(A), convert(T, val), sub2ind(istransposed(A), ($M,$N), i, j) )

        unsafe_store!(point(A), convert(T, val), $linear_ind_expr )
    end
end



###
### from src/arithmetic/gemm.jl
###

"""
Use `q, qa = create_quote(); mul_quote!(qa, M, N, P); q` to get an idea of what the
generated Julia code looks like.
"""
function mul_quote!(qa, M, N, P, ta=false, A = :A, tb=false, B = :B, tc=false, C = :C,
                        extract = extract_symbol, insert = extract_symbol, eq = :(=))

    eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
    eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
    eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
    chunk = 4
    N4, Nr = divrem(N, chunk)
    for j = 1:P, i = 1:M
        C_ij = eC(i, j)
        push!(qa, Expr(eq, C_ij, :(+$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1:chunk]...) ) ))
        for k = 2:N4
            push!(qa, :($C_ij = $C_ij +$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+chunk*(k-1):chunk*k ]...) ) )
        end
        Nr > 0 && push!(qa, :($C_ij = $C_ij +$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+chunk*k:N ]...) ) )
    end
end

function mul_kernel(M, N, P, tA=false, tB=false, tC=false, eq = :(=), LA = M*N, LB = N*P, LC = M*P)
    q, qa = create_quote()
    push!(q.args, :(Base.@_inline_meta))
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LB, :B)
    mul_quote!(qa, M, N, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )
    insert_linear!(qa, LC, :C)
    # isa(dummy, Bool) || push!(q.args, :C) # worth considering?
    push!(q.args, :C)
    q
end

@generated function mul!(C::RecursiveMatrixOrTranpose{T,M,P,LC},
                        A::RecursiveMatrixOrTranpose{T,M,N,LA},
                        B::RecursiveMatrixOrTranpose{T,N,P,LB},
                        dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    mul_kernel(M, N, P, istransposed(A), istransposed(B), istransposed(C))
end

"""
Not type stable, but lazy.
"""
randsquare(n) = RecursiveMatrix{Float64,n,n,abs2(n)}(ntuple(i -> randn(), Val(abs2(n))))

a8 = randsquare(8);
b8 = randsquare(8);
c8 = randsquare(8);

mul!(c8, a8, b8)

using BenchmarkTools
@benchmark mul!($c8, $a8, $b8)

@code_llvm mul!(c8, a8, b8)
@code_native mul!(c8, a8, b8)

```


