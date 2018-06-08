

@static if VERSION < v"0.7-"
    struct Adjoint{T,A} <: AbstractMatrix{T}
        parent::A
    end
else
    import LinearAlgebra: Adjoint
end



abstract type AbstractRecursiveMatrix{T,M,N,L} <: StaticArrays.StaticArray{Tuple{M,N}, T, 2} end
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

randmat(m,n) = RecursiveMatrix{Float64,m,n,m*n}(ntuple(i -> randn(), Val(m*n)))
randmat(n) = randmat(n,n)

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

"""
The methods in this library assume that operations on the pointer matrices do not alias.
It is discouraged to use these externally at all. If you do, make sure you don't violate these assumptions.

For this reason, few methods are actually implemented on this type. It is used internally via direct llvmcalls.
"""
struct PointerRecursiveMatrix{T,M,N,L} <: MutableRecursiveMatrix{T,M,N,L}
    data::Ptr{T}
    PointerRecursiveMatrix(ptr::Ptr{T}, ::Val{M}, ::Val{N}, ::Val{L}) where {T,M,N,L} = new{T,M,N,L}(ptr)
end
function PointerRecursiveMatrix(parent::RecursiveMatrix{T,MP,NP,LP}, ::Val{M}, ::Val{N}, offset=0) where {T,M,N,MP,NP,LP}
    PointerRecursiveMatrix(
        Base.unsafe_convert(Ptr{T}, pointer_from_objref(parent) + offset * sizeof(T)),
        Val{M}(), Val{N}(), ValProd(Val{M}(), Val{N}())
    )
end



const RecursiveMatrixOrTranpose{T,M,N,L} = Union{
    RecursiveMatrix{T,M,N,L},
    Adjoint{T,RecursiveMatrix{T,N,M,L}}
}
const RecursivePointerMatrixOrTranpose{T,M,N,L} = Union{
    PointerRecursiveMatrix{T,M,N,L},
    Adjoint{T,PointerRecursiveMatrix{T,N,M,L}}
}
const RecurseOrTranpose{T,M,N,L} = Union{
    RecursiveMatrix{T,M,N,L},
    Adjoint{T,RecursiveMatrix{T,N,M,L}},
    PointerRecursiveMatrix{T,M,N,L},
    Adjoint{T,PointerRecursiveMatrix{T,N,M,L}}
}
function dereftype(::Type{PointerRecursiveMatrix{T,M,N,L}}) where {T,M,N,L}
    RecursiveMatrix{T,M,N,L}
end
function dereftype(::Type{Adjoint{T,PointerRecursiveMatrix{T,N,M,L}}}) where {T,M,N,L}
    Adjoint{T,RecursiveMatrix{T,N,M,L}}
end


# abstract type transposition end
# struct t <: transposition end
# struct n <: transposition end
# tranpose(::t) = n()
# tranpose(::n) = t()
istransposed(::AbstractRecursiveMatrix) = false #n()
istransposed(::Type{<:AbstractRecursiveMatrix}) = false #n()
istransposed(::Adjoint{T,<:AbstractRecursiveMatrix{T}}) where T = true #t()
istransposed(::Type{<:Adjoint{T,<:AbstractRecursiveMatrix{T}}}) where T = true #t()

sub2ind(t, dims, i, j) = t ? j + (i-1)*dims[2] : i + (j-1)*dims[1]


Base.size(::RecurseOrTranpose{T,M,N}) where {T,M,N} = (M,N)
# @static if VERSION < v"0.7-"
@inline Base.getindex(A::Adjoint{T,<:MutableRecursiveMatrix{T}}, i) where T = A.parent[i]
@inline Base.setindex!(A::Adjoint{T,<:MutableRecursiveMatrix{T}}, v, i) where T = (A.parent[i] = v)
# end


# function forloopmul!(C, A, B)
#     @inbounds for i = 1:size(B,2)
#         for j = 1:size(A,1)
#             t = A[j,1] * B[1,j]
#             for k = 1:size(A,2)
#                 t += A[j,k] * B[k,j]
#             end
#             C[j,i] = t
#         end
#     end
# end

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
# @generated function Base.getindex(A::RecursiveMatrixOrTranpose{T,M,N,L}, i::Int, ::Type{K}) where {T,M,N,L,K}
#     bounds_error_string = "($M, $N) array at index (\$i,\$j)."
#     quote #Add recursion.
#         Base.@_inline_meta
#         # @boundscheck begin
#         #     $L < i && throw(BoundsError($bounds_error_string))
#         # end
#         unsafe_load(Base.unsafe_convert(Ptr{$K}, pointer_from_data(A.data) + i*$(sizeof(T)) ), 1 )
#         # A.data[][i]
#     end
# end


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

# rmatA8 = RecursiveMatrix{Float64,8,8,abs2(8)}(Ref(ntuple(i -> randn(), Val(abs2(8)))));
# rmatB8 = RecursiveMatrix{Float64,8,8,abs2(8)}(Ref(ntuple(i -> randn(), Val(abs2(8)))));
# rmatC8 = RecursiveMatrix{Float64,8,8,abs2(8)}(Ref{NTuple{64,Float64}}());
# rmatA12 = RecursiveMatrix{Float64,12,12,abs2(12)}(Ref(ntuple(i -> randn(), Val(abs2(12)))));
# rmatB12 = RecursiveMatrix{Float64,12,12,abs2(12)}(Ref(ntuple(i -> randn(), Val(abs2(12)))));
# rmatC12 = RecursiveMatrix{Float64,12,12,abs2(12)}(Ref{NTuple{144,Float64}}());
# rmatA16 = RecursiveMatrix{Float64,16,16,abs2(16)}(Ref(ntuple(i -> randn(), Val(abs2(16)))));
# rmatB16 = RecursiveMatrix{Float64,16,16,abs2(16)}(Ref(ntuple(i -> randn(), Val(abs2(16)))));
# rmatC16 = RecursiveMatrix{Float64,16,16,abs2(16)}(Ref{NTuple{256,Float64}}());