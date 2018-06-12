

@static if VERSION < v"0.7-"
    struct Adjoint{T,A} <: AbstractMatrix{T}
        parent::A
    end
else
    import LinearAlgebra: Adjoint
end



abstract type AbstractRecursiveMatrix{T,M,N,L} <: StaticArrays.StaticArray{Tuple{M,N}, T, 2} end
abstract type MutableRecursiveMatrix{T,M,N,L} <: AbstractRecursiveMatrix{T,M,N,L} end

struct StaticRecursiveMatrix{T,M,N,L} <: AbstractRecursiveMatrix{T,M,N,L}
    data::NTuple{L,T}
    StaticRecursiveMatrix{T,M,N,L}(data::NTuple{L,T}) where {T,M,N,L} = new(data)
end
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

function randmat(::Val{M},::Val{N},::Val{L}) where {M,N,L}
    out = RecursiveMatrix{Float64,M,N,L}()
    @inbounds for i = 1:L
        out[i] = randn()
    end
    out
end
randmat(m,n) = randmat(Val(m),Val(n),Val(m*n))
randmat(n) = randmat(n,n)
srandmat(m,n) = StaticRecursiveMatrix{Float64,m,n,m*n}(ntuple(i -> randn(), Val(m*n)))
srandmat(n) = srandmat(n,n)

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
    data::Ptr{StaticRecursiveMatrix{T,cutoff,cutoff,cutoff2}}
    # PointerRecursiveMatrix(ptr::Ptr{StaticRecursiveMatrix{T,cutoff,cutoff,cutoff2}}, ::Val{M}, ::Val{N}, ::Val{L}) where {T,M,N,L} = new{T,M,N,L}(ptr)
end
function PointerRecursiveMatrix(parent::RecursiveMatrix{T,MP,NP,LP}, ::Val{M}, ::Val{N}, offset=0) where {T,M,N,MP,NP,LP}
    PointerRecursiveMatrix(
        Base.unsafe_convert(Ptr{StaticRecursiveMatrix{T,cutoff,cutoff,cutoff2}}, pointer_from_objref(parent) + offset * sizeof(T)),
        Val{M}(), Val{N}(), ValProd(Val{M}(), Val{N}())
    )
end

# Needs to be made faster.
function Base.Matrix(A::RecursiveMatrix{T,M,N,L}) where {T,M,N,L}
    out = Matrix{T}(undef, M, N)
    @inbounds for i = 1:N, j = 1:M
        out[j,i] = A[j,i]
    end
    out
end

struct Block{m,n} end
struct HalfBlock{m,n,h} end
Base.@pure Block(m,n) = Block{m,n}()
Base.@pure HalfBlock(m,n,h) = HalfBlock{m,n,h}()
struct BlockIndex{T,M,N,L}
    i::UInt
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
    StaticRecursiveMatrix{T,M,N,L},
    Adjoint{T,StaticRecursiveMatrix{T,M,N,L}},
    PointerRecursiveMatrix{T,M,N,L},
    Adjoint{T,PointerRecursiveMatrix{T,N,M,L}}
}
const MutableRecurseOrTranpose{T,M,N,L} = Union{
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
vistransposed(::AbstractRecursiveMatrix) = Val{false}() #n()
vistransposed(::Type{<:AbstractRecursiveMatrix}) = Val{false}() #n()
vistransposed(::Adjoint{T,<:AbstractRecursiveMatrix{T}}) where T = Val{true}() #t()
vistransposed(::Type{<:Adjoint{T,<:AbstractRecursiveMatrix{T}}}) where T = Val{true}() #t()

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

@inline float_point(A::RecursiveMatrix{T}) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
@inline float_point(A::PointerRecursiveMatrix{T}) where T = Base.unsafe_convert(Ptr{T}, A.data)
@inline point(A::RecursiveMatrix) = pointer_from_objref(A)
@inline point(A::PointerRecursiveMatrix) = A.data
@inline convert_point(A::RecursiveMatrix, ::Type{T}) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
@inline convert_point(A::RecursiveMatrix{T}) where T = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{T,cutoff,cutoff,cutoff2}}, pointer_from_objref(A))
@inline convert_point(A::PointerRecursiveMatrix, ::Type{T}) where T = Base.unsafe_convert(Ptr{T}, A.data)
@inline convert_point(A::PointerRecursiveMatrix) = A.data
@inline point(A::Adjoint{T,<:MutableRecurseOrTranpose{T}}) where T = point(A.parent)

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
@inline function Base.getindex(A::StaticRecursiveMatrix{T,M,N,L}, i::Int) where {T,M,N,L}
    @boundscheck i > L && throw(BoundsError())
    @inbounds out = A.data[i]
    out
end
@inline function Base.getindex(A::StaticRecursiveMatrix{T,M,N,L}, i::Int, j::Int) where {T,M,N,L}
    @boundscheck ( i > M || j > N ) && throw(BoundsError())
    @inbounds out = A.data[ (j-1)*M + i ]
    out
end

Base.getindex(A::RecursiveMatrixOrTranpose{T,M,N}, i::Int, j::Int) where {T,M,N} = A.data[dense_sub2ind((M,N),i,j)]
# @generated function Base.getindex(A::RecursiveMatrixOrTranpose{T,M,N,L}, i::Int, j::Int) where {T,M,N,L}
#     bounds_error_string = "($M, $N) array at index (\$i,\$j)."
#     if istransposed(A)
#         linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), j, i))
#     else
#         linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), i, j))
#     end
#     quote 
#         # Base.@_inline_meta
#         @boundscheck ($M < i || $N < j) && throw(BoundsError($bounds_error_string))
#         A.data[$linear_ind_expr]
#     end
# end
# function getindex_block_quote(::Type{T}, full_m_blocks, mdim, m_r, full_n_blocks, ndim, n_r, m, n) where T
#     if m <= full_m_blocks && n <= full_n_blocks
#         Lout = mdim*ndim
#         if mdim == ndim == cutoff
#             q = quote
#                     Base.@_inline_meta
#                     Base.unsafe_load(convert_point(A) + $(sizeof(T)*( (m-1)*Lout + (n-1)*M*ndim) )   )
#                 end
#         else
#             q = quote
#                     Base.@_inline_meta
#                     Base.unsafe_load(convert_point(A,
#                         Ptr{StaticRecursiveMatrix{$T,$mdim,$ndim,$Lout}}) +
#                         $(sizeof(T)*( (m-1)*Lout + (n-1)*M*ndim) )   )
#             end
#         end

#     elseif m <= full_m_blocks # we're on the far right
#         Lout = mdim * n_r
#         q = quote
#                 Base.@_inline_meta
#                 Base.unsafe_load(convert_point( A,
#                     Ptr{StaticRecursiveMatrix{$T,$mdim,$n_r,$Lout}}) +
#                     $(sizeof(T)*( (m-1)*Lout + (n-1)*M*ndim) )  )
#             end

#     elseif n <= full_n_blocks
#         Lout = m_r * ndim
#         q = quote
#                 Base.@_inline_meta
#                 Base.unsafe_load(convert_point( A,
#                     Ptr{StaticRecursiveMatrix{$T,$m_r,$ndim,$Lout}}) +
#                     $(sizeof(T)*( (m-1)*ndim*mdim + (n-1)*M*ndim) ) )
#             end

#     else
#         Lout = m_r*n_r
#         q = quote
#                 Base.@_inline_meta
#                 Base.unsafe_load(convert_point( A,
#                     Ptr{StaticRecursiveMatrix{$T,$m_r,$n_r,$Lout}}) +
#                     $(sizeof(T)*( L - Lout ))  )
#             end

#     end    
#     q
# end

# @generated function Base.getindex(A::MutableRecursiveMatrix{T,M,N,L}, ::Block{m,n}) where {T,M,N,L,m,n}
#     if max(M,N) <= 3cutoff # We are returning a StaticRecursiveMatrix
#         full_m_blocks, incomplete_m_block = divrem(M, cutoff)
#         full_n_blocks, incomplete_n_block = divrem(N, cutoff)

#         q = getindex_block_quote(T, full_m_blocks, cutoff, incomplete_m_block, full_n_blocks, cutoff, incomplete_n_block, m, n)

#     else # We are returning a pointer matrix.

#         if reasonably_square(minmax(M,N)...) # split both dims
#             M1, M2 = split_dim(M)
#             N1, N2 = split_dim(N)

#             Mout = (M1, M2)[m]
#             Nout = (N1, N2)[n]
#             if n == 1
#                 offset = sizeof(T) * (m-1)*M1*N1
#             else # n == 2
#                 offset = sizeof(T) * ( (m-1)*M1*N2 + M*N1 )
#             end
#             q = :( PointerRecursiveMatrix{$T,$Mout,$Nout,$(Mout*Nout)}( convert_point(A) + $offset ) )

#         elseif M > N
#             M1, M2 = split_dim(M)

#             Mout = (M1, M2)[m]
#             offset = sizeof(T) * (m-1)*M1*N
#             q = :( PointerRecursiveMatrix{$T,$Mout,$N,$(Mout*N)}( convert_point(A) + $offset ) )

#         else # M < N
#             N1, N2 = split_dim(N)

#             Nout = (N1, N2)[n]
#             offset = sizeof(T) * (n-1)*M*N1
#             q = :( PointerRecursiveMatrix{$T,$M,$Nout,$(M*Nout)}( convert_point(A) + $offset ) )

#         end

#     end
#     q
# end

# @generated function Base.getindex(A::MutableRecursiveMatrix{T,M,N,L}, ::HalfBlock{m,n}) where {T,M,N,L,m,n}
#     full_m_blocks, incomplete_m_block = divrem(M, cutoff)
#     full_n_blocks, incomplete_n_block = divrem(N, halfcutoff)
#     getindex_block_quote(T, full_m_blocks, cutoff, incomplete_m_block, full_n_blocks, halfcutoff, incomplete_n_block, m, n)
# end

@generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val::T, i::Int) where {T,M,N,L}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    quote #Add recursion.
        Base.@_inline_meta
        @boundscheck $L < i && throw(BoundsError($bounds_error_string))
        unsafe_store!(float_point(A), val, i )
    end
end
@generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val::I, i::Int) where {T,M,N,L,I<:Integer}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    quote #Add recursion.
        # Base.@_inline_meta
        @boundscheck begin
            $L < i && throw(BoundsError($bounds_error_string))
        end
        unsafe_store!(float_point(A), convert(T,val), i )
    end
end
@generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val::V, i::Int) where {T,M,N,L,V}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    quote #Add recursion.
        Base.@_inline_meta
        @boundscheck begin
            $L < i && throw(BoundsError($bounds_error_string))
        end
        unsafe_store!(Base.unsafe_convert(Ptr{V}, pointer_from_objref(A) + sizeof(T)*(i-1) ), val)#, 1 )
    end
end

@generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val::T, i::Int, j::Int) where {T,M,N,L}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    if istransposed(A)
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), j, i))
    else
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), i, j))
    end
    quote #Add recursion.
        # Base.@_inline_meta
        @boundscheck begin
            ($M < i || $N < j) && throw(BoundsError($bounds_error_string))
        end
        # unsafe_store!(point(A), convert(T, val), sub2ind(istransposed(A), ($M,$N), i, j) )

        unsafe_store!(point(A), val, $linear_ind_expr )
    end
end


@generated function Base.setindex!(A::RecursiveMatrixOrTranpose{T,M,N,L}, val::V, i::Int, j::Int) where {T,M,N,L,V}
    bounds_error_string = "($M, $N) array at index (\$i,\$j)."
    if istransposed(A)
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), j, i))
    else
        linear_ind_expr = :(dense_sub2ind(Val{$M}(), Val{$N}(), i, j))
    end
    quote #Add recursion.
        # Base.@_inline_meta
        @boundscheck begin
            ($M < i || $N < j) && throw(BoundsError($bounds_error_string))
        end
        # unsafe_store!(point(A), convert(T, val), sub2ind(istransposed(A), ($M,$N), i, j) )
        unsafe_store!(Base.unsafe_convert(Ptr{V}, pointer_from_objref(A) + sizeof(T)*(($linear_ind_expr)-1) ), val)#, 1 )
    end
end

# """
# The second block indicates that the number of dimensions we are splitting things into.
# """
# @generated function BlockIndex(::Type{MutableRecursiveMatrix{T,M,N,L}}, ::Block{m,n}, ::Block{MC,NC}) where {T,M,N,L,m,n,MC,NC}
#     full_m_blocks, m_r = divrem(M, cutoff)
#     full_n_blocks, n_r = divrem(N, cutoff)

#     if m <= full_m_blocks && n <= full_n_blocks
#         return BlockIndex{T,cutoff,cutoff,cutoff2}(sizeof(T)*( (m-1)*cutoff2 + (n-1)*M*cutoff) )
#     elseif m <= full_m_blocks # we're on the far right
#         return BlockIndex{T,cutoff,n_r,cutoff*n_r}(sizeof(T)*( (m-1)*cutoff*n_r + (n-1)*M*cutoff))
#     elseif n <= full_n_blocks
#         return BlockIndex{T,m_r,cutoff,m_r*cutoff}(sizeof(T)*( (m-1)*cutoff2 + (n-1)*M*cutoff ) )
#     else
#         return BlockIndex{T,m_r,n_r,m_r*n_r}(sizeof(T)*( L - m_r*n_r ) )
#     end
# end
#
# This function shouldn't really be typed on m,n, or...
#
@generated function BlockIndex(::Type{RecursiveMatrix{T,M,N,L}}, ::Block{m,n}) where {T,M,N,L,m,n}
    full_n_blocks, n_r = divrem(N, cutoff)
    full_m_blocks, m_r = divrem(M, cutoff)

    offset = sizeof(T) * ( dense_sub2ind((M,N), (m-1)*cutoff + 1, (n-1)*cutoff + 1 ) - 1 )

    if m <= full_m_blocks && n <= full_n_blocks
        return BlockIndex{T,cutoff,cutoff,cutoff2}( offset )
    elseif m <= full_m_blocks # we're on the far right
        return BlockIndex{T,cutoff,n_r,cutoff*n_r}( offset )
    elseif n <= full_n_blocks
        return BlockIndex{T,m_r,cutoff,m_r*cutoff}( offset )
    else
        return BlockIndex{T,m_r,n_r,m_r*n_r}( offset )
    end
end


@generated function BlockIndex(::Type{RecursiveMatrix{T,M,N,L}}, ::HalfBlock{m,n,h}) where {T,M,N,L,m,n,h}
    full_n_blocks, n_r = divrem(N, cutoff)
    full_m_blocks, m_r = divrem(M, cutoff)
    
    offset = sizeof(T) * ( dense_sub2ind((M,N), (m-1)*cutoff + 1, (n-1)*cutoff + (h-1)*halfcutoff + 1 ) - 1 )

    if m <= full_m_blocks && n <= full_n_blocks
        return BlockIndex{T,cutoff,halfcutoff,cutoff*halfcutoff}( offset )
    elseif m <= full_m_blocks # we're on the far right
        ndim = h == 1 ? min(n_r, halfcutoff) : n_r - halfcutoff
        return BlockIndex{T,cutoff,ndim,cutoff*ndim}( offset )
    elseif n <= full_n_blocks # we're on the bottom
        return BlockIndex{T,m_r,halfcutoff,m_r*halfcutoff}( offset )
    else
        ndim = h == 1 ? min(n_r, halfcutoff) : n_r - halfcutoff
        return BlockIndex{T,m_r,ndim,m_r*ndim}( offset )
    end
end



@generated function Base.getindex(A::MutableRecursiveMatrix{T}, i::BlockIndex{T,m,n,mn}) where {T,m,n,mn}
    quote 
        Base.@_inline_meta
        Base.unsafe_load(convert_point(A, StaticRecursiveMatrix{$T,$m,$n,$mn} ) + i.i)
    end
end
function Base.getindex(A::MutableRecursiveMatrix{T,M,N,L}, i::BlockIndex{T,cutoff,cutoff,cutoff2}) where {T,M,N,L}
    Base.@_inline_meta
    Base.unsafe_load(convert_point(A) + i.i)
end

"""
The setindex methods don't work too well...
"""
@generated function Base.setindex!(A::MutableRecursiveMatrix{T},
                                    val::StaticRecursiveMatrix{T,m,n,mn},
                                    i::BlockIndex{T,m,n,mn}) where {T,m,n,mn}
    quote
        Base.@_inline_meta
        Base.unsafe_store!(convert_point(A, StaticRecursiveMatrix{$T,$m,$n,$mn} ) + i.i, val); nothing
    end
end

function Base.setindex!(A::MutableRecursiveMatrix{T,M,N,L},
                                    val::StaticRecursiveMatrix{T,cutoff,cutoff,cutoff2},
                                    i::BlockIndex{T,cutoff,cutoff,cutoff2}) where {T,M,N,L}
    Base.@_inline_meta
    Base.unsafe_store!(convert_point(A) + i.i, val); nothing
end

@generated function point(A::MutableRecursiveMatrix{T}, i::BlockIndex{T,m,n,mn}) where {T,m,n,mn}
    :(convert_point(A, StaticRecursiveMatrix{$T,$m,$n,$mn} ) + i.i)
end
function point(A::MutableRecursiveMatrix{T}, i::BlockIndex{T,cutoff,cutoff,cutoff2}) where T
    convert_point(A) + i.i
end

# function setindex_block_quote(::Type{T}, ::Type{V}, full_m_blocks, mdim, m_r, full_n_blocks, ndim, n_r, m, n) where {T,V}
#     if m <= full_m_blocks && n <= full_n_blocks
#         Lout = mdim*ndim
#         if mdim == ndim == cutoff
#             q = quote
#                     Base.@_inline_meta
#                     Base.unsafe_store!(convert_point(A) + $(sizeof(T)*( (m-1)*Lout + (n-1)*M*ndim)), val )
#                 end
#         else
#             q = quote
#                     Base.@_inline_meta
#                     Base.unsafe_store!(convert_point(A, Ptr{V}) + $(sizeof(T)*( (m-1)*Lout + (n-1)*M*ndim) ), val )
#             end
#         end

#     elseif m <= full_m_blocks # we're on the far right
#         Lout = mdim * n_r
#         q = quote
#                 Base.@_inline_meta
#                 Base.unsafe_store!(convert_point( A, Ptr{V}) + $(sizeof(T)*( (m-1)*Lout + (n-1)*M*ndim) ). val )
#             end

#     elseif n <= full_n_blocks
#         Lout = m_r * ndim
#         q = quote
#                 Base.@_inline_meta
#                 Base.unsafe_store!(convert_point( A, Ptr{V}) + $(sizeof(T)*( (m-1)*ndim*mdim + (n-1)*M*ndim) ), val )
#             end

#     else
#         Lout = m_r*n_r
#         q = quote
#                 Base.@_inline_meta
#                 Base.unsafe_store!(convert_point( A, Ptr{V}) + $(sizeof(T)*( L - Lout )), val )
#             end

#     end    
#     q
# end

# @generated function Base.setindex!(A::MutableRecursiveMatrix{T,M,N,L}, val::V, ::HalfBlock{m,n}) where {T,M,N,L,V,m,n}
#     full_m_blocks, incomplete_m_block = divrem(M, cutoff)
#     full_n_blocks, incomplete_n_block = divrem(N, halfcutoff)
#     setindex_block_quote(T, V, full_m_blocks, cutoff, incomplete_m_block, full_n_blocks, halfcutoff, incomplete_n_block, m, n)
# end

# @generated function Base.setindex!(A::MutableRecursiveMatrix{T,M,N,L}, val::V, ::Block{m,n}) where {T,M,N,L,V,m,n}
#     full_m_blocks, incomplete_m_block = divrem(M, cutoff)
#     full_n_blocks, incomplete_n_block = divrem(N, cutoff)
#     setindex_block_quote(T, V, full_m_blocks, cutoff, incomplete_m_block, full_n_blocks, cutoff, incomplete_n_block, m, n)
# end

# rmatA8 = RecursiveMatrix{Float64,8,8,abs2(8)}(Ref(ntuple(i -> randn(), Val(abs2(8)))));
# rmatB8 = RecursiveMatrix{Float64,8,8,abs2(8)}(Ref(ntuple(i -> randn(), Val(abs2(8)))));
# rmatC8 = RecursiveMatrix{Float64,8,8,abs2(8)}(Ref{NTuple{64,Float64}}());
# rmatA12 = RecursiveMatrix{Float64,12,12,abs2(12)}(Ref(ntuple(i -> randn(), Val(abs2(12)))));
# rmatB12 = RecursiveMatrix{Float64,12,12,abs2(12)}(Ref(ntuple(i -> randn(), Val(abs2(12)))));
# rmatC12 = RecursiveMatrix{Float64,12,12,abs2(12)}(Ref{NTuple{144,Float64}}());
# rmatA16 = RecursiveMatrix{Float64,16,16,abs2(16)}(Ref(ntuple(i -> randn(), Val(abs2(16)))));
# rmatB16 = RecursiveMatrix{Float64,16,16,abs2(16)}(Ref(ntuple(i -> randn(), Val(abs2(16)))));
# rmatC16 = RecursiveMatrix{Float64,16,16,abs2(16)}(Ref{NTuple{256,Float64}}());