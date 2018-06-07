
"""
A: M x N
B: N x P
C: M x P

C = A * B
"""

# function simd_mul_quote!(qa, ::Type{T}, M, N, P, ta=n(), A = :A, tb=n(), B = :B, tc=n(), C = :C, extract = extract_symbol, insert = extract_symbol) where T
#     eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     chunk = 64 ÷ sizeof(T)
#     N4, Nr = divrem(N, chunk)
#     if Nr == 0
#         Nt = N4
#     else
#         Nt = N4+1
#     end

#     # for m = 1:M
#     #     for s = 1:N4
#     #         push!(qa, :( $(Symbol(A, :_simd, s, :_, m))  = Vec{$chunk,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
#     #     end
#     #     Nr == 0 || push!(qa, :( $(Symbol(A, :_simd, N4+1, :_, m))  = Vec{$Nr,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*N4:N ]... )) ) ))
#     # end

#     j = 1
#     for s = 1:N4
#         push!(qa, :( $(Symbol(B, :_simd, s))  = Vec{$chunk,$T}( $(Expr(:tuple, [eB(i,j) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
#     end
#     Nr == 0 || push!(qa, :( $(Symbol(B, :_simd, N4+1))  = Vec{$Nr,$T}( $(Expr(:tuple, [eB(i,j) for i = 1+chunk*N4:N ]... )) ) ))
#     for m = 1:M
#         for s = 1:N4
#             push!(qa, :( $(Symbol(A, :_simd, s, :_, m))  = Vec{$chunk,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
#         end
#         Nr == 0 || push!(qa, :( $(Symbol(A, :_simd, N4+1, :_, m))  = Vec{$Nr,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*N4:N ]... )) ) ))
#         if Nt == 1
#             push!(qa, :( $(eC(m, j)) =  sum( $(Symbol(A, :_simd, 1, :_, m)) * $(Symbol(B, :_simd, 1))  ) )   )
#         else
#             push!(qa, :( $(eC(m, j)) =  +$( [ :(sum( $(Symbol(A, :_simd, s, :_, m)) * $(Symbol(B, :_simd, s))  ) ) for s = 1:Nt]... )  ) )
#         end
#     end

#     for j = 2:P
#         for s = 1:N4
#             # @show j, s
#             # expr = Expr(:tuple, [eB(i,j) for i = 1+chunk*(s-1):chunk*s ]... )
#             # @show expr
#             push!(qa, :( $(Symbol(B, :_simd, s))  = Vec{$chunk,$T}( $(Expr(:tuple, [eB(i,j) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
# #            @show qa[end]
#         end
#         Nr == 0 || push!(qa, :( $(Symbol(B, :_simd, N4+1))  = Vec{$Nr,$T}( $(Expr(:tuple, [eB(i,j) for i = 1+chunk*N4:N ]... )) ) ))
#         for m = 1:M
#             if Nt == 1
#                 push!(qa, :( $(eC(m, j)) =  sum( $(Symbol(A, :_simd, 1, :_, m)) * $(Symbol(B, :_simd, 1))  ) )   )
#             else
#                 push!(qa, :( $(eC(m, j)) =  +$( [ :(sum( $(Symbol(A, :_simd, s, :_, m)) * $(Symbol(B, :_simd, s))  ) ) for s = 1:Nt]... )  ) )
#             end
#         end
#     end
# end
# function simd_mul_quote!(qa, ::Type{T}, M, N, P, ta, A, tb::n, B = :B, tc=n(), C = :C) where T
#     eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     chunk = 32 ÷ sizeof(T)
#     N4, Nr = divrem(N, chunk)
#     if Nr == 0
#         Nt = N4
#     else
#         Nt = N4+1
#     end

#     # for m = 1:M
#     #     for s = 1:N4
#     #         push!(qa, :( $(Symbol(A, :_simd, s, :_, m))  = Vec{$chunk,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
#     #     end
#     #     Nr == 0 || push!(qa, :( $(Symbol(A, :_simd, N4+1, :_, m))  = Vec{$Nr,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*N4:N ]... )) ) ))
#     # end

#     j = 1
#     for s = 0:N4-1
#         push!(qa, :( $(Symbol(B, :_simd, s+1))  =  $B[$(1+chunk*s), Vec{$chunk,$T}]) )
#     end
#     Nr == 0 || push!(qa, :( $(Symbol(B, :_simd, N4+1))  =  $B[$(1+chunk*N4), Vec{$Nr,$T}] ))
#     for m = 1:M
#         for s = 1:N4
#             push!(qa, :( $(Symbol(A, :_simd, s, :_, m))  = Vec{$chunk,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
#         end
#         Nr == 0 || push!(qa, :( $(Symbol(A, :_simd, N4+1, :_, m))  = Vec{$Nr,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*N4:N ]... )) ) ))
#         if Nt == 1
#             push!(qa, :( $(eC(m, j)) =  sum( $(Symbol(A, :_simd, 1, :_, m)) * $(Symbol(B, :_simd, 1))  ) )   )
#         else
#             push!(qa, :( $(eC(m, j)) =  +$( [ :(sum( $(Symbol(A, :_simd, s, :_, m)) * $(Symbol(B, :_simd, s))  ) ) for s = 1:Nt]... )  ) )
#         end
#     end

#     for j = 2:P
#         for s = 0:N4-1
#             # @show j, s
#             # expr = Expr(:tuple, [eB(i,j) for i = 1+chunk*(s-1):chunk*s ]... )
#             # @show expr
#             push!(qa, :( $(Symbol(B, :_simd, s+1))  = Vec{$chunk,$T}( $B[$(1+chunk*s), Val{$chunk}()]) ))
# #            @show qa[end]
#         end
#         Nr == 0 || push!(qa, :( $(Symbol(B, :_simd, N4+1))  = Vec{$Nr,$T}( $B[$(1+chunk*N4), Val{$Nr}()])  ))
#         for m = 1:M
#             if Nt == 1
#                 push!(qa, :( $(eC(m, j)) =  sum( $(Symbol(A, :_simd, 1, :_, m)) * $(Symbol(B, :_simd, 1))  ) )   )
#             else
#                 push!(qa, :( $(eC(m, j)) =  +$( [ :(sum( $(Symbol(A, :_simd, s, :_, m)) * $(Symbol(B, :_simd, s))  ) ) for s = 1:Nt]... )  ) )
#             end
#         end
#     end
# end

# function svec_mul_quote!(qa, ::Type{T}, M, N, P, ta=n(), A = :A, tb=n(), B = :B, tc=n(), C = :C, extract = extract_symbol, insert = extract_symbol) where T
#     eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     chunk = 64 ÷ sizeof(T)
#     N4, Nr = divrem(N, chunk)
#     if Nr == 0
#         Nt = N4
#     else
#         Nt = N4+1
#     end

#     for m = 1:M
#         for s = 1:N4
#             push!(qa, :( $(Symbol(A, :_simd, s, :_, m))  = SVector{$chunk,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*(s-1):chunk*s ]... )) )' ))
#         end
#         Nr == 0 || push!(qa, :( $(Symbol(A, :_simd, N4+1, :_, m))  = SVector{$Nr,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*N4:N ]... )) )' ))
#     end

#     # j = 1
#     # for s = 1:N4
#     #     push!(qa, :( $(Symbol(B, :_simd, s))  = SVector{$chunk,$T}( $(Expr(:tuple, [eB(i,j) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
#     # end
#     # Nr == 0 || push!(qa, :( $(Symbol(B, :_simd, N4+1))  = SVector{$Nr,$T}( $(Expr(:tuple, [eB(i,j) for i = 1+chunk*N4:N ]... )) ) ))
#     # for m = 1:M
#     #     for s = 1:N4
#     #         push!(qa, :( $(Symbol(A, :_simd, s, :_, m))  = SVector{$chunk,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*(s-1):chunk*s ]... )) )' ))
#     #     end
#     #     Nr == 0 || push!(qa, :( $(Symbol(A, :_simd, N4+1, :_, m))  = SVector{$Nr,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*N4:N ]... )) )' ))
#     #     if Nt == 1
#     #         push!(qa, :( $(eC(m, j)) =  $(Symbol(A, :_simd, 1, :_, m)) * $(Symbol(B, :_simd, 1))  )   )
#     #     else
#     #         push!(qa, :( $(eC(m, j)) =  +$( [ :( $(Symbol(A, :_simd, s, :_, m)) * $(Symbol(B, :_simd, s))  ) for s = 1:Nt]... )  ) )
#     #     end
#     # end

#     for j = 1:P
#         for s = 1:N4
#             # @show j, s
#             # expr = Expr(:tuple, [eB(i,j) for i = 1+chunk*(s-1):chunk*s ]... )
#             # @show expr
#             push!(qa, :( $(Symbol(B, :_simd, s))  =  B[$(1+chunk*(s-1)), SVector{$chunk,$T}] ))
# #            @show qa[end]
#         end
#         Nr == 0 || push!(qa, :( $(Symbol(B, :_simd, N4+1))  = SVector{$Nr,$T}( $(Expr(:tuple, [eB(i,j) for i = 1+chunk*N4:N ]... )) ) ))
#         for m = 1:M
#             if Nt == 1
#                 push!(qa, :( $(eC(m, j)) =  $(Symbol(A, :_simd, 1, :_, m)) * $(Symbol(B, :_simd, 1))   )   )
#             else
#                 push!(qa, :( $(eC(m, j)) =  +$( [ :( $(Symbol(A, :_simd, s, :_, m)) * $(Symbol(B, :_simd, s))   ) for s = 1:Nt]... )  ) )
#             end
#         end
#     end
# end

# function tup_mul_quote!(qa, ::Type{T}, M, N, P, ta=n(), A = :A, tb=n(), B = :B, tc=n(), C = :C, extract = extract_symbol, insert = extract_symbol) where T
#     eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     chunk = 64 ÷ sizeof(T)
#     N4, Nr = divrem(N, chunk)
#     if Nr == 0
#         Nt = N4
#     else
#         Nt = N4+1
#     end

#     VT = Base.VecElement{T}
#     # for m = 1:M
#     #     for s = 1:N4
#     #         push!(qa, :( $(Symbol(A, :_simd, s, :_, m))  = Vec{$chunk,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*(s-1):chunk*s ]... )) ) ))
#     #     end
#     #     Nr == 0 || push!(qa, :( $(Symbol(A, :_simd, N4+1, :_, m))  = Vec{$Nr,$T}( $(Expr(:tuple, [eA(m,i) for i = 1+chunk*N4:N ]... )) ) ))
#     # end

#     j = 1
#     for s = 1:N4
#         push!(qa, :( ($(Symbol(B, :_simd, s)))::NTuple{$chunk,$VT}  = ($(Expr(:tuple, [ :($VT($(eB(i,j)))) for i = 1+chunk*(s-1):chunk*s ]... )))::NTuple{$chunk,$VT} ) )
#     end
#     Nr == 0 || push!(qa, :( ($(Symbol(B, :_simd, N4+1)))::NTuple{$Nr,$VT}  =  ($(Expr(:tuple, [:($VT($(eB(i,j)))) for i = 1+chunk*N4:N ]... )) )::NTuple{$Nr,$VT}) )
#     for m = 1:M
#         for s = 1:N4
#             push!(qa, :( $(Symbol(A, :_simd, s, :_, m))::NTuple{$chunk,$VT}  = ($(Expr(:tuple, [eA(m,i)))) for i = 1+chunk*(s-1):chunk*s ]... )) )::NTuple{$chunk,$VT} ) )
#         end
#         Nr == 0 || push!(qa, :( ($(Symbol(A, :_simd, N4+1, :_, m)))::NTuple{$Nr,$VT}  =  ($(Expr(:tuple, [eA(m,i)))) for i = 1+chunk*N4:N ]... )) )::NTuple{$Nr,$VT}) )
#         if Nt == 1
#             push!(qa, :( $(eC(m, j)) =  sum( $(Symbol(A, :_simd, 1, :_, m)) .* $(Symbol(B, :_simd, 1))  ) )   )
#         else
#             push!(qa, :( $(eC(m, j)) =  +$( [ :(sum( $(Symbol(A, :_simd, s, :_, m)) .* $(Symbol(B, :_simd, s))  ) ) for s = 1:Nt]... )  ) )
#         end
#     end

#     for j = 2:P
#         for s = 1:N4
#             # @show j, s
#             # expr = Expr(:tuple, [eB(i,j) for i = 1+chunk*(s-1):chunk*s ]... )
#             # @show expr
#             push!(qa, :( $(Symbol(B, :_simd, s))  = $(Expr(:tuple, [eB(i,j)))) for i = 1+chunk*(s-1):chunk*s ]... )) ) )
# #            @show qa[end]
#         end
#         Nr == 0 || push!(qa, :( $(Symbol(B, :_simd, N4+1))  =  $(Expr(:tuple, [eB(i,j)))) for i = 1+chunk*N4:N ]... )) ) )
#         for m = 1:M
#             if Nt == 1
#                 push!(qa, :( $(eC(m, j)) =  sum( $(Symbol(A, :_simd, 1, :_, m)) .* $(Symbol(B, :_simd, 1))  ) )   )
#             else
#                 push!(qa, :( $(eC(m, j)) =  +$( [ :(sum( $(Symbol(A, :_simd, s, :_, m)) .* $(Symbol(B, :_simd, s))  ) ) for s = 1:Nt]... )  ) )
#             end
#         end
#     end
# end
# function chunk_mul_quote!(qa, M, N, P, ta=n(), A = :A, tb=n(), B = :B, tc=n(), C = :C, extract = extract_symbol, insert = extract_symbol)
#     eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     chunk = 4
#     N4, Nr = divrem(N, chunk)
#     for j = 1:P, i = 1:M
#         C_ij = eC(i, j)
#         if Nr > 0
#             push!(qa, :($C_ij = $( reduce((ex1,ex2) -> :(+($ex1,$ex2)), vcat(
#                 [ :(+$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1:Nr]...)) ],
#                 [ :(+$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+Nr+chunk*(k-1):Nr + chunk*k ]...)) for k = 1:N4 ] )  ) )  )  )
#         else
#             push!(qa, :($C_ij = $( reduce((ex1,ex2) -> :(+($ex1,$ex2)), 
#                 [ :(+$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+chunk*(k-1):chunk*k ]...)) for k = 1:N4 ] )  ) )  )
#         end
#     end
# end
# function mul_quote!(qa, M, N, P, ta=n(), A = :A, tb=n(), B = :B, tc=n(), C = :C, extract = extract_symbol, insert = extract_symbol)
#     eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     for j = 1:P, i = 1:M
#         push!(qa, :( $(eC(i, j)) = $(reduce((ex1,ex2) -> :(+($ex1,$ex2)), [ :( $(eA(i,k))*$(eB(k,j)) ) for k = 1:N ] )) ) )
#     end
# end
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
# function mul_quotev!(qa, ::Type{T}, M, N, P, ta=n(), A = :A, tb=n(), B = :B, tc=n(), C = :C,
#             extractA = extract_symbol, extractB = extract_symbol, insert = extract_symbol) where T
#     eA = (i,j) -> extractA(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extractB(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     # chunk = 4
#     # N4, Nr = divrem(N, chunk)
#     for j = 1:P
#         push!(qa, :( bv = SVector{$N,$T}( $(Expr(:tuple, [eB(n,j) for n ∈ 1:N]...)) )' ))
#         for i = 1:M
#             C_ij = eC(i, j)
#             push!(qa, :( $(eC(i, j)) = bv * SVector{$N,$T}( $(Expr(:tuple, [eA(i,n) for n ∈ 1:N]...))  )) )
#         end
#     end
# # end
# function gemm_quote!(qa, M, N, P, ta=n(), A = :A, tb=n(), B = :B, tc=n(), C = :C, extract = extract_symbol)
#     # eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     # eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     # eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     # chunk = 4
#     # N4, Nr = divrem(N, chunk)
#     # for j = 1:P, i = 1:M
#     #     C_ij = eC(i, j)
#     #     if Nr > 0
#     #         push!(qa, :($C_ij += $( reduce((ex1,ex2) -> :(+($ex1,$ex2)), vcat(
#     #             [ :(+$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1:Nr]...)) ],
#     #             [ :(+$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+Nr+chunk*(k-1):Nr + chunk*k ]...)) for k = 1:N4 ] )  ) )  )  )
#     #     else
#     #         push!(qa, :($C_ij += $( reduce((ex1,ex2) -> :(+($ex1,$ex2)), 
#     #             [ :(+$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+chunk*(k-1):chunk*k ]...)) for k = 1:N4 ] )  ) )  )
#     #     end
#     # end
#     eA = (i,j) -> extract(A, sub2ind(ta, (M, N), i, j))
#     eB = (i,j) -> extract(B, sub2ind(tb, (N, P), i, j))
#     eC = (i,j) -> insert(C, sub2ind(tc, (M, P), i, j))
#     chunk = 4
#     N4, Nr = divrem(N, chunk)
#     for j = 1:P, i = 1:M
#         C_ij = eC(i, j)
#         if Nr > 0
#             push!(qa, :($C_ij += +$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1:Nr]...) ) )
#             for k = 1:N4
#                 push!(qa, :($C_ij += +$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+Nr+chunk*(k-1):Nr + chunk*k ]...) ) )
#             end
#         else
#             push!(qa, :($C_ij += +$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1:chunk ]...) ) )
#             for k = 2:N4
#                 push!(qa, :($C_ij += +$([:( $(eA(i,r)) * $(eB(r,j))) for r ∈ 1+chunk*(k-1):chunk*k ]...) ) )
#             end
#         end
#     end
# end

# @generated function mul_unrolled_chunks!(::Size{sc}, c::StaticMatrix, ::Size{sa}, ::Size{sb}, a::StaticMatrix, b::StaticMatrix) where {sa, sb, sc}
#     if sb[1] != sa[2] || sa[1] != sc[1] || sb[2] != sc[2]
#         throw(DimensionMismatch("Tried to multiply arrays of size $sa and $sb and assign to array of size $sc"))
#     end

#     #vect_exprs = [:($(Symbol("tmp_$k2")) = partly_unrolled_multiply(A, B[:, $k2])) for k2 = 1:sB[2]]

#     # Do a custom b[:, k2] to return a SVector (an isbits type) rather than a mutable type. Avoids allocation == faster
#     tmp_type = SVector{sb[1], eltype(c)}
#     vect_exprs = [:($(Symbol("tmp_$k2")) = partly_unrolled_multiply($(Size(sa)), $(Size(sb[1])), a, $(Expr(:call, tmp_type, [Expr(:ref, :b, LinearIndices(sb)[i, k2]) for i = 1:sb[1]]...)))) for k2 = 1:sb[2]]

#     exprs = [:(c[$(LinearIndices(sc)[k1, k2])] = $(Symbol("tmp_$k2"))[$k1]) for k1 = 1:sa[1], k2 = 1:sb[2]]

#     return quote
#         @_inline_meta
#         @inbounds $(Expr(:block, vect_exprs...))
#         @inbounds $(Expr(:block, exprs...))
#     end
# end
# @generated function partly_unrolled_multiply(::Size{sa}, ::Size{sb}, a::StaticMatrix{<:Any, <:Any, Ta}, b::StaticArray{<:Any, Tb}) where {sa, sb, Ta, Tb}
#     if sa[2] != sb[1]
#         throw(DimensionMismatch("Tried to multiply arrays of size $sa and $sb"))
#     end

#     if sa[2] != 0
#         exprs = [reduce((ex1,ex2) -> :(+($ex1,$ex2)), [:(a[$(LinearIndices(sa)[k, j])]*b[$j]) for j = 1:sa[2]]) for k = 1:sa[1]]
#     else
#         exprs = [:(zero(promote_op(matprod,Ta,Tb))) for k = 1:sa[1]]
#     end

#     return quote
#         $(Expr(:meta,:noinline))
#         @inbounds return SVector(tuple($(exprs...)))
#     end
# end

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

function block_kernel(M, N, P, tA=false, tB=false, tC=false, eq = :(=), LA = M*N, LB = N*P, LC = M*P)
    A_min_dim, A_max_dim = minmax(M,N)
    B_min_dim, B_max_dim = minmax(N,P)
    if (√2 * A_min_dim < A_max_dim) && (√2 * B_min_dim < B_max_dim) #Both matrices are divided into four blocks.

    elseif

    end

end
function blockmull4x4()


end
function blockmull2x4()


end
function blockmull4x2()


end
function blockmull2x2_out1()


end
function blockmull2x2_out4()


end


"""
The dummy argument always gets optimized out.
However, when using code_llvm on the function and filling it with a real type (eg, Bool), it causes the SSA names to get incremented by 1 (they start as %0 for the first argument, %1 for the second...). In llvmcall, these values get incremented by one between the first in the body and the argument list, while they do not in the function code returned by llvmcall.
Maybe there is a solution that is less of a hack, but just passing an extra dummy argument to the code_llvm call to get them to line up was a simple solution.
"""
@generated function mul!(C::RecursiveMatrixOrTranpose{T,M,P,LC},
                        A::RecursiveMatrixOrTranpose{T,M,N,LA},
                        B::RecursiveMatrixOrTranpose{T,N,P,LB},
                        dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    if max(M,N,P) < cutoff # No recursion; we multiply.
        return mul_kernel(M, N, P, istransposed(A), istransposed(B), istransposed(C))
    else # Recursion.
        return block_kernel(M,N,P, istransposed(A), istransposed(B), istransposed(C))
    end
end

"""
By default, Julia refuses to emit SIMD instructions where aliasing is possible.
However, the pointer matrices are only used internally where I can guarantee that they wont alias.

Unlike C/C++, which have `restrict` compiler hints (or Fortran which simply assumes you aren't aliasing),
there's no way for us to make that promise to the compiler. So, instead, the approach is to use `code_llvm`
for the corresponding types where aliasing isn't possible, and then use this code with `llvmcall`.
"""
@generated function mul!(C::RecursivePointerMatrixOrTranpose{T,M,P,LC},
                        A::RecursivePointerMatrixOrTranpose{T,M,N,LA},
                        B::RecursivePointerMatrixOrTranpose{T,N,P,LB},
                        dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    iobuffer = IOBuffer()
    code_llvm(iobuffer, mul!, (dereftype(C), dereftype(A), dereftype(B), Bool))
    kernel_code = String(iobuffer)
    codestart = search(mulstring,"{\ntop:\n ")[end]
    quote
        Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{$T}, Tuple{Ptr{$T},Ptr{$T},Ptr{$T}}, C.data, A.data, B.data)
        C
    end
end

@generated function gemm!(C::RecursiveMatrixOrTranpose{T,M,P,LC},
                        A::RecursiveMatrixOrTranpose{T,M,N,LA},
                        B::RecursiveMatrixOrTranpose{T,N,P,LB},
                        dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    if M*N < abs2(cutoff) && N*P < abs2(cutoff) # No recursion; we multiply.
        return mul_kernel(M, N, P, istransposed(A), istransposed(B), istransposed(C), :(+=) )
    else # Recursion.
        return block_kernel(M,N,P, istransposed(A), istransposed(B), istransposed(C), :(+=) )
    end
end
@generated function gemm!(C::RecursivePointerMatrixOrTranpose{T,M,P,LC},
                        A::RecursivePointerMatrixOrTranpose{T,M,N,LA},
                        B::RecursivePointerMatrixOrTranpose{T,N,P,LB},
                        dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    iobuffer = IOBuffer()
    code_llvm(iobuffer, gemm!, (dereftype(C), dereftype(A), dereftype(B), Bool))
    kernel_code = String(iobuffer)
    codestart = search(mulstring,"{\ntop:\n ")[end]
    quote
        Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{$T}, Tuple{Ptr{$T},Ptr{$T},Ptr{$T}}, C.data, A.data, B.data)
        C
    end
end


# const As = [:(n()), :(MutableRecursiveMatrix{T,M,P,LC}),
#             :(t()), :(Adjoint{T,MutableRecursiveMatrix{T,P,M,LA}})]

# const Bs = [:(n()), :(MutableRecursiveMatrix{T,P,N,LC}),
#             :(t()), :(Adjoint{T,MutableRecursiveMatrix{T,N,P,LB}})]

# const Cs = [:(n()), :(MutableRecursiveMatrix{T,M,N,LC}),
#             :(t()), :(Adjoint{T,MutableRecursiveMatrix{T,N,M,LC}})]


# for (atrans,atype) ∈ As, (btrans,btype) ∈ Bs, (ctrans,ctype) ∈ Cs
#     @eval @generated function mul!(C::$(ctype), A::$(atype), B::$(btype)) where {T,M,N,P,LA,LB,LC}
#             q, qa = create_quote()
#             extract_linear!(qa, LA, :A)
#             extract_linear!(qa, LB, :B)
#             mul_quote!(qa, M, N, P, $(atrans), :A, $(btrans), :B, $(ctrans), :C, id_symbol)
#             insert_linear!(qa, LC, :C)
#             push!(q, :C)
#             q
#         end


#     @eval @generated function gemm!(C::$(ctype), A::$(atype), B::$(btype)) where {T,M,N,P,LA,LB,LC}
#             q, qa = create_quote()
#             extract_linear!(qa, LA, :A)
#             extract_linear!(qa, LB, :B)
#             gemm_quote!(qa, M, N, P, $(atrans), :A, $(btrans), :B, $(ctrans), :C, id_symbol)
#             insert_linear!(qa, LC, :C)
#             push!(q, :C)
#             q
#         end

# end

# function gemm!(C, A, B)
#     q, qa = create_quote()

# end