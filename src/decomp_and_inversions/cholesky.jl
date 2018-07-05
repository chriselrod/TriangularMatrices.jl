

# function gen_ip_chol_quote!(qa,N,L,U = :U, S = :S, extract = extract_symbol)

#     for i ∈ 1:N
#         lti = small_triangle(i)
#         bti = lti+i
#         Ui_i = extract(U, bti)
#         U_i = extract(S, bti)
#         if i > 1
#             Uj_i = extract(U, lti+1)
#             Ui_i_expr = :($Ui_i = sqrt($U_i - $Uj_i*$Uj_i))
#             for j ∈ 2:i - 1
#                 Uj_i = extract(U, lti+j)
#                 push!(Ui_i_expr.args[2].args[2].args, :($Uj_i*$Uj_i))
#             end
#             push!(qa, Ui_i_expr)
#         else
#             push!(qa, :($Ui_i = sqrt($U_i)))
#         end


#         for j ∈ i+1:N
#             ltj = small_triangle(j)
#             Uj_i = extract(U, ltj+i)
#             U_i = extract(S, ltj+i)
#             if i > 1
#                 Ujk, Uik = extract(U, ltj+k), extract(U, lti+k)
#                 Uj_i_expr = :($Uj_i = ($U_i - $Ujk * $Uik) / $Ui_i) 
#                 for k ∈ 2:i - 1
#                     Ujk, Uik = extract(U, ltj+k), extract(U, lti+k)
#                     push!(Uj_i_expr.args[2].args[2].args, :($Ujk * $Uik))
#                 end
#                 push!(qa, Uj_i_expr )
#             else
#                 push!(qa, :($U_i /= $Ui_i) )
#             end
#         end
#     end

#     # push!(q.args, U)
#     qa
# end

# function gen_ip_chol_quote!(qa,N,L,U = :U, S = :S, extract = extract_symbol)

#     for i ∈ 1:N
#         lti = small_triangle(i)
#         bti = lti+i
#         Ui_i = extract(U, bti)
#         U_i = extract(S, bti)
#         if i == 1
#             push!(qa, :($Ui_i = sqrt($U_i)))
#         elseif i == 2
#             Uj_i = extract(U, lti+1)
#             push!(qa, :($Ui_i = sqrt($U_i - $Uj_i*$Uj_i)))
#         else
#             Uj_i = extract(U, lti+1)
#             push!(qa, :($Ui_i = $U_i - $Uj_i*$Uj_i))
#             for j ∈ 2:i - 1
#                 Uj_i = extract(U, lti+j)
#                 push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
#             end
#             push!(qa, :($Ui_i = sqrt($Ui_i)))
#         end


#         for j ∈ i+1:N
#             ltj = small_triangle(j)
#             Uj_i = extract(U, ltj+i)
#             U_i = extract(S, ltj+i)
#             push!(qa, :($Uj_i = $U_i))
#             for k ∈ 1:i - 1
#                 Ujk, Uik = extract(U, ltj+k), extract(U, lti+k)
#                 push!(qa, :($Uj_i -= $Ujk * $Uik))
#             end
#             push!(qa, :($Uj_i /= $Ui_i) )
#         end
#     end

#     # push!(q.args, U)
#     qa
# end
function gen_ip_chol_quote!(qa,::Type{T},N::Int,L,U = :U, S = :S, extract = extract_symbol) where T

    chunksize = 32 ÷ sizeof(T)
    # @show qa

    for i ∈ 1:N
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = extract(U, bti)
        U_i = extract(S, bti)
        if i == 1
            push!(qa, :($Ui_i = sqrt($U_i)))
        else
            nchunks, r = divrem(i-1, chunksize)
            if r == 0
                push!(qa, :($Ui_i = $U_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1:chunksize ]...)  ))
                for chunk ∈ 2:nchunks
                    push!(qa, :($Ui_i = $Ui_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1+(chunk-1)*chunksize:chunk*chunksize ]...)  ))
                end
                
            else
                push!(qa, :($Ui_i = $U_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1:r ]...)  ))
                for chunk ∈ 1:nchunks
                    push!(qa, :($Ui_i = $Ui_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1+r+(chunk-1)*chunksize:chunk*chunksize+r ]...)  ))
                end
            end
            # for j ∈ 1:i - 1
            #     Uj_i = extract(U, lti+j)
            #     push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
            # end
            push!(qa, :($Ui_i = sqrt($Ui_i)))
        end


        for j ∈ i+1:N
            ltj = small_triangle(j)
            Uj_i = extract(U, ltj+i)
            U_i = extract(S, ltj+i)
            if i == 1
                push!(qa, :($Uj_i = $U_i / $Ui_i) )
            else


                nchunks, r = divrem(i-1, chunksize)
                if r == 0
                    push!(qa, :($Uj_i = $U_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1:chunksize ]...)  ))
                    for chunk ∈ 2:nchunks
                        push!(qa, :($Uj_i= $Uj_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1+(chunk-1)*chunksize:chunk*chunksize ]...)  ))
                    end
                    
                else
                    push!(qa, :($Uj_i = $U_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1:r ]...)  ))
                    for chunk ∈ 1:nchunks
                        push!(qa, :($Uj_i = $Uj_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1+r+(chunk-1)*chunksize:chunk*chunksize+r ]...)  ))
                    end
                end

                # push!(qa, :($Uj_i = $U_i))
                # for k ∈ 1:i - 1
                #     Ujk, Uik = extract(U, ltj+k), extract(U, lti+k)
                #     push!(qa, :($Uj_i -= $Ujk * $Uik))
                # end
                push!(qa, :($Uj_i /= $Ui_i) )
            end
        end
    end

    # push!(q.args, U)
    qa
end
function gen_ip_choldet_quote!(qa,::Type{T},N::Int,L,U = :U, S = :S, extract = extract_symbol) where T

    chunksize = 32 ÷ sizeof(T)
    # @show qa
    push!(qa, :(root_det = $(one(T))))
    for i ∈ 1:N
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = extract(U, bti)
        U_i = extract(S, bti)
        if i == 1
            push!(qa, :($Ui_i = sqrt($U_i)))
        else
            nchunks, r = divrem(i-1, chunksize)
            if r == 0
                push!(qa, :($Ui_i = $U_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1:chunksize ]...)  ))
                for chunk ∈ 2:nchunks
                    push!(qa, :($Ui_i = $Ui_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1+(chunk-1)*chunksize:chunk*chunksize ]...)  ))
                end
                
            else
                push!(qa, :($Ui_i = $U_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1:r ]...)  ))
                for chunk ∈ 1:nchunks
                    push!(qa, :($Ui_i = $Ui_i - +$([:( $(extract(U, lti+j)) * $(extract(U, lti+j))) for j ∈ 1+r+(chunk-1)*chunksize:chunk*chunksize+r ]...)  ))
                end
            end
            # for j ∈ 1:i - 1
            #     Uj_i = extract(U, lti+j)
            #     push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
            # end
            push!(qa, :($Ui_i = sqrt($Ui_i)))
            push!(qa, :(root_det *= $Ui_i))
        end


        for j ∈ i+1:N
            ltj = small_triangle(j)
            Uj_i = extract(U, ltj+i)
            U_i = extract(S, ltj+i)
            if i == 1
                push!(qa, :($Uj_i = $U_i / $Ui_i) )
            else


                nchunks, r = divrem(i-1, chunksize)
                if r == 0
                    push!(qa, :($Uj_i = $U_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1:chunksize ]...)  ))
                    for chunk ∈ 2:nchunks
                        push!(qa, :($Uj_i= $Uj_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1+(chunk-1)*chunksize:chunk*chunksize ]...)  ))
                    end
                    
                else
                    push!(qa, :($Uj_i = $U_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1:r ]...)  ))
                    for chunk ∈ 1:nchunks
                        push!(qa, :($Uj_i = $Uj_i - +$([:( $(extract(U, ltj+k)) * $(extract(U, lti+k))) for k ∈ 1+r+(chunk-1)*chunksize:chunk*chunksize+r ]...)  ))
                    end
                end

                # push!(qa, :($Uj_i = $U_i))
                # for k ∈ 1:i - 1
                #     Ujk, Uik = extract(U, ltj+k), extract(U, lti+k)
                #     push!(qa, :($Uj_i -= $Ujk * $Uik))
                # end
                push!(qa, :($Uj_i /= $Ui_i) )
            end
        end
    end

    # push!(q.args, U)
    qa
end
# function gen_ip_chol_quote!(qa,N,L,U = :U, S = :S, extract = extract_symbol)

#     for i ∈ 1:N
#         lti = small_triangle(i)
#         bti = lti+i
#         Ui_i = extract(U, bti)
#         U_i = extract(S, bti)
#         push!(qa, :($Ui_i = $U_i))
#         for j ∈ 1:i - 1
#             Uj_i = extract(U, lti+j)
#             push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
#         end
#         push!(qa, :($Ui_i = sqrt($Ui_i)))


#         for j ∈ i+1:N
#             ltj = small_triangle(j)
#             Uj_i = extract(U, ltj+i)
#             U_i = extract(S, ltj+i)
#             push!(qa, :($Uj_i = $U_i))
#             for k ∈ 1:i - 1
#                 Ujk, Uik = extract(U, ltj+k), extract(U, lti+k)
#                 push!(qa, :($Uj_i -= $Ujk * $Uik))
#             end
#             push!(qa, :($Uj_i /= $Ui_i) )
#         end
#     end

#     # push!(q.args, U)
#     qa
# end
function gen_extracted_ip_chol_quote!(qa, ::Type{T}, N, L, U = :U, S = :S) where T
    for i ∈ 1:L
        push!(qa, :($(id_symbol(S, i)) = $(S)[$i]))
    end
    gen_ip_chol_quote!(qa, T, N, L, U, S, id_symbol)
    for i ∈ 1:L
        push!(qa, :($(U)[$i] = $(id_symbol(U, i))))
    end
end
# function chol_quote(S = :S, ::Type{T} = Float64) where T
#     q, qa = create_quote()
#     for i ∈ 1:L
#         push!(qa, :($(id_symbol(S, i)) = $(S)[$i]))
#     end
#     gen_ip_chol_quote!(qa, T, N, L, :U, :S, id_symbol)
#     q
# end
@generated function LinearAlgebra.chol(S::AbstractSymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    gen_extracted_ip_chol_quote!(qa, T, N, L, :S, :S)
    push!(q.args, :(StaticUpperTriangularMatrix{$T,$N,$L}( ( @ntuple $L S )  )))
    q
end
@generated function cholesky!(S::SymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    gen_extracted_ip_chol_quote!(qa, T, N, L, :S, :S)
    push!(q.args, :(UpperTriangularMatrix{$T,$N,$L}(S.data)))
    q
end
@generated function cholesky!(U::UpperTriangularMatrix{T,N,L}, S::AbstractSymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    gen_extracted_ip_chol_quote!(qa, T, N, L, :U, :S)
    push!(q.args, :U)
    q
end

@generated function choldet!(S::SymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    extract_linear!(qa, L, :S)
    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl
    insert_linear!(qa, L, :S, :Si)
    push!(q.args, :(root_det))
    q
end
@generated function choldet!(S::LinearAlgebra.Symmetric{RecursiveMatrix{T,N,N,N2}}) where {T,N,N2}
    q, qa = create_quote()
    L = big_triangle(N)
    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( $(Symbol(:S_, ind) = S.data[$nr, $nc]) ))
    end

    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl

    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( S.data[$nr, $nc] = $(Symbol(:Si_, ind)) ))
    end
    push!(q.args, :(root_det))
    q
end
@generated function choldet!(S::SymmetricMatrix{T,M,L2}, ::Val{N}) where {T,M,N,L2}
    @assert N <= M
    L = big_triangle(N)
    q, qa = create_quote()
    extract_linear!(qa, L, :S)
    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl
    insert_linear!(qa, L, :S, :Si)
    push!(q.args, :(root_det))
    q
end
@generated function choldet!(S::LinearAlgebra.Symmetric{RecursiveMatrix{T,M,M,N2}}, ::Val{N}) where {T,M,N,N2}
    @assert N <= M
    q, qa = create_quote()
    L = big_triangle(N)
    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( $(Symbol(:S_, ind) = S.data[$nr, $nc]) ))
    end

    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl

    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( S.data[$nr, $nc] = $(Symbol(:Si_, ind)) ))
    end
    push!(q.args, :(root_det))
    q
end
@generated function choldet!(Si::Union{SymmetricMatrix{T,N,L1},Symmetric{RecursiveMatrix{T,N,N,L1}}},
                            S::Union{SymmetricMatrix{T,N,L2},Symmetric{RecursiveMatrix{T,N,N,L2}}}) where {T,N,L1,L2}
    q, qa = create_quote()
    if S <: Symmetric
        L_triangle = big_triangle(N)
        ind = 0
        for nc ∈ 1:N, nr ∈ 1:nc
            ind += 1
            push!(qa, :( $(Symbol(:S_, ind) = S.data[$nr, $nc]) ))
        end
    else
        L_triangle = L1
        extract_linear!(qa, L1, :S)
    end

    gen_ip_choldet_quote!(qa, T, N, L_triangle, :S, :S, id_symbol) #in this file, cholesky.jl #creates root_det variable, root of the determinant of S.

    if Si <: Symmetric
        ind = 0
        for nc ∈ 1:N, nr ∈ 1:nc
            ind += 1
            push!(qa, :( Si.data[$nr, $nc] = $(Symbol(:Si_, ind)) ))
        end
    else
        insert_linear!(qa, L_triangle, :Si)
    end
    push!(q.args, :(root_det))
    q
end
@generated function inv!(S::SymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    extract_linear!(qa, L, :S)

    gen_ip_chol_quote!(qa, T, N, L, :S, :S, id_symbol)

    insert_linear!(qa, L, :S)
    push!(q.args, :(S))
    q
end
@generated function invdet!(S::SymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    extract_linear!(qa, L, :S)

    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl
    gen_ip_inv_quote!(qa, N, L, :S, :S, id_symbol) #in inverse.jl # triangular matrix inversion
    gen_upper_xxt_quote!(qa,T,N,L, :S, :Si, id_symbol) #in arithmetic.jl

    insert_linear!(qa, L, :S, :Si)
    push!(q.args, :(root_det))
    q
end
@generated function invdet!(S::SymmetricMatrix{T,M,L2}, ::Val{N}) where {T,N,L2,M}
    @assert N <= M
    L = big_triangle(N)
    q, qa = create_quote()
    extract_linear!(qa, L, :S)

    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl
    gen_ip_inv_quote!(qa, N, L, :S, :S, id_symbol) #in inverse.jl # triangular matrix inversion
    gen_upper_xxt_quote!(qa,T,N,L, :S, :Si, id_symbol) #in arithmetic.jl

    insert_linear!(qa, L, :S, :Si)
    push!(q.args, :(root_det))
    q
end
@generated function invdet!(S::LinearAlgebra.Symmetric{RecursiveMatrix{T,N,N,N2}}) where {T,N,N2}
    q, qa = create_quote()
    L = big_triangle(N)
    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( $(Symbol(:S_, ind) = S.data[$nr, $nc]) ))
    end
    # extract_linear!(qa, L, :S)

    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl
    gen_ip_inv_quote!(qa, N, L, :S, :S, id_symbol) #in inverse.jl # triangular matrix inversion
    gen_upper_xxt_quote!(qa,T,N,L, :S, :Si, id_symbol) #in arithmetic.jl

    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( S.data[$nr, $nc] = $(Symbol(:Si_, ind)) ))
    end
    push!(q.args, :(root_det))
    q
end
@generated function invdet!(S::LinearAlgebra.Symmetric{RecursiveMatrix{T,M,M,N2}}, ::Val{N}) where {T,N,N2,M}
    @assert N <= M
    q, qa = create_quote()
    L = big_triangle(N)
    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( $(Symbol(:S_, ind) = S.data[$nr, $nc]) ))
    end
    # extract_linear!(qa, L, :S)

    gen_ip_choldet_quote!(qa, T, N, L, :S, :S, id_symbol) #in this file, cholesky.jl
    gen_ip_inv_quote!(qa, N, L, :S, :S, id_symbol) #in inverse.jl # triangular matrix inversion
    gen_upper_xxt_quote!(qa,T,N,L, :S, :Si, id_symbol) #in arithmetic.jl

    ind = 0
    for nc ∈ 1:N, nr ∈ 1:nc
        ind += 1
        push!(qa, :( S.data[$nr, $nc] = $(Symbol(:Si_, ind)) ))
    end
    push!(q.args, :(root_det))
    q
end
@generated function invdet!(Si::Union{SymmetricMatrix{T,N,L1},Symmetric{RecursiveMatrix{T,N,N,L1}}},
                            S::Union{SymmetricMatrix{T,N,L2},Symmetric{RecursiveMatrix{T,N,N,L2}}}) where {T,N,L1,L2}
    q, qa = create_quote()
    if S <: Symmetric
        L_triangle = big_triangle(N)
        ind = 0
        for nc ∈ 1:N, nr ∈ 1:nc
            ind += 1
            push!(qa, :( $(Symbol(:S_, ind) = S.data[$nr, $nc]) ))
        end
    else
        L_triangle = L1
        extract_linear!(qa, L1, :S)
    end

    gen_ip_choldet_quote!(qa, T, N, L_triangle, :S, :S, id_symbol) #in this file, cholesky.jl #creates root_det variable, root of the determinant of S.
    gen_ip_inv_quote!(qa, N, L_triangle, :S, :S, id_symbol) #in inverse.jl # triangular matrix inversion
    gen_upper_xxt_quote!(qa,T,N, L_triangle, :S, :Si, id_symbol) #in arithmetic.jl

    if Si <: Symmetric
        ind = 0
        for nc ∈ 1:N, nr ∈ 1:nc
            ind += 1
            push!(qa, :( Si.data[$nr, $nc] = $(Symbol(:Si_, ind)) ))
        end
    else
        insert_linear!(qa, L_triangle, :Si)
    end
    push!(q.args, :(root_det))
    q
end
# @generated function cholesky!(U::UpperTriangularMatrix{T,N,L}, S::AbstractSymmetricMatrix{T,N,L}) where {T,N,L}
#     q, qa = create_quote()
#     gen_extracted_ip_chol_quote!(qa, T, N, L, :U, :S)
#     push!(q.args, :U)
#     q
# end

function gen_ip_revchol_quote!(qa, ::Type{T}, N,L, U = :U, S = :S, extract = extract_symbol) where T

    for i ∈ N:-1:1
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = extract(U, bti)
        U_i = extract(S, bti)
        if i == N
            push!(qa, :($Ui_i = sqrt($U_i)))
        else
            push!(qa, :($Ui_i = $U_i))
            for j ∈ i+1:N
                Uj_i = extract(U, small_triangle(j) + i)
                push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
            end
            push!(qa, :($Ui_i = sqrt($Ui_i)))
        end

        for j ∈ 1:i-1
            ltj = small_triangle(j)
            Uj_i = extract(U, lti+j)
            U_i  = extract(S, lti+j)
            if i == N
                push!(qa, :($Uj_i = $U_i / $Ui_i) )
            else
                push!(qa, :($Uj_i = $U_i))
                for k ∈ i+1:N
                    ltk = small_triangle(k)
                    Ujk, Uik = extract(U, ltk+j), extract(U, ltk+i)
                    push!(qa, :($Uj_i -= $Ujk * $Uik))
                end
                push!(qa, :($Uj_i /= $Ui_i) )
            end
        end
    end
    # push!(q.args, :($T( ( @ntuple $L Ui ) )))
    qa

end

"""
This directly calculates the inverse of the Cholesky decomposition of the inverse.
Or, for positive definite input matrix S, it calculates the upper triangular matrix U such that
S = U * U'
"""
@generated function reverse_cholesky(S::SymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    for i ∈ 1:L
        push!(qa, :($(id_symbol(:S, i)) = S[$i]))
    end
    gen_ip_revchol_quote!(qa, T, N, L, :U, :S, id_symbol)
    push!(q.args, :(UpperTriangularMatrix{$T,$N,$L}( ( @ntuple $L U ) )))
    q
end

@generated function reverse_cholesky!(U::UpperTriangularMatrix{T,N,L}, S::AbstractSymmetricMatrix{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    for i ∈ 1:L
        push!(qa, :($(id_symbol(:S, i)) = S[$i]))
    end
    gen_ip_revchol_quote!(qa, T, N, L, :U, :S, id_symbol)
    for i ∈ 1:L
        push!(qa, :(U[$i] = $(id_symbol(:U, i))))
    end
    push!(q.args, :U)
    q
end
