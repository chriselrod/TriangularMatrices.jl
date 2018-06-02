

function gen_chol_quote(T,N,N2)
    q, qa = initial_quote(N2, :U)

    for i ∈ 1:N
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = Symbol(:Ui_, bti)
        U_i = Symbol(:U_, bti)
        push!(qa, :($Ui_i = $U_i))
        for j ∈ 1:i - 1
            Uj_i = Symbol(:Ui_, lti + j)
            push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
        end
        push!(qa, :($Ui_i = sqrt($Ui_i)))


        for j ∈ i+1:N
            ltj = small_triangle(j)
            Uj_i = Symbol(:Ui_, ltj+i)
            U_i = Symbol(:U_, ltj+i)
            push!(qa, :($Uj_i = $U_i))
            for k ∈ 1:i - 1
                Ujk, Uik = Symbol(:Ui_, ltj+k), Symbol(:Ui_, lti+k)
                push!(qa, :($Uj_i -= $Ujk * $Uik))
            end
            push!(qa, :($Uj_i /= $Ui_i) )
        end
    end

    push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    q
end


function gen_ip_chol_quote!(qa,N,N2,U = :U, S = :S, extract = extract_symbol)

    for i ∈ 1:N
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = extract(U, bti)
        U_i = extract(S, bti)
        push!(qa, :($Ui_i = $U_i))
        for j ∈ 1:i - 1
            Uj_i = extract(U, lti+j)
            push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
        end
        push!(qa, :($Ui_i = sqrt($Ui_i)))


        for j ∈ i+1:N
            ltj = small_triangle(j)
            Uj_i = extract(U, ltj+i)
            U_i = extract(S, ltj+i)
            push!(qa, :($Uj_i = $U_i))
            for k ∈ 1:i - 1
                Ujk, Uik = extract(U, ltj+k), extract(U, lti+k)
                push!(qa, :($Uj_i -= $Ujk * $Uik))
            end
            push!(qa, :($Uj_i /= $Ui_i) )
        end
    end

    # push!(q.args, U)
    qa
end
function gen_extracted_ip_chol_quote!(qa, N, N2, U = :U, S = :S)
    for i ∈ 1:N2
        push!(qa, :($(id_symbol(S, i)) = $(S)[$i]))
    end
    gen_ip_chol_quote!(qa, N, N2, U, S, id_symbol)
    for i ∈ 1:N2
        push!(qa, :($(U)[$i] = $(id_symbol(U, i))))
    end
end


@generated function LinearAlgebra.chol(U::AbstractSymmetricMatrix{T,N,L2}) where {T,N,L2}
    gen_chol_quote(UpperTriangularMatrix{T,N,L2},N,L2)
end
@generated function cholesky!(S::SymmetricMMatrix{T,N,N2}) where {T,N,N2}
    q, qa = create_quote()
    gen_extracted_ip_chol_quote!(qa, N, N2, :S, :S)
    push!(q.args, :(UpperTriangularMatrix{$T,$N,$N2}(S.data)))
    q
end
@generated function cholesky!(U::UpperTriangularMMatrix{T,N,N2}, S::AbstractSymmetricMatrix{T,N,N2}) where {T,N,N2}
    q, qa = create_quote()
    gen_extracted_ip_chol_quote!(qa, N, N2, :U, :S)
    push!(q.args, :U)
    q
end

function gen_revchol_quote(T,N,N2)
    q, qa = initial_quote(N2, :U)

    for i ∈ N:-1:1
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = Symbol(:Ui_, bti)
        U_i = Symbol(:U_, bti)
        push!(qa, :($Ui_i = $U_i))
        for j ∈ i+1:N
            Uj_i = Symbol(:Ui_, small_triangle(j) + i)
            push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
        end
        push!(qa, :($Ui_i = sqrt($Ui_i)))


        for j ∈ 1:i-1
            ltj = small_triangle(j)
            Uj_i = Symbol(:Ui_, lti+j)
            U_i  = Symbol(:U_,  lti+j)
            push!(qa, :($Uj_i = $U_i))
            for k ∈ i+1:N
                ltk = small_triangle(k)
                Ujk, Uik = Symbol(:Ui_, ltk+j), Symbol(:Ui_, ltk+i)
                push!(qa, :($Uj_i -= $Ujk * $Uik))
            end
            push!(qa, :($Uj_i /= $Ui_i) )
        end
    end

    push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    q

end
function gen_ip_revchol_quote!(qa,N,N2, U = :U, S = :S, extract = extract_symbol)


    for i ∈ N:-1:1
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = extract(U, bti)
        U_i = extract(S, bti)
        push!(qa, :($Ui_i = $U_i))
        for j ∈ i+1:N
            Uj_i = extract(U, small_triangle(j) + i)
            push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
        end
        push!(qa, :($Ui_i = sqrt($Ui_i)))


        for j ∈ 1:i-1
            ltj = small_triangle(j)
            Uj_i = extract(U, lti+j)
            U_i  = extract(S, lti+j)
            push!(qa, :($Uj_i = $U_i))
            for k ∈ i+1:N
                ltk = small_triangle(k)
                Ujk, Uik = extract(U, ltk+j), extract(U, ltk+i)
                push!(qa, :($Uj_i -= $Ujk * $Uik))
            end
            push!(qa, :($Uj_i /= $Ui_i) )
        end
    end

    # push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    qa

end

@generated function reverse_cholesky(U::SymmetricMatrix{T,N,N2}) where {T,N,N2}
    # q = gen_revchol_quote(UpperTriangularMatrix{T,N,N2},N,N2)
    # q
    gen_revchol_quote(UpperTriangularMatrix{T,N,N2},N,N2)
end

@generated function reverse_cholesky!(U::UpperTriangularMMatrix{T,N,N2}, S::AbstractSymmetricMatrix{T,N,N2}) where {T,N,N2}
    q, qa = create_quote()
    for i ∈ 1:N2
        push!(qa, :($(id_symbol(:S, i)) = S[$i]))
    end
    gen_ip_revchol_quote!(qa, N, N2, :U, :S, id_symbol)
    for i ∈ 1:N2
        push!(qa, :(U[$i] = $(id_symbol(:U, i))))
    end
    push!(q.args, :U)
    q
end
