

function gen_chol_quote(T,N,N2)
    q, qa = initial_quote(N2, :U)

    for i ∈ 1:N
        lti = ltriangle(i)
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
            ltj = ltriangle(j)
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

@generated function LinearAlgebra.chol(U::SymmetricMatrix{T,N,N2}) where {T,N,N2}
    gen_chol_quote(UpperTriangularMatrix{T,N,N2},N,N2)
end

function gen_revchol_quote(T,N,N2)
    q, qa = initial_quote(N2, :U)

    for i ∈ N:-1:1
        lti = ltriangle(i)
        bti = lti+i
        Ui_i = Symbol(:Ui_, bti)
        U_i = Symbol(:U_, bti)
        push!(qa, :($Ui_i = $U_i))
        for j ∈ i+1:N
            Uj_i = Symbol(:Ui_, ltriangle(j) + i)
            push!(qa, :($Ui_i -= $Uj_i*$Uj_i))
        end
        push!(qa, :($Ui_i = sqrt($Ui_i)))


        for j ∈ 1:i-1
            ltj = ltriangle(j)
            Uj_i = Symbol(:Ui_, lti+j)
            U_i  = Symbol(:U_,  lti+j)
            push!(qa, :($Uj_i = $U_i))
            for k ∈ i+1:N
                ltk = ltriangle(k)
                Ujk, Uik = Symbol(:Ui_, ltk+j), Symbol(:Ui_, ltk+i)
                push!(qa, :($Uj_i -= $Ujk * $Uik))
            end
            push!(qa, :($Uj_i /= $Ui_i) )
        end
    end

    push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    q

end

@generated function revchol(U::SymmetricMatrix{T,N,N2}) where {T,N,N2}
    # q = gen_revchol_quote(UpperTriangularMatrix{T,N,N2},N,N2)
    # q
    gen_revchol_quote(UpperTriangularMatrix{T,N,N2},N,N2)
end
