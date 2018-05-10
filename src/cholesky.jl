

function gen_chol_quote(T,N,N2)
    q = quote
        @fastmath @inbounds begin
            U_1 = U.data[1]
        end
    end
    @static if VERSION > v"0.6.9"
        qa = q.args[2].args[3].args[3].args
    else
        qa = q.args[2].args[2].args[2].args
    end
    for i ∈ 2:N2
        U_i = Symbol(:U_, i)
        push!(qa, :($U_i = U.data[$i]))
    end

    for i ∈ 1:N
        lti = ltriangle(i)
        bti = lti+i
        Ui_i = Symbol(:Ui_, bti)
        U_i = Symbol(:U_, bti)
        push!(qa, :($Ui_i = $U_i))
        for j ∈ 1:i - 1
            Uj_i = Symbol(:Ui_, lti + j)
            push!(qa, :($Ui_i -= $Uj_i*Uj_i))
        end
        push!(qa, :($Ui_i = sqrt($Ui_i)))


        for j ∈ i+1:N
            ltj = ltriangle(j)
            Uj_i = Symbol(:Ui_, ltj+i)
            U_i = Symbol(:U_, ltj+i)
            push!(qa, :($Uj_i = $U_i))
            for k = 1:i - 1
                push!(qa, :($Uj_i -= $Symbol(:Ui_, ltj+k) * $Symbol(:Ui_, lti+k)))
            end
            push!(qa, :($Uj_i /= $Ui_i) )
        end
    end

    push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    q

end

@generated function chol(A::SymmetricMatrix{T,N,N2}) where {T,N,N2}
    gen_chol_quote(UpperTriangularMatrix{T,N,N2},N,N2)
end

function gen_revchol_quote(T,N,N2)
    q = quote
        @fastmath @inbounds begin
            U_1 = U.data[1]
        end
    end
    @static if VERSION > v"0.6.9"
        qa = q.args[2].args[3].args[3].args
    else
        qa = q.args[2].args[2].args[2].args
    end
    for i ∈ 2:N2
        U_i = Symbol(:U_, i)
        push!(qa, :($U_i = U.data[$i]))
    end

    for i ∈ N:-1:1
        lti = ltriangle(i)
        bti = lti+i
        Ui_i = Symbol(:Ui_, bti)
        U_i = Symbol(:U_, bti)
        push!(qa, :($Ui_i = $U_i))
        for j ∈ i+1:N
            Uj_i = Symbol(:Ui_, lti + j)
            push!(qa, :($Ui_i -= $Uj_i*Uj_i))
        end
        push!(qa, :($Ui_i = sqrt($Ui_i)))


        for j ∈ i-1:-1:1
            ltj = ltriangle(j)
            Uj_i = Symbol(:Ui_, ltj+i)
            U_i  = Symbol(:U_,  ltj+i)
            push!(qa, :($Uj_i = $U_i))
            for k = i+1:N
                push!(qa, :($Uj_i -= $Symbol(:Ui_, ltj+k) * $Symbol(:Ui_, lti+k)))
            end
            push!(qa, :($Uj_i /= $Ui_i) )
        end
    end

    push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    q

end

@generated function revchol(A::SymmetricMatrix{T,N,N2}) where {T,N,N2}
    gen_revchol_quote(UpperTriangularMatrix{T,N,N2},N,N2)
end
