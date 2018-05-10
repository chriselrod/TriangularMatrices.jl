
function gen_inv_quote(T, N, N2)
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
        push!(qa, :($Ui_i = 1 / $U_i))
        for j ∈ i+1:N
            ltj = ltriangle(j)
            Ui_i = Symbol(:Ui_, bti)
            U_i = Symbol(:U_, ltj+i)
            push!(qa, :(Utemp = $U_i * $Ui_i))
            for k ∈ i+1:j-1
                Ui_i = Symbol(:Ui_, ltriangle(k) + i)
                U_i = Symbol(:U_, ltj+k)
                push!(qa, :(Utemp += $U_i * $Ui_i))
            end
            Ui_i = Symbol(:Ui_, ltj+i)
            U_i = Symbol(:U_, ltj+j)
            push!(qa, :($Ui_i = -Utemp / $U_i))
        end
    end

    push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    q
end


@generated function LinearAlgebra.inv(U::UpperTriangularMatrix{T,N,N2}) where {T,N,N2}
    gen_inv_quote(UpperTriangularMatrix{T,N,N2}, N, N2)
end
@generated function LinearAlgebra.inv(U::LowerTriangularMatrix{T,N,N2}) where {T,N,N2}
    gen_inv_quote(LowerTriangularMatrix{T,N,N2}, N, N2)
end