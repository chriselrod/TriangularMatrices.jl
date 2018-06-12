
function gen_inv_quote(T, N, N2)
    q, qa = initial_quote(N2::Int, :U)

    for i ∈ 1:N
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = Symbol(:Ui_, bti)
        U_i = Symbol(:U_, bti)
        push!(qa, :($Ui_i = 1 / $U_i))
        for j ∈ i+1:N
            ltj = small_triangle(j)
            Ui_i = Symbol(:Ui_, bti)
            U_i = Symbol(:U_, ltj+i)
            push!(qa, :(Utemp = $U_i * $Ui_i))
            for k ∈ i+1:j-1
                Ui_i = Symbol(:Ui_, small_triangle(k) + i)
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
function gen_ip_inv_quote!(qa, N, N2, Ui = :Ui, U = :U, extract = extract_symbol)

    temp_num = 0

    for i ∈ 1:N
        lti = small_triangle(i)
        bti = lti+i
        Ui_i = extract(Ui, bti)
        U_i = extract(U, bti)
        push!(qa, :($Ui_i = 1 / $U_i))
        for j ∈ i+1:N
            ltj = small_triangle(j)
            Ui_i = extract(Ui, bti)
            U_i = extract(U, ltj+i)
            Utemp = Symbol(:t_, temp_num)
            temp_num += 1 #Comment this line on / off to change.
            push!(qa, :($Utemp = $U_i * $Ui_i)) # I'm curious. Woud renaming Utemps help optimization, or just be wasteful?
            for k ∈ i+1:j-1
                Ui_i = extract(Ui, small_triangle(k) + i)
                U_i = extract(U, ltj+k)
                push!(qa, :($Utemp += $U_i * $Ui_i))
            end
            Ui_i = extract(Ui, ltj+i)
            U_i = extract(U, ltj+j)
            push!(qa, :($Ui_i = -$Utemp / $U_i))
        end
    end

    # push!(q.args, :($T( ( @ntuple $N2 Ui ) )))
    qa
end


@generated function LinearAlgebra.inv(U::AbstractUpperTriangular{T,N,N2}) where {T,N,N2}
    gen_inv_quote(UpperTriangularMatrix{T,N,N2}, N, N2)
end
@generated function LinearAlgebra.inv(U::AbstractLowerTriangular{T,N,N2}) where {T,N,N2}
    gen_inv_quote(LowerTriangularMatrix{T,N,N2}, N, N2)
end

function gen_extracted_ip_inv_quote!(qa, N, L, Ui = :Ui, U = :U)
    for i ∈ 1:L
        push!(qa, :($(id_symbol(U, i)) = $(U)[$i]))
    end
    gen_ip_inv_quote!(qa, N, L, Ui, U, id_symbol)
    for i ∈ 1:L
        push!(qa, :($(Ui)[$i] = $(id_symbol(Ui, i))))
    end
end


@generated function LinearAlgebra.inv!(U::Union{UpperTriangularMatrix{T,N,L},LowerTriangularMatrix{T,N,L}}) where {T,N,L}
    q, qa = create_quote()
    gen_extracted_ip_inv_quote!(qa, N, L, :U, :U)
    push!(q.args, :U)
    q
end
@generated function LinearAlgebra.inv!(Ui::UpperTriangularMatrix{T,N,L}, U::AbstractUpperTriangular{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    gen_extracted_ip_inv_quote!(qa, N, L, :Ui, :U)
    push!(q.args, :Ui)
    q
end
@generated function LinearAlgebra.inv!(Ui::LowerTriangularMatrix{T,N,L},U::AbstractLowerTriangular{T,N,L}) where {T,N,L}
    q, qa = create_quote()
    gen_extracted_ip_inv_quote!(qa, N, L, :Ui, :U)
    push!(q.args, :Ui)
    q
end