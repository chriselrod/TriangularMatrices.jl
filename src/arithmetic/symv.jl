


@generated function LinearAlgebra.mul!(d::RecursiveVector{T,N}, S::SymmetricMatrix{T,N,L}, x::RecursiveVector{T,N}) where {T,N,L}
    q, qa = create_quote()
    extract_linear!(qa, L, :S)
    extract_linear!(qa, N, :x)

    dnr = Symbol(:d_, 1)
    push!(qa, :( $dnr = S[1] * x[1] ) )
    for c ∈ 2:N
        push!(qa, :( $dnr += S[ $(small_triangle(c)+1) ] * x[$c] ) )
    end

    for r ∈ 2:N
        dnr = Symbol(:d_, r)
        push!(qa, :( $dnr = S[ $(small_triangle(r)+1) ] * x[1] ) )
        for c ∈ 2:r-1
            push!(qa, :( $dnr += S[ $(small_triangle(r)+c) ] * x[$c] ) )
        end
        # push!(qa, :() )
        for c ∈ r:N
            push!(qa, :( $dnr += S[ $(small_triangle(c)+r) ] * x[$c] ) )
        end
    end
    insert_linear!(qa, N, :d)
    push!(q.args, :(d))
    q
end