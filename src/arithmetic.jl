

function gen_upper_mul_quote(T, N, N2)
    q, qa = initial_quote(N2)
    for i ∈ 1:N2
        push!(qa, :($B_i = B.data[$i]))
    end
    for i ∈ 1:N
        lti = ltriangle(i)
        for j ∈ 1:i
            C_i = Symbol(:C_, lti + j)
            B_i = Symbol(:B_, lti + j)
            A_i = Symbol(:A_, ltriangle(j) + j)
            push!(qa, :($C_i = $A_i * $B_i) )
            for k ∈ j+1:i
                A_i = Symbol(:A_, ltriangle(k) + j)
                B_i = Symbol(:B_, lti + k)
                push!(qa, :($C_i += $A_i * $B_i) )
            end
        end
    end
    push!(q.args, :($T( ( @ntuple $N2 C ) )))
end

@generated function Base.:*(A::UpperTriangularMatrix{T,N,N2}, B::UpperTriangularMatrix{T,N,N2}) where {T,N,N2}
    gen_upper_mul_quote(UpperTriangularMatrix{T,N,N2}, N, N2)
end



function gen_lower_mul_quote(T, N, N2)
    q, qa = initial_quote(N2)
    for i ∈ 1:N2
        push!(qa, :($B_i = B.data[$i]))
    end
    for i ∈ 1:N
        lti = ltriangle(i)
        for j ∈ i:N
            ltj = ltriangle(j)
            C_i = Symbol(:C_, ltj + i)
            B_i = Symbol(:B_, lti + i)
            A_i = Symbol(:A_, ltj + i)
            push!(qa, :($C_i = $A_i * $B_i) )
            for k ∈ i+1:j
                A_i = Symbol(:A_, ltj + k)
                B_i = Symbol(:B_, ltriangle(k) + i)
                push!(qa, :($C_i += $A_i * $B_i) )
            end
        end
    end
    push!(q.args, :($T( ( @ntuple $N2 C ) )))
end

@generated function Base.:*(A::LowerTriangularMatrix{T,N,N2}, B::LowerTriangularMatrix{T,N,N2}) where {T,N,N2}
    gen_lower_mul_quote(LowerTriangularMatrix{T,N,N2}, N, N2)
end



function gen_upper_xxt_quote(T,N,N2)
    q, qa = initial_quote(N2)
    for i ∈ 1:N
        lti = ltriangle(i)
        for j ∈ 1:i
            B_i = Symbol(:B_, lti + j)
            A1_i = Symbol(:A_, lti + j)
            A2_i = Symbol(:A_, lti + i)
            push!(qa, :($B_i = $A1_i * $A2_i) )
            for k ∈ i+1:N
                ltk = ltriangle(k)
                A1_i = Symbol(:A_, ltk + j)
                A2_i = Symbol(:A_, ltk + i)
                push!(qa, :($B_i += $A1_i * $A2_i) )
            end
        end
    end
    push!(q.args, :(SymmetricMatrix{T,N,N2}( ( @ntuple $N2 B ) )))
end

function gen_upper_xtx_quote(T,N,N2)
    q, qa = initial_quote(N2)
    for i ∈ 1:N
        lti = ltriangle(i)
        for j ∈ 1:i
            ltj = ltriangle(j)
            B_i = Symbol(:B_, lti + j)
            A1_i = Symbol(:A_, ltj + 1)
            A2_i = Symbol(:A_, lti + 1)
            push!(qa, :($B_i = $A1_i * $A2_i) )
            for k ∈ 2:j
                A1_i = Symbol(:A_, ltj + k)
                A2_i = Symbol(:A_, lti + k)
                push!(qa, :($B_i += $A1_i * $A2_i) )
            end
        end
    end
    push!(q.args, :(SymmetricMatrix{T,N,N2}( ( @ntuple $N2 B ) )))
end


function gen_xxt_quote(T,N,N2)
    q, qa = initial_quote(N*N2)
    for i ∈ 1:N2
        lti = (i-1)*N
        for j ∈ 1:i
            B_i = Symbol(:B_, lti + j)
            A1_i = Symbol(:A_, j)
            A2_i = Symbol(:A_, i)
            push!(qa, :($B_i = $A1_i * $A2_i) )
            for k ∈ 1:N-1
                Nk = N*k
                A1_i = Symbol(:A_, Nk + j)
                A2_i = Symbol(:A_, Nk + i)
                push!(qa, :($B_i += $A1_i * $A2_i) )
            end
        end
    end
    push!(q.args, :(SymmetricMatrix{T,N,N2}( ( @ntuple $N2 B ) )))
end
function gen_xtx_quote(T,N,N2)
    q, qa = initial_quote(N*N2)
    for i ∈ 1:N2
        lti = (i-1)*N
        for j ∈ 1:i
            ltj = (j-1)*N
            B_i = Symbol(:B_, lti + j)
            A1_i = Symbol(:A_, ltj + 1)
            A2_i = Symbol(:A_, lti + 1)
            push!(qa, :($B_i = $A1_i * $A2_i) )
            for k ∈ 2:N
                A1_i = Symbol(:A_, ltj + k)
                A2_i = Symbol(:A_, lti + k)
                push!(qa, :($B_i += $A1_i * $A2_i) )
            end
        end
    end
    push!(q.args, :(SymmetricMatrix{T,N,N2}( ( @ntuple $N2 B ) )))
end



@generated function xxt(A::UpperTriangularMatrix{T,N,N2}) where {T,N,N2}
    gen_upper_xxt_quote(T,N,N2)
end
@generated function xtx(A::UpperTriangularMatrix{T,N,N2}) where {T,N,N2}
    gen_upper_xtx_quote(T,N,N2)
end
@generated function xxt(A::SMatrix{N,N2,T}) where {T,N,N2}
    gen_xxt_quote(T,N,N2)
end
@generated function xtx(A::SMatrix{N,N2,T}) where {T,N,N2}
    gen_xxt_quote(T,N,N2)
end

xxt(A::LowerTriangularMatrix{T,N,N2}) where {T,N,N2} = xtx(UpperTriangularMatrix{T,N,N2}(A.data))
xtx(A::LowerTriangularMatrix{T,N,N2}) where {T,N,N2} = xxt(UpperTriangularMatrix{T,N,N2}(A.data))


function gen_AL_quote()

end


@generated function Base.:*(A::SMatrix{N3,N}, B::LowerTriangularMatrix{T,N,N2}) where {T,N,N2,N3}
    gen_lower_mul_quote(LowerTriangularMatrix{T,N,N2}, N, N2)
end
@generated function Base.:*(A::SMatrix{N3,N}, B::UpperTriangularMatrix{T,N,N2}) where {T,N,N2,N3}
    gen_lower_mul_quote(LowerTriangularMatrix{T,N,N2}, N, N2)
end
@generated function Base.:*(A::LowerTriangularMatrix{T,N,N2}, B::SMatrix{N3,N}) where {T,N,N2,N3}
    gen_lower_mul_quote(LowerTriangularMatrix{T,N,N2}, N, N2)
end
@generated function Base.:*(A::UpperTriangularMatrix{T,N,N2}, B::SMatrix{N3,N}) where {T,N,N2,N3}
    gen_lower_mul_quote(LowerTriangularMatrix{T,N,N2}, N, N2)
end
