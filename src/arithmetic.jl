

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


# output is an N x N2 SMatrix
function AU_quote!(qa, N, N2)
    for i ∈ 1:N2
        lmi = (i-1)*N
        lti = ltriangle(i)
        for j ∈ 1:N
            C_i = Symbol(:C_, j + lmi)
            A_i = Symbol(:A_, j)
            U_i = Symbol(:U_, lti + 1)
            push!(qa, :($C_i = $A_i * $U_i) )
            for k ∈ 2:i
                A_i = Symbol(:A_, j + (k-1)*N)
                U_i = Symbol(:U_, lti + k)
                push!(qa, :($C_i += $A_i * $U_i) )
            end
        end
    end
end
function UA_quote!(qa, N, N2)
    for i ∈ 1:N2
        lmi = (i-1)*N
        # lti = ltriangle(i)
        for j ∈ 1:N
            # ltj = (j-1)*N
            C_i = Symbol(:C_, j + lmi)
            A_i = Symbol(:A_, j + lmi)
            U_i = Symbol(:U_, btriangle(j)) # equivalent to j + ltriangle(j)
            push!(qa, :($C_i = $A_i * $U_i) )
            for k ∈ j+1:N
                A_i = Symbol(:A_, k + lmi)
                U_i = Symbol(:U_, j + ltriangle(k))
                push!(qa, :($C_i += $A_i * $U_i) )
            end
        end
    end
end



# Shoot, L is row-major storate.
# Must fix, or write row major extractor? #Writing for SMatrices anyway
function LU_quote!(N)
    for i ∈ 1:N
        lmi = (i-1)*N
        lti = ltriangle(i)
        for j ∈ 1:N
            ltj = ltriangle(j)
            C_i = Symbol(:C_, j + lmi)
            L_i = Symbol(:L_, 1 + ltj)
            U_i = Symbol(:U_, 1 + lti)
            push!(qa, :($C_i = $L_i * $U_i) )
            for k ∈ 2:min(i,j)
                L_i = Symbol(:L_, k + ltj)
                U_i = Symbol(:U_, k + lti)
                push!(qa, :($C_i += $L_i * $U_i) )
            end
        end
    end
end
function UL_quote!(N)
    for i ∈ 1:N
        lmi = (i-1)*N
        lti = ltriangle(i)
        for j ∈ 1:N
            k_start = max(i,j)
            ltk = ltriangle(ltk)
            C_i = Symbol(:C_, j + lmi)
            U_i = Symbol(:U_, ltk + j)
            L_i = Symbol(:L_, ltk + i)
            push!(qa, :($C_i = $A1_i * $A2_i) )
            for k ∈ k_start+1:N
                ltk = ltriangle(k)
                U_i = Symbol(:U_, ltk + j)
                L_i = Symbol(:L_, ltk + i)
                push!(qa, :($C_i += $U_i * $L_i) )
            end
        end
    end
end
function AS_quote!(N,N2)
    for i ∈ 1:N
        lmi = (i-1)*N2
        lti = ltriangle(i)
        for j ∈ 1:N2
            C_i = Symbol(:C_, j + lmi)
            A_i = Symbol(:A_, j)
            S_i = Symbol(:S_, 1 + lti)
            push!(qa, :($C_i = $A_i * $S_i) )
            for k ∈ 2:i #we're in upper triangle of S
                A_i = Symbol(:A_, j + (k-1)*N2)
                S_i = Symbol(:S_, k + lti)
                push!(qa, :($C_i += $A_i * $S_i) )
            end
            for k ∈ 1+i:N
                #we're in mirrored lower triangle
                #thus, k is column, i is row
                A_i = Symbol(:A_, j + (k-1)*N2)
                S_i = Symbol(:S_, i + ltriangle(k))
                push!(qa, :($C_i += $A_i * $S_i) )
            end
        end
    end
end
function US_quote!()

end
function SU_quote!()

end
function USU_quote!()

end
function LSL_quote!()

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
