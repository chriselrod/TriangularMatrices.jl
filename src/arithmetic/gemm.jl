
"""
A: M x N
B: N x P
C: M x P

C = A * B
"""

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
    q, qa = create_quote()
    # Because this gets called when max(M,N,P) <= 2cutoff, we go ahead and split all dims > cutoff
    if M > cutoff && N > cutoff && P > cutoff
        Mh, Ml = splitint(M)
        Nh, Nl = splitint(N)
        Ph, Pl = splitint(P)
        extract_linear!(qa, Mh*Nh, :A)#, ind_offset = 0, label_offset = 0) #A11
        extract_linear!(qa, Nh*Ph, :B)#, ind_offset = 0, label_offset = 0) #B11
        mul_quote!(qa, Mh, Nh, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, Mh*Nl, :A, ifelse(tA, Mh * Nh, M  * Nh) ) #A12
        extract_linear!(qa, Nl*Ph, :B, ifelse(tB, N  * Ph, Nh * Ph) ) #B21
        mul_quote!(qa, Mh, Nl, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, Mh * Ph, :C) # C11

        extract_linear!(qa, Ml*Nl, :A, M*N - Ml*Nl ) #A22
        mul_quote!(qa, Ml, Nl, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, Ml*Nh, :A, ifelse(tA, M  * Nh, Mh * Nh) )  #A21
        extract_linear!(qa, Nh*Ph, :B)                                 #B11
        mul_quote!(qa, Ml, Nh, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, Ml * Ph, :C, ifelse(tC, M  * Ph, Mh * Ph) ) #C21

        extract_linear!(qa, Nh*Pl, :B, ifelse(tB, Nh * Ph, N * Ph))   #B12
        mul_quote!(qa, Ml, Nh, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, Ml*Nl, :A, M*N - Ml*Nl ) #A22
        extract_linear!(qa, Nl*Pl, :B, N*P - Nl*Pl ) #B22
        mul_quote!(qa, Ml, Nl, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, Ml * Pl, :C, M*P - Ml*Pl ) #C22

        extract_linear!(qa, Mh*Nl, :A, ifelse(tA, Mh * Nh, M  * Nh) ) #A12
        mul_quote!(qa, Mh, Nl, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, Mh*Nh, :A ) #A22
        extract_linear!(qa, Nh*Pl, :B, ifelse(tB,  Nh * Ph, N * Ph)) #B12
        mul_quote!(qa, Mh, Nh, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, Mh * Pl, :C, ifelse(tC, Mh * Ph, M  * Ph) ) #C12

    elseif M > cutoff && N > cutoff # don't split P
        Mh, Ml = splitint(M)
        Nh, Nl = splitint(N)

        extract_linear!(qa, Mh*Nh, :A)#, ind_offset = 0, label_offset = 0) #A11
        extract_linear!(qa, Nh*P, :B)#, ind_offset = 0, label_offset = 0) #B11
        mul_quote!(qa, Mh, Nh, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, Mh*Nl, :A, ifelse(tA, Mh * Nh, M  * Nh) ) #A12
        extract_linear!(qa, Nl*P, :B, Nh * P ) #B21
        mul_quote!(qa, Mh, Nl, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, Mh * P, :C) # C11

        extract_linear!(qa, Ml*Nl, :A, M*N - Ml*Nl ) #A22
        mul_quote!(qa, Ml, Nl, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, Ml*Nh, :A, ifelse(tA, M  * Nh, Mh * Nh) )  #A21
        extract_linear!(qa, Nh*P, :B)                                 #B11
        mul_quote!(qa, Ml, Nh, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, Ml * Ph, :C, ifelse(tC, M  * Ph, Mh * Ph) ) #C21
    elseif M > cutoff && P > cutoff # don't split N
        Mh, Ml = splitint(M)
        Ph, Pl = splitint(P)

        extract_linear!(qa, Mh*N, :A)#, ind_offset = 0, label_offset = 0) #A1
        extract_linear!(qa, N*Ph, :B)#, ind_offset = 0, label_offset = 0) #B1
        mul_quote!(qa, Mh, N, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, Mh * Ph, :C) # C11

        extract_linear!(qa, Ml*N, :A, Mh * N ) #A2
        mul_quote!(qa, Ml, N, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, Ml * Ph, :C, ifelse(tC, M  * Ph, Mh * Ph) ) #C21

        extract_linear!(qa, N*Pl, :B, N * Ph ) #B2
        mul_quote!(qa, Ml, N, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, Ml * Pl, :C, M*P - Ml*Pl ) #C22

        extract_linear!(qa, Mh*N, :A) #A1
        mul_quote!(qa, Mh, N, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, Mh * Pl, :C, ifelse(tC, Mh * Ph, M  * Ph) ) #C12
    elseif N > cutoff && P > cutoff # don't split M
        Nh, Nl = splitint(N)
        Ph, Pl = splitint(P)

        extract_linear!(qa, M*Nh, :A)#, ind_offset = 0, label_offset = 0) #A1
        extract_linear!(qa, Nh*Ph, :B)#, ind_offset = 0, label_offset = 0) #B11
        mul_quote!(qa, M, Nh, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, M*Nl, :A, M * Nh ) #A2
        extract_linear!(qa, Nl*Ph, :B, ifelse(tB, N  * Ph, Nh * Ph) ) #B21
        mul_quote!(qa, M, Nl, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, M * Ph, :C) # C11

        extract_linear!(qa, Nl*Pl, :B, N*P - Nl*Pl ) #B22
        mul_quote!(qa, M, Nl, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        extract_linear!(qa, M*Nh, :A)
        extract_linear!(qa, Nh*Pl, :B, ifelse(tB, Nh * Ph, N * Ph)) #B12
        mul_quote!(qa, M, Nh, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, :(+=) )

        insert_linear!(qa, M * Pl, :C, M*Ph ) #C2


    elseif M > cutoff # only split M
        Mh, Ml = splitint(M)

        extract_linear!(qa, Mh*N, :A)#, ind_offset = 0, label_offset = 0) #A11
        extract_linear!(qa, N*P, :B)#, ind_offset = 0, label_offset = 0) #B11
        mul_quote!(qa, Mh, N, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, Mh * P, :C) # C11

        extract_linear!(qa, Ml*N, :A, Mh*N ) #A2
        mul_quote!(qa, Ml, N, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, Ml * P, :C, Mh * P ) #C2
    elseif N > cutoff # only split N
        Nh, Nl = splitint(N)
    else # only split P
        Ph, Pl = splitint(P)

        extract_linear!(qa, M*N, :A)#, ind_offset = 0, label_offset = 0) #A1
        extract_linear!(qa, N*Ph, :B)#, ind_offset = 0, label_offset = 0) #B11
        mul_quote!(qa, M, N, Ph, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, M * Ph, :C) # C11

        extract_linear!(qa, N*Pl, :B, N*Ph ) #B2
        mul_quote!(qa, M, N, Pl, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )

        insert_linear!(qa, M * Pl, :C, M*Ph ) #C2
    end
    q
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

function recursion_mul(M, N, P, tA=false, tB=false, tC=false, eq = :(=), LA = M*N, LB = N*P, LC = M*P)

    A_blocks = block_dim(M, N)
    B_blocks = block_dim(N, P)
    while size(A_blocks,2) < size(B_blocks,1) # should be equal, but under certain "pathological" cases they may not be.
        A_blocks = split_col(A_blocks)
    end
    while size(A_blocks,2) > size(B_blocks,1)
        B_blocks = split_row(B_blocks)
    end
    # The matrices should now line up.
    # Is it possible that they don't?
    
end


"""
The dummy argument always gets optimized out.
However, when using code_llvm on the function and filling it with a real type (eg, Bool), it causes the SSA names to get incremented by 1 (they start as %0 for the first argument, %1 for the second...). In llvmcall, these values get incremented by one between the first in the body and the argument list, while they do not in the function code returned by llvmcall.
Maybe there is a solution that is less of a hack, but just passing an extra dummy argument to the code_llvm call to get them to line up was a simple solution.
"""
@generated function mul!(C::RecursiveMatrixOrTranpose{T,M,P,LC},
                        A::RecursiveMatrixOrTranpose{T,M,N,LA},
                        B::RecursiveMatrixOrTranpose{T,N,P,LB},
                        dummy = true) where {T,M,N,P,LA,LB,LC}
                        #dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    maxdim = max(M,N,P)
    if maxdim <= cutoff # No recursion; we multiply.
        return mul_kernel(M,N,P, istransposed(A), istransposed(B), istransposed(C))
    elseif maxdim <= 2cutoff
        return block_kernel(M,N,P, istransposed(A), istransposed(B), istransposed(C))
    else # Recursion.
        return recursion_mul(M,N,P, istransposed(A), istransposed(B), istransposed(C))
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

const r64_1 = Ref{NTuple{64,Float64}}()
const r64_2 = Ref{NTuple{64,Float64}}()
const r64_3 = Ref{NTuple{64,Float64}}()

function create_mul_method(T,M,N,P,LA = M*N,LB = N*P,LC = M*P)
    iobuffer = IOBuffer()
    code_llvm(iobuffer, mul!, (RecursiveMatrix{T,M,P,LC}, RecursiveMatrix{T,M,N,LC}, RecursiveMatrix{T,N,P,LC}, Bool))
    kernel_code = String(iobuffer)
    codestart = search(kernel_code,"{\ntop:\n ")[end]
    @eval function mul!(C:: PointerRecursiveMatrix{$T,$M,$P,$LC},
                        A:: PointerRecursiveMatrix{$T,$M,$N,$LA},
                        B:: PointerRecursiveMatrix{$T,$N,$P,$LB})
            # Base.llvmcall($(kernel_code[codestart:end-3]), Base.RefValue{NTuple{$LC,$T}}, Tuple{Base.RefValue{NTuple{$LC,$T}},Base.RefValue{NTuple{$LA,$T}},Base.RefValue{NTuple{$LB,$T}}}, Base.unsafe_convert(Base.RefValue{NTuple{$LC,$T}}, Base.unsafe_convert(Ptr{NTuple{$LC,$T}}, C.data)), Base.unsafe_convert(Base.RefValue{NTuple{$LA,$T}}, Base.unsafe_convert(Ptr{NTuple{$LA,$T}}, A.data)), Base.unsafe_convert(Base.RefValue{NTuple{$LB,$T}}, Base.unsafe_convert(Ptr{NTuple{$LB,$T}}, B.data)))
            
            #works, but slow
            # Base.llvmcall($(kernel_code[codestart:end-3]), Base.RefValue{NTuple{$LC,$T}}, Tuple{Base.RefValue{NTuple{$LC,$T}},Base.RefValue{NTuple{$LA,$T}},Base.RefValue{NTuple{$LB,$T}}}, Ref(Base.unsafe_load(Base.unsafe_convert(Ptr{NTuple{$LC,$T}}, C.data),1)), Ref(Base.unsafe_load(Base.unsafe_convert(Ptr{NTuple{$LA,$T}}, A.data),1)), Ref(Base.unsafe_load(Base.unsafe_convert(Ptr{NTuple{$LB,$T}}, B.data),1)))

            #works, but also pretty slow
            # r64_1[] = Base.unsafe_load(Base.unsafe_convert(Ptr{NTuple{$LC,$T}}, C.data),1)
            # r64_2[] = Base.unsafe_load(Base.unsafe_convert(Ptr{NTuple{$LA,$T}}, A.data),1)
            # r64_3[] = Base.unsafe_load(Base.unsafe_convert(Ptr{NTuple{$LB,$T}}, B.data),1)
            
            # Base.llvmcall($(kernel_code[codestart:end-3]), Base.RefValue{NTuple{$LC,$T}}, Tuple{Base.RefValue{NTuple{$LC,$T}},Base.RefValue{NTuple{$LA,$T}},Base.RefValue{NTuple{$LB,$T}}}, r64_1, r64_2, r64_3)


            
            Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{Ptr{Int8}}, Tuple{Ptr{Ptr{Int8}},Ptr{Ptr{Int8}},Ptr{Ptr{Int8}}}, pointer(Base.unsafe_convert(Ptr{Int8}, C.data)), pointer(Base.unsafe_convert(Ptr{Int8}, A.data)), pointer(Base.unsafe_convert(Ptr{Int8}, B.data)))
            
            # Base.llvmcall($(kernel_code[codestart:end-3]), Base.RefValue{NTuple{$LC,$T}}, Tuple{Base.RefValue{NTuple{$LC,$T}},Base.RefValue{NTuple{$LA,$T}},Base.RefValue{NTuple{$LB,$T}}}, C.data, A.data, B.data)

            # Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{NTuple{$LC,$T}}, Tuple{Ptr{NTuple{$LC,$T}},Ptr{NTuple{$LA,$T}},Ptr{NTuple{$LB,$T}}}, Base.unsafe_convert(Ptr{NTuple{$LC,$T}}, C.data), Base.unsafe_convert(Ptr{NTuple{$LA,$T}}, A.data), Base.unsafe_convert(Ptr{NTuple{$LB,$T}}, B.data))
            # Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{$T}, Tuple{Ptr{$T},Ptr{$T},Ptr{$T}}, C.data, A.data, B.data)
            C
        end

end
function create_gemm_method(T,M,N,P,LA = M*N,LB = N*P,LC = M*P)
    iobuffer = IOBuffer()
    code_llvm(iobuffer, gemm!, (RecursiveMatrix{T,M,P,LC}, RecursiveMatrix{T,M,N,LC}, RecursiveMatrix{T,N,P,LC}, Bool))
    kernel_code = String(iobuffer)
    codestart = search(kernel_code,"{\ntop:\n ")[end]
    @eval function gemm!(C:: PointerRecursiveMatrix{$T,$M,$P,$LC},
                        A:: PointerRecursiveMatrix{$T,$M,$N,$LA},
                        B:: PointerRecursiveMatrix{$T,$N,$P,$LB})
            Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{NTuple{$LC,$T}}, Tuple{Ptr{NTuple{$LC,$T}},Ptr{NTuple{$LA,$T}},Ptr{NTuple{$LB,$T}}}, Base.unsafe_convert(NTuple{$LC,$T}, C.data), Base.unsafe_convert(NTuple{$LA,$T}, A.data), Base.unsafe_convert(NTuple{$LB,$T}, B.data))
            C
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