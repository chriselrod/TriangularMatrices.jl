
"""
A: M x N
B: N x P
C: M x P

C = A * B
"""
function broadcast_mul_quote!(qa, ::Type{T}, M, N, P,
                ta::Val{false}=Val(false), A = :A, tb::Val{false}=Val(false), B = :B, C = :C,
                extract = extract_symbol, insert = extract_symbol, equalop = :(=) ) where T

    eA = (i,j) -> extract(A, sub2ind(false, (M, N), i, j))
    eB = (i,j) -> extract(B, sub2ind(false, (N, P), i, j))
    eC = (i,j) -> insert(C, sub2ind(false, (M, P), i, j))
    chunk = max(1,32 ÷ sizeof(T))
    Num_chunks, Num_remaining = divrem(M, chunk)
    total_chunked = chunk * Num_chunks # = M - Num_remaining
    # The goal is to broadcast B 
    # The remainers will be done after, so we group similar computations together
    # Prioritize that over memory access order?
    # This is a kernel for small-scale computations, so I think that makes sense.

    
    
    # Handling chunked and non-chunked separately, in hopes of encouraging vectorization...
    for p = 1:P
        B_np = eB(1,p)
        # for c = 0:Num_chunks-1
        #     for m = (1:chunk) + c*chunk
        #         push!(qa, Expr(equalop, eC(m, p), :( $B_np * $(eA(m,1)) ) ) )
        #     end
        # end
        for m = 1:total_chunked
            push!(qa, Expr( equalop, eC(m, p), :( $B_np * $(eA(m,1)) ) ) )
        end
        for n = 2:N
            B_np = eB(n,p)
            # for c = 0:Num_chunks-1
            #     for m = (1:chunk) + c*chunk
            #         push!(qa, Expr( :(+=), eC(m, p), :( $B_np * $(eA(m,n)) ) ) )
            #     end
            # end
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eC(m, p), :( $B_np * $(eA(m,n)) ) ) )
            end
        end
    end

    if Num_remaining > 0
        # Here, we follow a different pattern, trying to vectorize across N instead of across M

        N4, Nr = divrem(N, chunk)
        if N4 > 0
            for p = 1:P, m = 1+total_chunked:M
                C_mp = eC(m, p)
                push!(qa, Expr(equalop, C_mp, :(+$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1:chunk]...) ) ))
                for k = 2:N4
                    push!(qa, :($C_mp += +$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1+chunk*(k-1):chunk*k ]...) ) )
                end
                Nr > 0 && push!(qa, :($C_mp += +$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1+chunk*N4:N ]...) ) )
            end
        else
            for p = 1:P, m = 1+total_chunked:M
                C_mp = eC(m, p)
                push!(qa, Expr(equalop, C_mp, :( +$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1:Nr ]...) ) ))
            end
        end

        # Approach used for the first block.
        # for p = 1:P, n = 1:N
        #     B_np = eB(n,p)
        #     for m = 1+total_chunked:M
        #         push!(qa, Expr(equalop, eC(m, p), :( $B_np * $(eA(m,n)) ) ) )
        #     end
        # end
    end
end
function broadcast_mul_lazyextract_quote!(qa, ::Type{T}, M, N, P,
                ta::Val{false}=Val(false), A = :A, tb::Val{false}=Val(false), B = :B, C = :C,
                extract = extract_symbol ) where T

    eA = (i,j) -> extract(A, sub2ind(false, (M, N), i, j))
    eB = (i,j) -> extract(B, sub2ind(false, (N, P), i, j))
    exC = (i,j) -> extract_symbol(C, sub2ind(false, (M, P), i, j))
    eC = (i,j) -> id_symbol(C, sub2ind(false, (M, P), i, j))
    chunk = max(1,32 ÷ sizeof(T))
    Num_chunks, Num_remaining = divrem(M, chunk)
    total_chunked = chunk * Num_chunks # = M - Num_remaining
    # The goal is to broadcast B 
    # The remainers will be done after, so we group similar computations together
    # Prioritize that over memory access order?
    # This is a kernel for small-scale computations, so I think that makes sense.

    
    
    # Handling chunked and non-chunked separately, in hopes of encouraging vectorization...
    for p = 1:P
        B_np = eB(1,p)
        # for c = 0:Num_chunks-1
        #     for m = (1:chunk) + c*chunk
        #         push!(qa, Expr(equalop, eC(m, p), :( $B_np * $(eA(m,1)) ) ) )
        #     end
        # end
        for m = 1:total_chunked
            push!(qa, :( $(eC(m, p)) = $(exC(m, p)) + $B_np * $(eA(m,1)) ) )
        end
        for n = 2:N
            B_np = eB(n,p)
            # for c = 0:Num_chunks-1
            #     for m = (1:chunk) + c*chunk
            #         push!(qa, Expr( :(+=), eC(m, p), :( $B_np * $(eA(m,n)) ) ) )
            #     end
            # end
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eC(m, p), :( $B_np * $(eA(m,n)) ) ) )
            end
        end
    end

    if Num_remaining > 0
        # Here, we follow a different pattern, trying to vectorize across N instead of across M

        N4, Nr = divrem(N, chunk)
        if N4 > 0
            for p = 1:P, m = 1+total_chunked:M
                C_mp = eC(m, p)
                push!(qa, :( $C_mp = $(exC(m, p)) + +$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1:chunk]...) ) )
                for k = 2:N4
                    push!(qa, :($C_mp += +$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1+chunk*(k-1):chunk*k ]...) ) )
                end
                Nr > 0 && push!(qa, :($C_mp += +$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1+chunk*N4:N ]...) ) )
            end
        else
            for p = 1:P, m = 1+total_chunked:M
                push!(qa, :( $(eC(m, p)) = $(exC(m, p)) +$([:( $(eA(m,n)) * $(eB(n,p))) for n ∈ 1:Nr ]...) ) )
            end
        end

    end
end
function AX_plus_BY_quote!(qa, ::Type{T}, M, Nax, Nby, P,
                ta::Val{false}=Val(false), A = :A, tx::Val{false}=Val(false), X = :X, B = :B, Y = :Y, C = :C,
                extract = extract_symbol, insert = extract_symbol, equalop = :(=) ) where T

    eA = (i,j) -> extract(A, sub2ind(false, (M, Nax), i, j))
    eB = (i,j) -> extract(B, sub2ind(false, (M, Nby), i, j))
    eX = (i,j) -> extract(X, sub2ind(false, (Nax, P), i, j))
    eY = (i,j) -> extract(Y, sub2ind(false, (Nby, P), i, j))
    eC = (i,j) -> insert(C, sub2ind(false, (M, P), i, j))
    chunk = max(1,32 ÷ sizeof(T))
    Num_chunks, Num_remaining = divrem(M, chunk)
    total_chunked = chunk * Num_chunks # = M - Num_remaining
    # The goal is to broadcast B 
    # The remainers will be done after, so we group similar computations together
    # Prioritize that over memory access order?
    # This is a kernel for small-scale computations, so I think that makes sense.

    
    # The intention of this generated code is to
    # broadcast X and Y across a ymm register (or zmm on Skylake!?! -- may have to change chunk; would be good to test)
    # and leave these and aC_mp in a register for as long as possible, while repeatedly FMAing it.
    # we loop through the matrices A and B as doing we do it.

    # I'm assuming we're on a computer with 16 ymm registers.
    # That covers basically all modern X86 systems.
    #
    # For each p in 1:P,
    # X and Y each take up 1 register at a time.
    # we have Num_chunks worth of registers of C. Should be 2.
    # This leaves 12 registers devoted to A and B.
    #
    # Because we cycle through A and B per column,
    # Perhaps we should iteratively update C, choosing to reload just those 2 registers each time
    # As is, we are constantly swapping them out each time anyway.
    #
    # so we can reuse A and B?
    # Assuming 64 elements in A and B, they would require 16 registers each.
    # This means we could try breaking it into 3 separate loops, so we dedicate 11,11,10 registers to the both of them per loop.
    #
    for p = 1:P
        X_np = eX(1,p)
        Y_np = eY(1,p)
        for c = 1:Num_chunks, m = 1+(c-1)*chunk:c*chunk
            C_mp = eC(m, p)
            push!(qa, Expr( equalop, C_mp, :( $X_np * $(eA(m,1)) ) ) ) # we might be overwriting C the first time.
            push!(qa, Expr( :(+=),   C_mp, :( $Y_np * $(eB(m,1)) ) ) )
        end
        for n = 2:min(Nax,Nby)
            X_np = eX(n,p)
            Y_np = eY(n,p)
            for c = 1:Num_chunks, m = 1+(c-1)*chunk:c*chunk
                C_mp = eC(m, p)
                push!(qa, Expr( :(+=), C_mp, :( $X_np * $(eA(m,n)) ) ) )
                push!(qa, Expr( :(+=), C_mp, :( $Y_np * $(eB(m,n)) ) ) )
            end
        end
        # We ran out of rows for at least one of X and Y. If we only ran out of rows for one of them,
        # lets keep going with the other.
        if Nax != Nby
            if Nax > Nby
                for n = 1+min(Nax,Nby):Nax
                    X_np = eX(n,p)
                    for m = 1:total_chunked
                        push!(qa, Expr( :(+=), eC(m, p), :( $X_np * $(eA(m,n)) ) ) )
                    end
                end
            else
                for n = 1+min(Nax,Nby):Nby
                    Y_np = eY(n,p)
                    for m = 1:total_chunked
                        push!(qa, Expr( :(+=), eC(m, p), :( $Y_np * $(eB(m,n)) ) ) )
                    end
                end
            end
        end
    end

    # Now, after all the painstaking work, if things don't divide evenly we just
    # hope for the best on the rest of it.
    if Num_remaining > 0
        # Here, we follow a different pattern, trying to vectorize across N instead of across M

        N4ax, Nrax = divrem(Nax, chunk)
        N4by, Nrby = divrem(Nby, chunk)
        for p = 1:P, m = 1+total_chunked:M
            C_mp = eC(m, p)
            if N4ax > 0
                push!(qa, Expr(equalop, C_mp, :(+$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1:chunk]...) ) ))
                for k = 2:N4ax
                    push!(qa, :($C_mp += +$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*(k-1):chunk*k ]...) ) )
                end
                # Perhaps, try and combine excesses?
                if Nrax > 0 && N4by == 0
                    Nrax > 0 && push!(qa, :($C_mp += +$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ]...) ) )
                elseif Nrax > 0
                    push!(qa, Expr( :(+=) , C_mp, :(+$(vcat(
                        [:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ],
                        [:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]
                        )...) ) ) )
                end
            elseif N4by == 0
                push!(qa, Expr(equalop, C_mp, :(+$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ]...) ) ) )
            else
                push!(qa, Expr(equalop, C_mp, :(+$(vcat(
                    [:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ],
                    [:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]
                    )...) ) ) )
            end
            if N4by > 0
                push!(qa, :($C_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1:chunk]...) ) )
                for k = 2:N4by
                    push!(qa, :($C_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*(k-1):chunk*k ]...) ) )
                end
                if Nrby > 0 && Nrax == 0
                    push!(qa, :($C_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]...) ) )
                end
            elseif Nrax == 0
                push!(qa, :($C_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]...) ) ) 
            end
        end
        

        # Approach used for the first block.
        # for p = 1:P, n = 1:N
        #     B_np = eB(n,p)
        #     for m = 1+total_chunked:M
        #         push!(qa, Expr(equalop, eC(m, p), :( $B_np * $(eA(m,n)) ) ) )
        #     end
        # end
    end
end


### Seems slightly faster than AX_plus_BY2 on 0.7, but slower on 0.6
@generated function AX_plus_BY(A::StaticRecursiveMatrix{T,M,Nax,LA}, X::StaticRecursiveMatrix{T,Nax,P,LX},
                                B::StaticRecursiveMatrix{T,M,Nby,LB}, Y::StaticRecursiveMatrix{T,Nby,P,LY}) where {T,M,Nax,Nby,P,LA,LX,LB,LY}
    q, qa = create_quote(true)
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LX, :X)
    extract_linear!(qa, LB, :B)
    extract_linear!(qa, LY, :Y)
    AX_plus_BY_quote!(qa, T, M, Nax, Nby, P, Val(false), :A, Val(false), :X, :B, :Y, :D, id_symbol, id_symbol, :(=) )
    # insert_linear!(qa, LC, :C)
    # push!(q.args, :C)
    LD = M*P
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LD}( Base.Cartesian.@ntuple $LD D ) ))
    q

end

@generated function AX_plus_BY_plus_D(D::StaticRecursiveMatrix{T,M,P,LD},
                                A::StaticRecursiveMatrix{T,M,Nax,LA}, X::StaticRecursiveMatrix{T,Nax,P,LX},
                                B::StaticRecursiveMatrix{T,M,Nby,LB}, Y::StaticRecursiveMatrix{T,Nby,P,LY}) where {T,M,Nax,Nby,P,LA,LX,LB,LY,LD}
    q, qa = create_quote(true)
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LX, :X)
    extract_linear!(qa, LB, :B)
    extract_linear!(qa, LY, :Y)
    extract_linear!(qa, LD, :D)
    AX_plus_BY_quote!(qa, T, M, Nax, Nby, P, Val(false), :A, Val(false), :X, :B, :Y, :D, id_symbol, id_symbol, :(+=) )
    # insert_linear!(qa, LC, :C)
    # push!(q.args, :C)
    
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LD}( Base.Cartesian.@ntuple $LD D ) ))
    q

end

# While faster on 0.6, I will focus on 0.7.
# @generated function AX_plus_BY2(A::StaticRecursiveMatrix{T,M,Nax,LA}, X::StaticRecursiveMatrix{T,Nax,P,LX},
#                                 B::StaticRecursiveMatrix{T,M,Nby,LB}, Y::StaticRecursiveMatrix{T,Nby,P,LY}) where {T,M,Nax,Nby,P,LA,LX,LB,LY}
#     q, qa = create_quote(true)
#     @static if VERSION > v"0.7-"
#         pushfirst!(q.args, :(Base.@_inline_meta))
#     else
#         insert!(q.args, 1, :(Base.@_inline_meta))
#     end
#     extract_linear!(qa, LA, :A) # second arg is the name of the variable we're extracting.
#     extract_linear!(qa, LB, :A, :B, 0, LA)
#     for p = 1:P
#         extract_linear!(qa, Nax, :X, :X, (p-1)*Nax, (p-1)*(Nax+Nby)      )
#         extract_linear!(qa, Nby, :X, :Y, (p-1)*Nby, (p-1)*(Nax+Nby) + Nax)
#     end
#     broadcast_mul_quote!(qa, T, M, Nax + Nby, P, Val(false), :A, Val(false), :X, Val(false), :C, id_symbol, id_symbol, :(=) )
#     LC = M*P
#     # insert_linear!(qa, LC, :C)
    
#     # push!(q.args, :C)
#     # q
#     push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LC}( Base.Cartesian.@ntuple $LC C ) ))
#     q
# end

@generated function AX_plus_BY_plus_CZ(A::StaticRecursiveMatrix{T,M,Nax,LA}, X::StaticRecursiveMatrix{T,Nax,P,LX},
                                B::StaticRecursiveMatrix{T,M,Nby,LB}, Y::StaticRecursiveMatrix{T,Nby,P,LY},
                                C::StaticRecursiveMatrix{T,M,Ncz,LC}, Z::StaticRecursiveMatrix{T,Ncz,P,LZ}
                                ) where {T,M,Nax,Nby,Ncz,P,LA,LX,LB,LY,LC,LZ}
    q, qa = create_quote(true)
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A) # second arg is the name of the variable we're extracting.
    extract_linear!(qa, LB, :A, :B, 0, LA)
    extract_linear!(qa, LB, :A, :C, 0, LA + LB)
    N = Nax + Nby + Ncz
    for p = 1:P
        extract_linear!(qa, Nax, :X, :X, (p-1)*Nax, (p-1)*N             )
        extract_linear!(qa, Nby, :X, :Y, (p-1)*Nby, (p-1)*N + Nax       )
        extract_linear!(qa, Nby, :X, :Z, (p-1)*Ncz, (p-1)*N + Nax + Nby )
    end
    broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :X, :D, id_symbol, id_symbol, :(=) )
    LD = M*P
    # insert_linear!(qa, LD, :D)
    
    # push!(q.args, :D)
    # q
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LD}( Base.Cartesian.@ntuple $LD D ) ))
    q
end
@generated function AX_plus_BY_plus_CZ_plus_D(D::StaticRecursiveMatrix{T,M,P,LD},
                                A::StaticRecursiveMatrix{T,M,Nax,LA}, X::StaticRecursiveMatrix{T,Nax,P,LX},
                                B::StaticRecursiveMatrix{T,M,Nby,LB}, Y::StaticRecursiveMatrix{T,Nby,P,LY},
                                C::StaticRecursiveMatrix{T,M,Ncz,LC}, Z::StaticRecursiveMatrix{T,Ncz,P,LZ}
                                ) where {T,M,Nax,Nby,Ncz,P,LA,LX,LB,LY,LC,LZ,LD}
    q, qa = create_quote(true)
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A) # second arg is the name of the variable we're extracting.
    extract_linear!(qa, LB, :A, :B, 0, LA)
    extract_linear!(qa, LB, :A, :C, 0, LA + LB)
    N = Nax + Nby + Ncz
    for p = 1:P
        extract_linear!(qa, Nax, :X, :X, (p-1)*Nax, (p-1)*N             )
        extract_linear!(qa, Nby, :X, :Y, (p-1)*Nby, (p-1)*N + Nax       )
        extract_linear!(qa, Nby, :X, :Z, (p-1)*Ncz, (p-1)*N + Nax + Nby )
    end
    extract_linear!(qa, LD, :D)
    broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :X, :D, id_symbol, id_symbol, :(+=) )
    # insert_linear!(qa, LD, :D)
    
    # push!(q.args, :D)
    # q
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LD}( Base.Cartesian.@ntuple $LD D ) ))
    q
end
@generated function AX_plus_BY_plus_CZ_plus_Dv2(D::StaticRecursiveMatrix{T,M,P,LD},
                                       A::StaticRecursiveMatrix{T,M,Nax,LA}, X::StaticRecursiveMatrix{T,Nax,P,LX},
                                       B::StaticRecursiveMatrix{T,M,Nby,LB}, Y::StaticRecursiveMatrix{T,Nby,P,LY},
                                       C::StaticRecursiveMatrix{T,M,Ncz,LC}, Z::StaticRecursiveMatrix{T,Ncz,P,LZ}
                                       ) where {T,M,Nax,Nby,Ncz,P,LA,LX,LB,LY,LC,LZ,LD}
    q, qa = create_quote(true)
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A) # second arg is the name of the variable we're extracting.
    extract_linear!(qa, LB, :A, :B, 0, LA)
    extract_linear!(qa, LB, :A, :C, 0, LA + LB)
    N = Nax + Nby + Ncz
    for p = 1:P
        extract_linear!(qa, Nax, :X, :X, (p-1)*Nax, (p-1)*N             )
        extract_linear!(qa, Nby, :X, :Y, (p-1)*Nby, (p-1)*N + Nax       )
        extract_linear!(qa, Nby, :X, :Z, (p-1)*Ncz, (p-1)*N + Nax + Nby )
    end
    #extract_linear!(qa, LD, :D)
    broadcast_mul_lazyextract_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :X, :D, id_symbol )
    # insert_linear!(qa, LD, :D)
    
    # push!(q.args, :D)
    # q
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LD}( Base.Cartesian.@ntuple $LD D ) ))
    q
end

@generated function Base.:*(A::StaticRecursiveMatrix{T,M,N,LA}, X::StaticRecursiveMatrix{T,N,P,LX}) where {T,M,N,P,LA,LX}
    q, qa = create_quote()
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LX, :X)
    broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :X, :D, id_symbol, id_symbol, :(=) )
    # insert_linear!(qa, LC, :C)
    
    # push!(q.args, :C)
    # q
    LD = M*P
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LD}( Base.Cartesian.@ntuple $LD D ) ))
    q
end


@generated function AX_plus_D(D::StaticRecursiveMatrix{T,M,P,LD}, A::StaticRecursiveMatrix{T,M,N,LA}, X::StaticRecursiveMatrix{T,N,P,LX}) where {T,M,N,P,LA,LX,LD}
    q, qa = create_quote()
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LX, :X)
    extract_linear!(qa, LD, :D)
    broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :X, :D, id_symbol, id_symbol, :(+=) )
    # insert_linear!(qa, LC, :C)
    
    # push!(q.args, :C)
    # q
    # LC = M*P
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LC}( Base.Cartesian.@ntuple $LD D ) ))
    q
end

# @generated function muladd(A::StaticRecursiveMatrix{T,M,N,LA},
#                             B::StaticRecursiveMatrix{T,N,P,LB},
#                             C::StaticRecursiveMatrix{T,M,P,LC}) where {T,M,N,P,LA,LB,LC}

#     q, qa = create_quote()
#     @static if VERSION > v"0.7-"
#         pushfirst!(q.args, :(Base.@_inline_meta))
#     else
#         insert!(q.args, 1, :(Base.@_inline_meta))
#     end
#     extract_linear!(qa, LA, :A)
#     extract_linear!(qa, LB, :B)
#     # extract_linear!(qa, LC, :C)
#     broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :B, Val(false), :C, id_symbol, extract_symbol, :(+=) )
#     # insert_linear!(qa, LC, :C)
    
#     # push!(q.args, :C)
#     # q
#     push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LC}( Base.Cartesian.@ntuple $LC C ) ))
#     q
# end


function AX_plus_BY_plus_CZ_quote!(qa, ::Type{T}, M, Nax, Nby, Ndz, P,
                ta::Val{false}=Val(false), A = :A, tx::Val{false}=Val(false), X = :X, B = :B, Y = :Y, C = :C, Z = :Z, D = :D,
                extract = extract_symbol, insert = extract_symbol, equalop = :(=) ) where T

    eA = (i,j) -> extract(A, sub2ind(false, (M, Nax), i, j))
    eB = (i,j) -> extract(B, sub2ind(false, (M, Nby), i, j))
    eC = (i,j) -> extract(C, sub2ind(false, (M, Ndz), i, j))
    eX = (i,j) -> extract(X, sub2ind(false, (Nax, P), i, j))
    eY = (i,j) -> extract(Y, sub2ind(false, (Nby, P), i, j))
    eZ = (i,j) -> extract(Z, sub2ind(false, (Ndz, P), i, j))
    eD = (i,j) -> insert(D, sub2ind(false, (M, P), i, j))
    chunk = max(1,32 ÷ sizeof(T))
    Num_chunks, Num_remaining = divrem(M, chunk)
    total_chunked = chunk * Num_chunks # = M - Num_remaining
    # The goal is to broadcast B 
    # The remainers will be done after, so we group similar computations together
    # Prioritize that over memory access order?
    # This is a kernel for small-scale computations, so I think that makes sense.

    
    
    # Handling chunked and non-chunked separately, in hopes of encouraging vectorization...
    for p = 1:P
        XYZ_np = eX(1,p)
        for m = 1:total_chunked
            push!(qa, Expr( equalop, eD(m, p), :( $XYZ_np * $(eA(m,1)) ) ) )
        end
        for n = 2:Nax
            XYZ_np = eX(n,p)
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eD(m, p), :( $XYZ_np * $(eA(m,n)) ) ) )
            end
        end
        for n = 1:Nby
            XYZ_np = eY(n,p)
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eD(m, p), :( $XYZ_np * $(eB(m,n)) ) ) )
            end
        end
        for n = 1:Ndz
            XYZ_np = eZ(n,p)
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eD(m, p), :( $XYZ_np * $(eC(m,n)) ) ) )
            end
        end
    end

    if Num_remaining > 0
        # Here, we follow a different pattern, trying to vectorize across N instead of across M

        N4ax, Nrax = divrem(Nax, chunk)
        N4by, Nrby = divrem(Nby, chunk)
        for p = 1:P, m = 1+total_chunked:M
            D_mp = eD(m, p)
            if N4ax > 0
                push!(qa, Expr(equalop, D_mp, :(+$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1:chunk]...) ) ))
                for k = 2:N4ax
                    push!(qa, :($D_mp += +$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*(k-1):chunk*k ]...) ) )
                end
                # Perhaps, try and combine excesses?
                if Nrax > 0 && N4by == 0
                    Nrax > 0 && push!(qa, :($D_mp += +$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ]...) ) )
                elseif Nrax > 0
                    push!(qa, Expr( :(+=) , D_mp, :(+$(vcat(
                        [:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ],
                        [:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]
                        )...) ) ) )
                end
            elseif N4by == 0
                push!(qa, Expr(equalop, D_mp, :(+$([:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ]...) ) ) )
            else
                push!(qa, Expr(equalop, D_mp, :(+$(vcat(
                    [:( $(eA(m,n)) * $(eX(n,p))) for n ∈ 1+chunk*N4ax:Nax ],
                    [:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]
                    )...) ) ) )
            end
            if N4by > 0
                push!(qa, :($D_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1:chunk]...) ) )
                for k = 2:N4by
                    push!(qa, :($D_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*(k-1):chunk*k ]...) ) )
                end
                if Nrby > 0 && Nrax == 0
                    push!(qa, :($D_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]...) ) )
                end
            elseif Nrax == 0
                push!(qa, :($D_mp += +$([:( $(eB(m,n)) * $(eY(n,p))) for n ∈ 1+chunk*N4by:Nby ]...) ) ) 
            end
        end
        

        # Approach used for the first block.
        # for p = 1:P, n = 1:N
        #     B_np = eB(n,p)
        #     for m = 1+total_chunked:M
        #         push!(qa, Expr(equalop, eC(m, p), :( $B_np * $(eA(m,n)) ) ) )
        #     end
        # end
    end
end


function mul_kernel(::Type{T}, M, N, P, tA=Val{false}(), tX=Val{false}(), eq = :(=), LA = M*N, LX = N*P, LD = M*P) where T
    q, qa = create_quote()
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LX, :X)
    broadcast_mul_quote!(qa, T, M, N, P, tA, :A, tX, :X, :D, id_symbol, id_symbol, eq )
    insert_linear!(qa, LD, :D)
    # isa(dummy, Bool) || push!(q.args, :C) # worth considering?
    push!(q.args, :D)
    q
end



# """
# Making this generated makes this faster for some reason on 0.6???
# """
# @generated function AXv2!(D, di::BlockIndex{T,mad,pxd,mp}, A, ai::BlockIndex{T,mad,nax,mn}, X, xi::BlockIndex{T,nax,pxd,np}) where {T,mad,pxd,mp,nax,mn,np}
#     # mulandset(D, di::BlockIndex, , )
#     # D[di] = A[ai] * X[xi]
#     quote
#         pa = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{$T,$mad,$nax,$mn}}, pointer_from_objref(A)) + ai.i
#         px = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{$T,$nax,$pxd,$np}}, pointer_from_objref(X)) + xi.i
#         pd = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{$T,$mad,$pxd,$mp}}, pointer_from_objref(D)) #+ di.i

#         Base.unsafe_store!(pd, unsafe_load(pa) * unsafe_load(px) )
#         nothing
#     end
# end
# function AXv3!(D, di::BlockIndex{Float64,cutoff,cutoff,cutoff2}, A, ai::BlockIndex{Float64,cutoff,cutoff,cutoff2}, X, xi::BlockIndex{Float64,cutoff,cutoff,cutoff2})
#     # mulandset(D, di::BlockIndex, , )
#     # D[di] = A[ai] * X[xi]
#     # pa = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(A)) 
#     # px = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(X))
#     # pd = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(D))

#     # Base.unsafe_store!(pd + di.i, unsafe_load(pa + ai.i) * unsafe_load(px + xi.i) )
#     pa = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(A)) + ai.i
#     px = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(X)) + xi.i
#     pd = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(D))  #+ di.i

#     Base.unsafe_store!(pd, unsafe_load(pa) * unsafe_load(px) )
#     nothing
# end
# function AXv4!(ptr, A, ai::BlockIndex{Float64,cutoff,cutoff,cutoff2}, X, xi::BlockIndex{Float64,cutoff,cutoff,cutoff2})
#     # mulandset(D, di::BlockIndex, , )
#     # D[di] = A[ai] * X[xi]
#     # pa = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(A)) 
#     # px = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(X))
#     # pd = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(D))

#     # Base.unsafe_store!(pd + di.i, unsafe_load(pa + ai.i) * unsafe_load(px + xi.i) )
#     pa = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(A)) + ai.i
#     px = Base.unsafe_convert(Ptr{StaticRecursiveMatrix{Float64,cutoff,cutoff,cutoff2}}, pointer_from_objref(X)) + xi.i

#     Base.unsafe_store!(ptr, unsafe_load(pa) * unsafe_load(px) )
#     nothing
# end


"""
These are a list of wrappers.
The idea is that all the calls inside them get inlined, and the compiler should be able to optimize them well.
But they probably shouldn't be inlined to the calling function.
Currently, I have not annotated them with `@noinline`, because there are use cases where maybe it is appropriate to inline them.
Right now, I am trusting the compiler not to.

Passing in a pointer that we store the results into seems to satisfy the aliasing checks.
Calculating the pointer within this function, however fails.

Now, I think proper aliasing checks should be perfectly happy with both versions, because
even were the separate pointers to alias, the instructions quite clearly indicate to load
them separately onto the stack, then perform the calculations with the copies.

Perhaps it matters if the compiler wants the freedom to rearrange loading and storing data?
"""
function AX!(ptrd, A, ai::BlockIndex, X, xi::BlockIndex)
    # mulandset(D, di::BlockIndex, , )
    Base.unsafe_store!(ptrd, A[ai] * X[xi])
    nothing
end
function AX_plus_D!(ptrd, D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex)
    Base.unsafe_store!(ptrd, AX_plus_D(D[di], A[ai], X[xi]))
    nothing
end
function AX_plus_BY!(ptrd, A, ai::BlockIndex, X, xi::BlockIndex,
                                        B, bi::BlockIndex, Y, yi::BlockIndex)
    Base.unsafe_store!(ptrd, AX_plus_BY(A[ai], X[xi], B[bi], Y[yi]))
    nothing
end
function AX_plus_BY_plus_D!(ptrd, D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex,
                                                B, bi::BlockIndex, Y, yi::BlockIndex)
    Base.unsafe_store!(ptrd, AX_plus_BY_plus_D(D[di], A[ai], X[xi], B[bi], Y[yi]))
    nothing
end
function AX_plus_BY_plus_CZ!(ptrd, A, ai::BlockIndex, X, xi::BlockIndex,
                                                B, bi::BlockIndex, Y, yi::BlockIndex,
                                                C, ci::BlockIndex, Z, zi::BlockIndex)
    Base.unsafe_store!(ptrd, AX_plus_BY_plus_CZ(A[ai], X[xi], B[bi], Y[yi], C[ci], Z[zi]))
    nothing
end
function AX_plus_BY_plus_plus_CZ_D!(ptrd, D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex,
                                                        B, bi::BlockIndex, Y, yi::BlockIndex,
                                                        C, ci::BlockIndex, Z, zi::BlockIndex)
    Base.unsafe_store!(ptrd, AX_plus_BY_plus_CZ_D(D[di], A[ai], X[xi], B[bi], Y[yi], C[ci], Z[zi]))
    nothing
end
# function AX!(D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex)
#     # mulandset(D, di::BlockIndex, , )
#     D[di] = A[ai] * X[xi]
#     nothing
# end
# function AX_plus_D!(D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex)
#     D[di] = AX_plus_D(D[di], A[ai], X[xi])
#     nothing
# end
# function AX_plus_BY!(D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex,
#                                         B, bi::BlockIndex, Y, yi::BlockIndex)
#     D[di] = AX_plus_BY(A[ai], X[xi], B[bi], Y[yi])
#     nothing
# end
# function AX_plus_BY_plus_D!(D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex,
#                                                 B, bi::BlockIndex, Y, yi::BlockIndex)
#     D[di] = AX_plus_BY_plus_D(D[di], A[ai], X[xi], B[bi], Y[yi])
#     nothing
# end
# function AX_plus_BY_plus_CZ!(D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex,
#                                                 B, bi::BlockIndex, Y, yi::BlockIndex,
#                                                 C, ci::BlockIndex, Z, zi::BlockIndex)
#     D[di] = AX_plus_BY_plus_CZ(A[ai], X[xi], B[bi], Y[yi], C[ci], Z[zi])
#     nothing
# end
# function AX_plus_BY_plus_plus_CZ_D!(D, di::BlockIndex, A, ai::BlockIndex, X, xi::BlockIndex,
#                                                         B, bi::BlockIndex, Y, yi::BlockIndex,
#                                                         C, ci::BlockIndex, Z, zi::BlockIndex)
#     D[di] = AX_plus_BY_plus_CZ_D(D[di], A[ai], X[xi], B[bi], Y[yi], C[ci], Z[zi])
#     nothing
# end

"""
This kernel does a blocked multiplication.
It gets called when size is <= 3cutoff

A problem I see at the moment is that we are repeatedly loading the same StaticRecursiveMatrices from the heap onto the stack.
So I need to find a clever way of addressing this without having to copy too much in memory between function calls,
    and without letting any functions grow too large.

An alternative formulation would be to load all the sub-blocks once, and then use a series of let blocks?


Also, I do need to type on which D block I am inserting elements into...
"""
function block_kernel(::Type{T}, M, N, P, tA::Val{false}, tX::Val{false}, eq = :(=), LA = M*N, LX = N*P, LD = M*P) where T
    q, qa = create_quote()

    halfcutoff = cutoff ÷ 2

    Mblocks = cld(M, cutoff)
    Nblocks = cld(N, cutoff)
    Pblocks = cld(P, cutoff) #We split these down to 
    
    # 
    # I want this sequence of calls to be for the same method
    # ( with the exception of edge blocks that aren't full size )
    # We do not want to pollute the instruction cache with
    # calls to different functions based on block number.
    # 
    DT = MutableRecursiveMatrix{T,M,P,LD}
    DA = MutableRecursiveMatrix{T,M,N,LA}
    DX = MutableRecursiveMatrix{T,N,P,LX}

    if Nblocks == 1
        if eq == :(=)
            for p = 1:Pblocks, m = 1:Mblocks, h = 1:2
                # dynamic dispatch to (Half)Block and (Half)BlockIndex...
                # Should come up with a better way of doing this.
                push!(qa, :(AX!(point(D, $(BlockIndex(DT, HalfBlock(m,p,h)))),
                                A, $(BlockIndex(DA, Block(m,1))),
                                X, $(BlockIndex(DX, HalfBlock(1,p,h))))) )
            end
        else #eq == :(+=) # assuming; anything else would have to be implemented
            for p = 1:Pblocks, m = 1:Mblocks, h = 1:2
                push!(qa, :(AX_plus_D!(point(D, $(BlockIndex(DT, HalfBlock(m,p,h)))),
                                        D, $(BlockIndex(DT, HalfBlock(m,p,h))),
                                       A, $(BlockIndex(DA, Block(m,1))),
                                       X, $(BlockIndex(DX, HalfBlock(1,p,h))))) )
            end
        end
    elseif Nblocks == 2
        if eq == :(=)
            for p = 1:Pblocks, m = 1:Mblocks, h = 1:2
                push!(qa, :(AX_plus_BY!(point(D, $(BlockIndex(DT, HalfBlock(m,p,h)))),
                                A, $(BlockIndex(DA, Block(m,1))),
                                X, $(BlockIndex(DX, HalfBlock(1,p,h))),
                                A, $(BlockIndex(DA, Block(m,2))),
                                X, $(BlockIndex(DX, HalfBlock(2,p,h))) )) )
            end
        else #eq == :(+=) # assuming; anything else would have to be implemented
            for p = 1:Pblocks, m = 1:Mblocks, h = 1:2
                push!(qa, :(AX_plus_BY_plus_D!(point(D, $(BlockIndex(DT, HalfBlock(m,p,h)))),
                                D, $(BlockIndex(DT, HalfBlock(m,p,h))),
                                A, $(BlockIndex(DA, Block(m,1))),
                                X, $(BlockIndex(DX, HalfBlock(1,p,h))),
                                A, $(BlockIndex(DA, Block(m,2))),
                                X, $(BlockIndex(DX, HalfBlock(2,p,h))) )) )
            end
        end
    else #Nblocks == 3
        if eq == :(=)
            for p = 1:Pblocks, m = 1:Mblocks, h = 1:2
                push!(qa, :(AX_plus_BY_plus_CZ!(point(D, $(BlockIndex(DT, HalfBlock(m,p,h)))),
                                A, $(BlockIndex(DA, Block(m,1))),
                                X, $(BlockIndex(DX, HalfBlock(1,p,h))),
                                A, $(BlockIndex(DA, Block(m,2))),
                                X, $(BlockIndex(DX, HalfBlock(2,p,h))),
                                A, $(BlockIndex(DA, Block(m,3))),
                                X, $(BlockIndex(DX, HalfBlock(3,p,h))) )) )
            end
        else #eq == :(+=) # assuming; anything else would have to be implemented
            for p = 1:Pblocks, m = 1:Mblocks, h = 1:2
                push!(qa, :(AX_plus_BY_plus_CZ_plus_D!(point(D, $(BlockIndex(DT, HalfBlock(m,p,h)))),
                                D, $(BlockIndex(DT, HalfBlock(m,p,h))),
                                A, $(BlockIndex(DA, Block(m,1))),
                                X, $(BlockIndex(DX, HalfBlock(1,p,h))),
                                A, $(BlockIndex(DA, Block(m,2))),
                                X, $(BlockIndex(DX, HalfBlock(2,p,h))),
                                A, $(BlockIndex(DA, Block(m,3))),
                                X, $(BlockIndex(DX, HalfBlock(3,p,h))) )) )
            end
        end
    end
    q
end


"""
The chief problem here is getting the matrices to line up correctly on the N dimension.
In many cases, this should be trivial -- eg, whenever both matrices are approximately square.

When one of them is not "reasonably square", it may need to go through multiple splits before the matrices lign up correctly.

Woah -- hold up. Total N lines up.
So, if one of them isn't reasonably square...how about we just split it until it is, and don't split the other at all?
Doi.
"""
function recursion_mul(::Type{T}, M, N, P, tA::Val{false}=Val(false), tX::Val{false}=Val(false),
                                eq = :(=), LA = M*N, LX = N*P, LD = M*P) where T

    # M1, M2 = split_dim(M)
    # N1, N2 = split_dim(N)
    # P1, P2 = split_dim(P)

    q, qa = create_quote()

    A_square = reasonably_square(M, N)
    X_square = reasonably_square(N, P)

    if A_square && X_square
        # we split each once.
        M1, M2 = split_dim(M)
        N1, N2 = split_dim(N)
        P1, P2 = split_dim(P)
        s = sizeof(T)

        push!(qa, quote
            pd = convert_point(D)
            pa = convert_point(A)
            px = convert_point(X)
            ptr_d11 = PointerRecursiveMatrix{$T,$M1,$P1,$(M1*P1)}(pd)
            ptr_a11 = PointerRecursiveMatrix{$T,$M1,$N1,$(M1*N1)}(pa)
            ptr_x11 = PointerRecursiveMatrix{$T,$N1,$P1,$(N1*P1)}(px)
            mul!( ptr_d11, ptr_a11, ptr_x11 )
            ptr_a12 = PointerRecursiveMatrix{$T,$M1,$N2,$(M1*N2)}(pa + $(s*M *N1))
            ptr_x21 = PointerRecursiveMatrix{$T,$N2,$P1,$(N2*P1)}(px + $(s*N1*P1))
            gemm!( ptr_d11, ptr_a12, ptr_x21 )

            ptr_d21 = PointerRecursiveMatrix{$T,$M2,$P1,$(M2*P1)}(pd + $(s*M1*P1))
            ptr_a21 = PointerRecursiveMatrix{$T,$M2,$N1,$(M2*N1)}(pa + $(s*M1*N1))
            mul!( ptr_d21, ptr_a21, ptr_x11 )
            ptr_a22 = PointerRecursiveMatrix{$T,$M2,$N2,$(M2*N2)}(pa + $(s*(M *N1 + M1*N2)))
            gemm!( ptr_d21, ptr_a22, ptr_x21 )


            ptr_d12 = PointerRecursiveMatrix{$T,$M1,$P2,$(M1*P2)}(pd + $(s*M*P1))
            ptr_x12 = PointerRecursiveMatrix{$T,$N1,$P2,$(N1*P2)}(px + $(s*N*P1))
            mul!( ptr_d12, ptr_a11, ptr_x12 )
            ptr_x22 = PointerRecursiveMatrix{$T,$N2,$P2,$(N2*P2)}(px + $(s*(N*P1+N1*P2)))
            gemm!( ptr_d12, ptr_a12, ptr_x22 )



            ptr_d22 = PointerRecursiveMatrix{$T,$M2,$P2,$(M2*P2)}(pd + $(s*(M*P1 + M1*P2)))
            mul!( ptr_d22, ptr_a21, ptr_x12 )
            gemm!( ptr_d22, ptr_a22, ptr_x22 )
            nothing
        end)

    elseif A_square # X is not square, A is
        if N > P
            # then N is the dimension getting split.
            # So we split X, and proceed.

        else
            # P is the dimension getting split. This means we must split P, until X is reasonably square.
            # We can stop then, without actually splitting A.

        end
    elseif X_square # A is not square, X is
        if N > M
            # N is the dimension getting split.
            # So we split A, and proceed.
        else
            # We now split A until it becomes reasonably square.
            # We do not split X, but instead carry out the multiplication (recursion)
            # as soon A becomes square.

        end
    else # neither A or B are square
        if (M > N) && (N < P)
            # We can split both until they become square, or perhaps just once.
            # Then carry out the multiplication.

        elseif (M < N) && (N > P)
            # We can split until the first of them becomes square
            # then proceed with the recursion.

        elseif (M > N) && (N > P)
            # We split A until it becomes reasonably square.

        else # (M < N) && (N < P)
            # We split X until it becomes reasonably square.

        end
    end

    # A_blocks = block_dim(M, N)
    # B_blocks = block_dim(N, P)
    # while size(A_blocks,2) < size(B_blocks,1) # should be equal, but under certain "pathological" cases they may not be.
    #     A_blocks = split_col(A_blocks)
    # end
    # while size(A_blocks,2) > size(B_blocks,1)
    #     B_blocks = split_row(B_blocks)
    # end
    # # The matrices should now line up.
    # # Is it possible that they don't?
    
    # # Now, we multiply their blocks.

    return q

end


"""
The dummy argument always gets optimized out.
However, when using code_llvm on the function and filling it with a real type (eg, Bool), it causes the SSA names to get incremented by 1 (they start as %0 for the first argument, %1 for the second...). In llvmcall, these values get incremented by one between the first in the body and the argument list, while they do not in the function code returned by llvmcall.
Maybe there is a solution that is less of a hack, but just passing an extra dummy argument to the code_llvm call to get them to line up was a simple solution.
"""
@generated function mul!(D::MutableRecursiveMatrix{T,M,P,LD},
                        A::MutableRecursiveMatrix{T,M,N,LA},
                        X::MutableRecursiveMatrix{T,N,P,LX}) where {T,M,N,P,LA,LX,LD}
                        #dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    maxdim = max(M,N,P)
    # maxdim = max(M,P,N)
    if maxdim <= cutoff # No recursion; we multiply.
        return mul_kernel(T,M,N,P, vistransposed(A), vistransposed(X))
    elseif maxdim <= 3cutoff
        return block_kernel(T,M,N,P, vistransposed(A), vistransposed(X))
    else # Recursion.
        return recursion_mul(T,M,N,P, vistransposed(A), vistransposed(X))
    end
end

"""
By default, Julia refuses to emit SIMD instructions where aliasing is possible.
However, the pointer matrices are only used internally where I can guarantee that they wont alias.

Unlike C/C++, which have `restrict` compiler hints (or Fortran which simply assumes you aren't aliasing),
there's no way for us to make that promise to the compiler. So, instead, the approach is to use `code_llvm`
for the corresponding types where aliasing isn't possible, and then use this code with `llvmcall`.
"""
# @generated function mul!(C::RecursivePointerMatrixOrTranpose{T,M,P,LC},
#                         A::RecursivePointerMatrixOrTranpose{T,M,N,LA},
#                         B::RecursivePointerMatrixOrTranpose{T,N,P,LB},
#                         dummy = Nothing) where {T,M,N,P,LA,LB,LC}
#     iobuffer = IOBuffer()
#     code_llvm(iobuffer, mul!, (dereftype(C), dereftype(A), dereftype(B), Bool))
#     kernel_code = String(iobuffer)
#     codestart = search(mulstring,"{\ntop:\n ")[end]
#     quote
#         Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{$T}, Tuple{Ptr{$T},Ptr{$T},Ptr{$T}}, C.data, A.data, B.data)
#         C
#     end
# end

@generated function gemm!(D::MutableRecursiveMatrix{T,M,P,LD},
                        A::MutableRecursiveMatrix{T,M,N,LA},
                        X::MutableRecursiveMatrix{T,N,P,LX}) where {T,M,N,P,LA,LX,LD}
    # maxdim = max(M,N,P) # two orders, so we can easily force generated function to rerun =P
    maxdim = max(M,P,N)
    if maxdim <= cutoff # No recursion; we multiply.
        return mul_kernel(T,M,N,P, vistransposed(A), vistransposed(X), :(+=))
    elseif maxdim <= 3cutoff
        return block_kernel(T,M,N,P, vistransposed(A), vistransposed(X), :(+=))
    else # Recursion.
        return recursion_mul(T,M,N,P, vistransposed(A), vistransposed(X), :(+=))
    end
end


