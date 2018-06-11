
"""
A: M x N
B: N x P
C: M x P

C = A * B
"""
function broadcast_mul_quote!(qa, ::Type{T}, M, N, P,
                ta::Val{false}=Val(false), A = :A, tb::Val{false}=Val(false), B = :B, tc::Val{false}=Val(false), C = :C,
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
function AX_plus_BY_quote!(qa, ::Type{T}, M, Nax, Nby, P,
                ta::Val{false}=Val(false), A = :A, tx::Val{false}=Val(false), X = :X, B = :B, Y = :Y, tc::Val{false}=Val(false), C = :C,
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
function AX_plus_BY_quotev2!(qa, ::Type{T}, M, Nax, Nby, P,
                ta::Val{false}=Val(false), A = :A, tx::Val{false}=Val(false), X = :X, B = :B, Y = :Y, tc::Val{false}=Val(false), C = :C,
                extract = extract_symbol, insert = extract_symbol, equalop = :(=) ) where T

    # This version seems slower, so the reasoning about loading A and B onto registers and iterating does not seem quite correct.

    eA = (i,j) -> extract(A, sub2ind(false, (M, Nax), i, j))
    eB = (i,j) -> extract(B, sub2ind(false, (M, Nby), i, j))
    eX = (i,j) -> extract(X, sub2ind(false, (Nax, P), i, j))
    eY = (i,j) -> extract(Y, sub2ind(false, (Nby, P), i, j))
    eC = (i,j) -> insert(C, sub2ind(false, (M, P), i, j))
    # Chunk is assumed number we can fit into a register.
    # How can we increase this for avx-512 capable cpus?
    # Do we care about cutting for avx-incapable?
    chunk = max(1, 32 ÷ sizeof(T))
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
    #
    
    Ntotal = Nax + Nby
    
    registers_for_AB = 16 - Num_chunks - 2
    AB_register_requirement = Ntotal * Num_chunks
    n_per_iter = registers_for_AB ÷ Num_chunks
    shared_n_per_iter = n_per_iter >> 1
    n_per_iter = shared_n_per_iter << 1 # make it even, rounding down.

    iterations = cld( Ntotal, n_per_iter )
    # not M * Ntotal ÷ chunk, because we separate out the M ÷ chunk remainder

    # We will break up our computation
    # into this many iterations, so we can hopefully hold these chunks of A and B in memory.
    # iterations = cld(AB_register_requirement, registers_for_AB)

    smaller_N = min(Nax, Nby)
    min_solo_n, max_solo_n = smaller_N, 0
    last_n = 0
    # remainder = 0
    for iter = 1:iterations
        new_n = shared_n_per_iter * iter
        if new_n > smaller_N
            shared_n = smaller_N
            min_solo_n, max_solo_n = 1+max(smaller_N,last_n), min(2new_n, Ntotal)
        else
            shared_n = new_n
        end

        for p = 1:P
            X_np = eX(1,p)
            Y_np = eY(1,p)
            if iter == 1
                for c = 1:Num_chunks, m = 1+(c-1)*chunk:c*chunk
                    C_mp = eC(m, p)
                    push!(qa, Expr( equalop, C_mp, :( $X_np * $(eA(m,1)) ) ) ) # equalop -- we might be writing C the first time.
                    push!(qa, Expr( :(+=),   C_mp, :( $Y_np * $(eB(m,1)) ) ) )
                end
            end
            for n = max(2,last_n+1):shared_n
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
            if Nax != Nby && smaller_N == shared_n
  
                if Nax > Nby
                    for n = min_solo_n:max_solo_n
                        X_np = eX(n,p)
                        for m = 1:total_chunked
                            push!(qa, Expr( :(+=), eC(m, p), :( $X_np * $(eA(m,n)) ) ) )
                        end
                    end
                else
                    for n = min_solo_n:max_solo_n
                        Y_np = eY(n,p)
                        for m = 1:total_chunked
                            push!(qa, Expr( :(+=), eC(m, p), :( $Y_np * $(eB(m,n)) ) ) )
                        end
                    end
                end
            end
        end
        last_n = shared_n
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
    AX_plus_BY_quote!(qa, T, M, Nax, Nby, P,
    Val(false), :A, Val(false), :X, :B, :Y, Val(false), :C, id_symbol, id_symbol, :(=) )
    # insert_linear!(qa, LC, :C)
    # push!(q.args, :C)
    LC = M*P
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LC}( Base.Cartesian.@ntuple $LC C ) ))
    q

end
@generated function AX_plus_BY2(A::StaticRecursiveMatrix{T,M,Nax,LA}, X::StaticRecursiveMatrix{T,Nax,P,LX},
                                B::StaticRecursiveMatrix{T,M,Nby,LB}, Y::StaticRecursiveMatrix{T,Nby,P,LY}) where {T,M,Nax,Nby,P,LA,LX,LB,LY}
    q, qa = create_quote(true)
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A) # second arg is the name of the variable we're extracting.
    extract_linear!(qa, LB, :A, :B, 0, LA)
    for p = 1:P
        extract_linear!(qa, Nax, :X, :X, (p-1)*Nax, (p-1)*(Nax+Nby)      )
        extract_linear!(qa, Nby, :X, :Y, (p-1)*Nby, (p-1)*(Nax+Nby) + Nax)
    end
    broadcast_mul_quote!(qa, T, M, Nax + Nby, P, Val(false), :A, Val(false), :X, Val(false), :C, id_symbol, id_symbol, :(=) )
    LC = M*P
    # insert_linear!(qa, LC, :C)
    
    # push!(q.args, :C)
    # q
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LC}( Base.Cartesian.@ntuple $LC C ) ))
    q
end

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
    broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :X, Val(false), :D, id_symbol, id_symbol, :(=) )
    LD = M*P
    # insert_linear!(qa, LD, :D)
    
    # push!(q.args, :D)
    # q
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LD}( Base.Cartesian.@ntuple $LD D ) ))
    q
end

@generated function Base.:*(A::StaticRecursiveMatrix{T,M,N,LA}, B::StaticRecursiveMatrix{T,N,P,LB}) where {T,M,N,P,LA,LB}
    q, qa = create_quote()
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LB, :B)
    broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :B, Val(false), :C, id_symbol, id_symbol, :(=) )
    # insert_linear!(qa, LC, :C)
    
    # push!(q.args, :C)
    # q
    LC = M*P
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LC}( Base.Cartesian.@ntuple $LC C ) ))
    q
end
@generated function muladd(A::StaticRecursiveMatrix{T,M,N,LA},
                            B::StaticRecursiveMatrix{T,N,P,LB},
                            C::StaticRecursiveMatrix{T,M,P,LC}) where {T,M,N,P,LA,LB,LC}

    q, qa = create_quote()
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LB, :B)
    # extract_linear!(qa, LC, :C)
    broadcast_mul_quote!(qa, T, M, N, P, Val(false), :A, Val(false), :B, Val(false), :C, id_symbol, extract_symbol, :(+=) )
    # insert_linear!(qa, LC, :C)
    
    # push!(q.args, :C)
    # q
    push!(q.args, :( StaticRecursiveMatrix{$T,$M,$P,$LC}( Base.Cartesian.@ntuple $LC C ) ))
    q
end


function AX_plus_BY_plus_DZ_quote!(qa, ::Type{T}, M, Nax, Nby, Ndz, P,
                ta::Val{false}=Val(false), A = :A, tx::Val{false}=Val(false), X = :X, B = :B, Y = :Y, D = :D, Z = :Z,
                tc::Val{false}=Val(false), C = :C,
                extract = extract_symbol, insert = extract_symbol, equalop = :(=) ) where T

    eA = (i,j) -> extract(A, sub2ind(false, (M, Nax), i, j))
    eB = (i,j) -> extract(B, sub2ind(false, (M, Nby), i, j))
    eD = (i,j) -> extract(D, sub2ind(false, (M, Ndz), i, j))
    eX = (i,j) -> extract(X, sub2ind(false, (Nax, P), i, j))
    eY = (i,j) -> extract(Y, sub2ind(false, (Nby, P), i, j))
    eZ = (i,j) -> extract(Z, sub2ind(false, (Ndz, P), i, j))
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
        XYZ_np = eX(1,p)
        for m = 1:total_chunked
            push!(qa, Expr( equalop, eC(m, p), :( $XYZ_np * $(eA(m,1)) ) ) )
        end
        for n = 2:Nax
            XYZ_np = eX(n,p)
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eC(m, p), :( $XYZ_np * $(eA(m,n)) ) ) )
            end
        end
        for n = 1:Nby
            XYZ_np = eY(n,p)
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eC(m, p), :( $XYZ_np * $(eB(m,n)) ) ) )
            end
        end
        for n = 1:Ndz
            XYZ_np = eZ(n,p)
            for m = 1:total_chunked
                push!(qa, Expr( :(+=), eC(m, p), :( $XYZ_np * $(eD(m,n)) ) ) )
            end
        end
    end

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


function mul_kernel(::Type{T}, M, N, P, tA=Val{false}(), tB=Val{false}(), tC=Val{false}(), eq = :(=), LA = M*N, LB = N*P, LC = M*P) where T
    q, qa = create_quote()
    @static if VERSION > v"0.7-"
        pushfirst!(q.args, :(Base.@_inline_meta))
    else
        insert!(q.args, 1, :(Base.@_inline_meta))
    end
    extract_linear!(qa, LA, :A)
    extract_linear!(qa, LB, :B)
    broadcast_mul_quote!(qa, T, M, N, P, tA, :A, tB, :B, tC, :C, id_symbol, id_symbol, eq )
    insert_linear!(qa, LC, :C)
    # isa(dummy, Bool) || push!(q.args, :C) # worth considering?
    push!(q.args, :C)
    q
end

function block_kernel(M, N, P, tA=false, tB=false, tC=false, eq = :(=), LA = M*N, LB = N*P, LC = M*P)
    q, qa = create_quote()
    # Because this gets called when max(M,N,P) <= 2cutoff, we go ahead and split all dims > cutoff
    if min(M,N,P) > cutoff
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
                        B::RecursiveMatrixOrTranpose{T,N,P,LB}) where {T,M,N,P,LA,LB,LC}
                        #dummy = Nothing) where {T,M,N,P,LA,LB,LC}
    maxdim = max(M,N,P)
    if maxdim <= cutoff # No recursion; we multiply.
        return mul_kernel(T,M,N,P, vistransposed(A), vistransposed(B), vistransposed(C))
    elseif maxdim <= 2cutoff
        return block_kernel(M,N,P, vistransposed(A), vistransposed(B), vistransposed(C))
    else # Recursion.
        return recursion_mul(M,N,P, vistransposed(A), vistransposed(B), vistransposed(C))
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
        return mul_kernel(M, N, P, vistransposed(A), vistransposed(B), vistransposed(C), :(+=) )
    else # Recursion.
        return block_kernel(M,N,P, vistransposed(A), vistransposed(B), vistransposed(C), :(+=) )
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


            
            Base.llvmcall($(kernel_code[codestart:end-3]), Ptr{Ptr{Int8}}, Tuple{Ptr{Ptr{Int8}},Ptr{Ptr{Int8}},Ptr{Ptr{Int8}}}, pointer(Base.unsafe_convert(Ptr{Int8}, C.data)), Ref(Base.unsafe_convert(Ptr{Int8}, A.data)), Ref(Base.unsafe_convert(Ptr{Int8}, B.data)))
            
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