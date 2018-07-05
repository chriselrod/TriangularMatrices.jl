extract_symbol(A, i) = :($A[$i])
id_symbol(A, i) = Symbol(A, :_, i)


function create_quote(fast::Bool = true)
    if fast
        q = quote @fastmath @inbounds begin end end
        @static if VERSION > v"0.7-"
            qa = q.args[2].args[3].args[3].args
        else
            qa = q.args[2].args[2].args[2].args
        end
    else
        q = quote @inbounds begin end end
        @static if VERSION > v"0.7-"
            qa = q.args[2].args[3].args
        else
            qa = q.args[2].args[2].args
        end
    end
    q, qa
end

function create_quote(A_i, A)
    q, qa = create_quote()
    push!(qa, :($A_i = $A[1]))
    q, qa
end

function initial_quote(N2::Int, prefix::Symbol = :A)
    prefix_ = Symbol(prefix, :_)
    q, qa = create_quote(Symbol(prefix_, 1), prefix)
    for i ∈ 2:N2
        push!(qa, :($(Symbol(prefix_, i)) = $(prefix)[$i]))
    end
    q, qa
end

function extract_linear!(qa, N, prefix = :B, extractedvar = prefix, ind_offset = 0, label_offset = 0)
    prefix_ = Symbol(prefix, :_)
    for i ∈ 1:N
        push!(qa, :($(Symbol(prefix_, i + label_offset)) = $(extractedvar)[$(i+ind_offset)]))
    end
    qa
end
function insert_linear!(qa, N, prefix = :B, ind_offset::Int = 0, label_offset::Int = 0)
    prefix_ = Symbol(prefix, :_)
    for i ∈ 1:N
        push!(qa, :( $(prefix)[$(i+ind_offset)] = $(Symbol(prefix_, i + label_offset))) )
    end
    qa
end
function insert_linear!(qa, N, output::Symbol, prefix::Symbol = :B)
    prefix_ = Symbol(prefix, :_)
    for i ∈ 1:N
        push!(qa, :( $(output)[$(i)] = $(Symbol(prefix_, i))) )
    end
    qa
end

function extract_transpose!(qa, N, N2, prefix::Symbol = :A)
    prefix_ = Symbol(prefix, :_)
    A_i = Symbol(prefix_, 1)
    ind_extract = 0
    for i ∈ 1:N2, j ∈ 0:N-1
        ind_extract += 1
        A_i = Symbol(prefix_, i + j*N2)
        push!(qa, :($A_i = $(prefix)[$ind_extract]))
    end
    qa
end

function extract_transpose_direct_indexing!(qa, N, N2, prefix::Symbol = :A)
    prefix_ = Symbol(prefix, :_)
    A_i = Symbol(prefix_, 1)
    ind_extract = 0
    for i ∈ 1:N2, j ∈ 0:N-1
        ind_extract += 1
        A_i = Symbol(prefix_, i + j*N2)
        push!(qa, :($A_i = $(prefix)[$ind_extract]))
    end
    qa
end

function output_transpose(N, N2, T, prefix = :C)
    prefix_ = Symbol(prefix, :_)
    out = :( StaticArrays.SArray{Tuple{$N,$N2},$T,2,$(N*N2)}( ( $(Symbol(prefix_, 1)) , ) ) )
    outa = out.args[2].args
    for j ∈ 1:N-1
        push!(outa, Symbol(prefix_, 1 + j*N2)  )
    end
    for i ∈ 2:N2, j ∈ 0:N-1
        push!(outa, Symbol(prefix_, i + j*N2)  )
    end
    out
end

# @inline pointer_from_data(r::Base.RefValue) = pointer_from_objref(r)
# @inline pointer_from_data(r::Ptr) = r

# @inline point(r::Base.RefValue{NTuple{N,T}}) where {N,T} = Base.unsafe_convert(Ptr{T}, pointer_from_objref(r))
# @inline point(p::Ptr) = p


@generated ValProd(::Val{M}, ::Val{N}) where {M,N} = Val{M*N}()
@generated function ValDiv(::Val{L}, ::Val{M}) where {L,M}
    N, r = divrem(L, M)
    @assert r == 0
    Val{N}()
end
