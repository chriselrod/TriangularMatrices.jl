

function initial_quote(N2::Int, prefix::Symbol = :A_)
    A_i = Symbol(prefix, 1)
    q = quote
        @fastmath @inbounds begin
            $A_i = A.data[1]
        end
    end
    @static if VERSION > v"0.6.9"
        qa = q.args[2].args[3].args[3].args
    else
        qa = q.args[2].args[2].args[2].args
    end
    for i ∈ 2:N2
        A_i = Symbol(prefix, i)
        push!(qa, :($A_i = A.data[$i]))
    end
    q, qa
end