
"""
This function defines how a dimension is split.
It is used everywhere dimensions are split, therefore it (and the cutoff) define the recursion pattern.

The current behavior is to try and get dimensions into primarily blocks of 8.
"""
function split_dim(N)
    N <= cutoff && return N, 0
    number_cutoff_blocks, remainder = divrem(N,cutoff)

    half_of_cutoff_blocks, odd = divrem(number_cutoff_blocks, 2)

    8(half_of_cutoff_blocks+odd), 8half_of_cutoff_blocks + remainder
end

function reasonably_square(mindim, maxdim)
    # mindim, maxdim = minmax(M, N)
    round(Int, √2 * mindim, RoundUp) >= maxdim
end

# function blocks_per_dim(N::Integer)
#     blocks = 1
#     block_size = N
#     cuts = 0
#     # bit twiddling is 50% faster
#     while cutoff < block_size
#         cuts += 1
#         blocks <<= 1 # *=2
#         block_size >>= 1 #  ÷=2
#     end
#     blocks + min(max(0, N - cutoff << cuts), blocks) #N - cutoff * 2^cuts
# end
# @generated blocks_per_dim(::Val{N}) where N = blocks_per_dim(N) # Is this needed?


# ### Implemented so that the leading blocks are larger. This is because we tend to do just slightly less with them.
# function block_position(N::Integer, i::Integer)
#     num_blocks = blocks_per_dim(N)
#     1 + num_blocks - cld(num_blocks*(N+1-i),N)
# end
# function block_position(N::Integer, i::Integer, j::Integer) #Not faster on 0.6?
#     num_blocks = blocks_per_dim(N)
#     1+num_blocks - cld(num_blocks*(N+1-i),N), 1+num_blocks - cld(num_blocks*(N+1-j),N)
# end

# function block_position(N, i)
#     position = 1
#     while N > cutoff
#         N, r = divrem(N, 2)
#         Npr = N + r
# #        @show i, N, Npr
#         if i > Npr
#             position += blocks_per_dim(Npr)#round(Int, log2(Npr/cutoff), RoundUp)
#             i -= Npr
#         else
#             N = Npr
#         end
#     end
#     position
# end

function dense_sub2ind(dims::NTuple{2}, i, j, offset = 0)
    M, N = dims
    minind, maxind = minmax(M, N)
    maxind <= cutoff && return (j-1)*M + i + offset
    if reasonably_square(minind, maxind) && minind > cutoff #divide both
        M1, M2 = split_dim(M)
        N1, N2 = split_dim(N)
        if i <= M1 && j <= N1 # upper left block
            return dense_sub2ind( (M1, N1), i, j, offset )
        elseif j <= N1 # lower left block
            return dense_sub2ind( (M2, N1), i - M1, j, offset + M1 * N1 )
        elseif i <= M1 # upper right block
            return dense_sub2ind( (M1, N2), i, j - N1, offset + M * N1 )
        else # lower right block
            return dense_sub2ind( (M2, N2), i - M1, j - N1, offset + M * N1 + M1 * N2 )
        end
    elseif M > N #divide M
        M1, M2 = split_dim(M)
        if i <= M1
            return dense_sub2ind( (M1, N), i, j, offset)
        else
            return dense_sub2ind( (M2, N), i - M1, j, offset + M1 * N)
        end
    else #divide N
        N1, N2 = split_dim(N)
        if j <= N1
            return dense_sub2ind( (M, N1), i, j, offset)
        else
            return dense_sub2ind( (M, N2), i, j - N1, offset + M * N1)
        end
    end
end

"""
Generates a recursive quote.
The function keeps track of
    M, number of rows
    N, number of columns
    linear_offset, linear offset into the data of original matrix
    row_offset, Cartesian offset of rows
    column_offset, Cartesian offset of columns
Through the recursion. When the recurion terminates, it walks the expression and replaces all instances of those symbols with their current value.
"""
function dense_recursion_quote(q, M, N, linear_offset = 0, row_offset = 0, column_offset = 0)
    minind, maxind = minmax(M, N)

    # The recursion terminates, and we substitute expressions for their value.
    if maxind <= cutoff
        return MacroTools.postwalk(x -> begin
            if isa(x, Symbol)
                if x == :M
                    return M
                elseif x == :N
                    return N
                elseif x == :linear_offset
                    return linear_offset
                elseif x == :row_offset
                    return row_offset
                elseif x == :column_offset
                    return column_offset
                else
                    return x
                end
            else
                return x
            end
        end, q)
    end
    
    if reasonably_square(minind, maxind) && minind > cutoff # we divide the matrix into four blocks.
        M1, M2 = split_dim(M)
        N1, N2 = split_dim(N)
        return :(if i <= $(M1 + row_offset) && j <= $(N1 + column_offset) # enter block 1
                $(dense_recursion_quote(q, M1, N1, linear_offset, row_offset, column_offset))
            elseif j <= $(N1 + column_offset) # enter block 2
                $(dense_recursion_quote(q, M2, N1, linear_offset + M1 * N1, row_offset + M1, column_offset))
            elseif i <= $(M1 + row_offset) # enter block 3
                $(dense_recursion_quote(q, M1, N2, linear_offset + M * N1, row_offset, column_offset + N1))
            else # enter block 4
                $(dense_recursion_quote(q, M2, N2, linear_offset + M * N1 + M1 * N2, row_offset + M1, column_offset + N1))
            end)
    elseif M > N # we are splitting the matrix into two blocks stacked on top of one another.
        M1, M2 = split_dim(M)
        return :(if i <= $(M1 + row_offset)
                $(dense_recursion_quote(q, M1, N, linear_offset, row_offset, column_offset))
            else
                $(dense_recursion_quote(q, M2, N, linear_offset + M1 * N, row_offset + M1, column_offset))
            end)
    else # we are splitting the matrix into two blocks, side by side.
        N1, N2 = split_dim(N)
        return :(if j <= $(N1 + column_offset)
                $(dense_recursion_quote(q, M, N1, linear_offset, row_offset, column_offset))
            else
                $(dense_recursion_quote(q, M, N2, linear_offset + M * N1, row_offset, column_offset + N1))
            end)
    end
end


function dense_sub2ind_quote(M, N) 
    dense_recursion_quote(
        :( (j-(column_offset+1) ) * M + i + (linear_offset-row_offset) ),
        M, N)
end

"""
Cartesian indexing isn't recomended, but it is convenient for printing.
The approach for sub2ind here incurs branches, which would disable SIMD.
Therefore, Cartesian indexing is strongly discouraged for hot loops.
Still, it is reasonably fast -- only a handful of ns.
"""
# @generated dense_sub2ind(::Val{M}, ::Val{N}, i, j) where {N,M} = dense_sub2ind_quote(M, N)
@generated function dense_sub2ind(::Val{M}, ::Val{N}, i, j) where {N,M}
    dense_sub2ind_quote(M, N)
end



function triangle_sub2ind(N, i, j, offset = 0)
    N <= cutoff && return small_triangle(j) + i + offset
    N1, N2 = split_dim(N)
    if max(i,j) <= N1 # enter block 1 # triangle / symmetric block
        return triangle_sub2ind(N1, i, j, offset)
    elseif i <= N1 # enter block 3 # square block
        return dense_sub2ind((N1, N2), i, j - N1, offset + big_triangle(N1))
    elseif j <= N1 # enter block 2 #symmetric or 0 block # Not supported; check i, j before calling.
        throw("Indexing into this half not supported.")
    else # enter block 4 # triangle / symmetric block
        return triangle_sub2ind(N2, i - N1, j - N1, offset + big_triangle(N1) + N1*N2)
    end
end




"""
The triangular (or Symmetric) matrix is assumed to be square, with only upper triangular elements represented.
Lower triangular matrices are simply a view of an upper triangular matrix that pretends the storage is row-major.

This function takes in two quotes, to allow for different behavior in a triangular / diagonal block, and the offdiagonal blocks.

When the recursion terminates, the following symbols are replaced with their current value:
    M -- only used for the off diagonal blocks. See dense_recursion_quote for details.
    N -- the current size of the matrix block.
    linear_offset -- linear offset within the original data matrix.
    row_offset -- how many rows the block is offset
    column_offset -- how many columns the block has been offset
"""
function triangle_recursion_quote(triangle_q, rectangle_q, N, linear_offset = 0, row_offset = 0, column_offset = 0)
    if N <= cutoff
        return MacroTools.postwalk(x -> begin
            if isa(x, Symbol)
                if x == :N
                    return N
                elseif x == :linear_offset
                    return linear_offset
                elseif x == :row_offset
                    return row_offset
                elseif x == :column_offset
                    return column_offset
                else
                    return x
                end
            else
                return x
            end
        end, triangle_q)
    end

    N1, N2 = split_dim(N)
    :(if max( i - $row_offset, j - $column_offset) <= $N1 # enter block 1 # triangle / symmetric block
            $(triangle_recursion_quote(triangle_q, rectangle_q, N1, linear_offset, row_offset, column_offset))
        elseif i <= $(N1 + row_offset) # enter block 3 # square block
            $(dense_recursion_quote(rectangle_q, N1, N2, linear_offset + big_triangle(N1), row_offset, column_offset + N1))
        else # enter block 4 # triangle / symmetric block
            $(triangle_recursion_quote(triangle_q, rectangle_q, N2,
                linear_offset + big_triangle(N1) + N1*N2, row_offset + N1, column_offset + N1))
        end)

end

function triangle_sub2ind_quote(N)
    triangle_recursion_quote(
        :( small_triangle( j-column_offset ) + i + (linear_offset-row_offset) ),
        :( (j-(column_offset+1) ) * M + i + (linear_offset-row_offset) ),
        N
    )
end


@generated function triangle_sub2ind(::Val{N}, i, j) where N
    quote
        Base.@_inline_meta
        @boundscheck i > j && throw("Indexing into this half not supported.")
        $(triangle_sub2ind_quote(N))
    end
end



function block_dim(M,N)
    mindim, maxdim = minmax(M,N)
    # The "reasonably square" may need improvement
    # For example, a 91 x 64 breaks down into
    # 46 x 64 and 45 x 64; the 45 x 64 also fails
    # the criterion, and breaks into two 45 x 32 blocks
    # How bad is this wonky pattern?
    # Should I round in the comparison?
    if reasonably_square(mindim, maxdim) && mindim > cutoff
        M1, M2 = split_dim(M)
        N1, N2 = split_dim(N)
        out = [ (M1,N1) (M1, N2); (M2,N1) (M2,N2) ] 
    elseif maxdim <= cutoff
        out = Matrix{Tuple{Int,Int}}(undef,1,1)
        out[1] = (M,N)
    elseif M > N
        out = Matrix{Tuple{Int,Int}}(undef,2,1)
        M1, M2 = split_dim(M)
        out[1] = (M1, N)
        out[2] = (M2, N)
    else
        out = Matrix{Tuple{Int,Int}}(undef,1,2)
        N1, N2 = split_dim(N)
        out[1] = (M, N1)
        out[2] = (M, N2)
    end
    out
end

# function split_dim!(n_dim::Vector{Int}, i)
#     N1, N2 = split_dim(n_dim[i])
#     n_dim[i] = N2
#     insert!(n_dim, i, N1)
# end
# function split_row(dims::Matrix{Tuple{Int}}, i)
#     N1, N2 = split_dim(n_dim[i])
#     n_dim[i] = N2
#     insert!(n_dim, i, N1)
# end
function split_col(dims::Matrix{Tuple{Int,Int}})
    num_rowblocks, num_colblocks = size(dims)
    @assert num_colblocks == 1
    vcat([block_dim(dims[i]...) for i = 1:length(dims)]...)
end
function split_row(dims::Matrix{Tuple{Int,Int}})
    num_rowblocks, num_colblocks = size(dims)
    @assert num_rowblocks == 1
    hcat([block_dim(dims[i]...) for i = 1:length(dims)]...)
end



# function block_structure(N, offset = 0)
#     bounds = Int[]
#     block_structure!(bounds, N)
#     # sort!(bounds)
# end
# function block_structure!(bounds, N, offset = 0)
#     if N > cutoff
#         N, r = divrem(N, 2)
#         Npr = N + r
#         push!(bounds, Npr + offset)
#         block_structure!(bounds, Npr, offset)
#         block_structure!(bounds, N, offset + Npr)
#         return bounds
#     else
#         return bounds
#     end
# end

# function identify_block_placement_quote(N)
#     bounds = block_structure(N)
#     num_bounds = length(bounds)
#     num_blocks = num_bounds + 1
#     q = quote
#         b_0 = 0
#         $(Symbol(:b_, num_blocks)) = $N
#     end
#     for k ∈ 1:num_bounds
#         push!(q.args, :( @inbounds $(Symbol(:b_, k)) = $(bounds[k])) )
#     end
#     push!(q.args,
#         :( Base.Cartesian.@nif $num_blocks d -> (i <= b_d)  d -> ( row_block_id = d; row_lbound = b_{d-1}; row_ubound = b_d ) d -> ( row_block_id = d; row_lbound = b_{d-1}; row_ubound = b_d ) )
#     )
#     push!(q.args,
#         :( Base.Cartesian.@nif $num_blocks d -> (j <= b_d)  d -> ( col_block_id = d; col_lbound = b_{d-1}; col_ubound = b_d ) d -> ( col_block_id = d; col_lbound = b_{d-1}; col_ubound = b_d ) )
#     )
#     push!(q.args, :(block_i = i - row_lbound))
#     push!(q.args, :(block_j = j - col_lbound))
#     q
# end
# function identify_block_placement_quote(M, N)
#     rowbounds = block_structure(M)
#     colbounds = block_structure(N)
#     num_rowbounds = length(rowbounds)
#     num_colbounds = length(colbounds)
#     num_rowblocks = num_rowbounds + 1
#     num_colblocks = num_colbounds + 1
#     q = quote
#         r_0 = 0
#         c_0 = 0
#         $(Symbol(:r_, num_rowblocks)) = $M
#         $(Symbol(:c_, num_colblocks)) = $N
#     end
#     for k ∈ 1:num_rowbounds
#         push!(q.args, :( @inbounds $(Symbol(:r_, k)) = $(rowbounds[k])) )
#     end
#     for k ∈ 1:num_colbounds
#         push!(q.args, :( @inbounds $(Symbol(:c_, k)) = $(colbounds[k])) )
#     end
#     push!(q.args,
#         :( Base.Cartesian.@nif $num_rowblocks d -> (i <= r_d)  d -> ( row_lbound = r_{d-1}; row_ubound = r_d ) d -> (row_lbound = r_{d-1}; row_ubound = r_d) )
#     )
#     push!(q.args,
#         :( Base.Cartesian.@nif $num_colblocks d -> (j <= c_d)  d -> ( col_lbound = c_{d-1}; col_ubound = c_d ) d -> (col_lbound = c_{d-1}; col_ubound = c_d) )
#     )
#     push!(q.args, :(block_i = i - row_lbound))
#     push!(q.args, :(block_j = j - col_lbound))
#     q
# end



# @generated function triangle_sub2ind(::Val{N}, i, j) where N
#     N <= cutoff && return :(small_triangle(j)+i)
#     q = identify_block_placement_quote(N)
#     push!(q.args, :(ind = big_triangle(col_lbound) + (col_ubound - col_lbound) * row_lbound))
#     push!(q.args, :(if row_block_id == col_block_id # then we're in a symmetric portion
#         ind += small_triangle(block_j) + block_i
#     else #then we are in square block
#         ind += (row_ubound-row_lbound) * (block_j-1) + block_i
#     end))
#     push!(q.args, :ind)
#     q
# end

# @generated dense_sub2ind(::Val{M}, ::Val{N}, i, j) where {M,N} 
#     max(M,N) <= cutoff && return :(j*$(M-1)+i)
#     q = identify_block_placement_quote(M,N)
#     push!(q.args :( col_lbound*$N + row_lbound*(col_ubound-col_lbound) + (row_ubound - row_lbound)*(block_j-1) + block_i ))
#     q
# end

# @generated function triangle_sub2ind(::Val{N}, i, j) where N
#     N <= cutoff && return :(small_triangle(j)+i)
#     num_blocks = blocks_per_dim(N)
#     quote
#         row_block_id = $(1+num_blocks) - cld($num_blocks*($(N+1)-i),$N)
#         col_block_id = $(1+num_blocks) - cld($num_blocks*($(N+1)-j),$N)

#         row_block_lbound = cld((row_block_id-1)*$N,$num_blocks)
#         # row_block_ubound = cld( row_block_id   *$N,$num_blocks)
#         col_block_lbound = cld((col_block_id-1)*$N,$num_blocks)
#         col_block_ubound = cld( col_block_id   *$N,$num_blocks)

#         ind = big_triangle(col_block_lbound) + (col_block_ubound - col_block_lbound) * row_block_lbound

#         block_i = i - row_block_lbound
#         block_j = j - col_block_lbound

#         if row_block_id == col_block_id # then we're in a symmetric portion
#             ind += small_triangle(block_j) + block_i
#         else #then we are in square block
#             nrows_in_block = cld( row_block_id   *$N,$num_blocks)  - row_block_lbound
#             ind += nrows_in_block * (block_j-1) + block_i
#         end
#         ind
#     end
# end