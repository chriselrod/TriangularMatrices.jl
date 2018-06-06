
# struct SquareSummary
#     N::Int
#     blocks::Int
# end

# struct RectangleSummary
#     n_rows::Int
#     n_cols::Int
#     row_blocks::Int
#     col_blocks::Int
# end

# SquareSummary(N) = SquareSummary(N, blocks_per_dim(N))

# function RectangleSummary(n_rows, n_cols)
#     RectangleSummary(
#         n_rows, n_cols,
#         blocks_per_dim(n_rows),
#         blocks_per_dim(n_cols)
#     )
# end

function blocks_per_dim(N::Integer)
    blocks = 1
    block_size = N
    cuts = 0
    # bit twiddling is 50% faster
    while cutoff < block_size
        cuts += 1
        blocks <<= 1 # *=2
        block_size >>= 1 #  ÷=2
    end
    blocks + min(max(0, N - cutoff << cuts), blocks) #N - cutoff * 2^cuts
end
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

function block_position(N, i)
    position = 1
    while N > cutoff
        N, r = divrem(N, 2)
        Npr = N + r
#        @show i, N, Npr
        if i > Npr
            position += blocks_per_dim(Npr)#round(Int, log2(Npr/cutoff), RoundUp)
            i -= Npr
        else
            N = Npr
        end
    end
    position
end

function dense_sub2ind(dims::NTuple{2}, i, j, offset = 0, joffset = 1)
    M, N = dims
    minind, maxind = minmax(M, N)
    maxind*minind < abs2(cutoff) && return (j-joffset)*M + i + offset
    if √2*minind > maxind #divide both
        Mh, Mr = divrem(M, 2)
        Nh, Nr = divrem(N, 2)
        Mh_p_Mr = Mh + Mr
        Nh_p_Nr = Nh + Nr
        if i <= Mh_p_Mr && j <= Nh_p_Nr
            return dense_sub2ind( (Mh_p_Mr, Nh_p_Nr), i, j, offset, joffset )
        elseif j <= Nh_p_Nr
            return dense_sub2ind( (Mh, Nh_p_Nr), i, j, offset + Mh_p_Mr * (Nh_p_Nr-1), joffset )
        elseif i <= Mh_p_Mr
            return dense_sub2ind( (Mh_p_Mr, Nh), i, j, offset + M * Nh_p_Nr, joffset + Nh_p_Nr )
        else
            return dense_sub2ind( (Mh, Nh), i, j, offset + M * Nh_p_Nr + Mh_p_Mr * (Nh-1), joffset + Nh_p_Nr )
        end
    elseif M > N #divide M
        Mh, Mr = divrem(M, 2)
        Mh_p_Mr = Mh + Mr
        if i < Mh_p_Mr
            return dense_sub2ind( (Mh_p_Mr, N), i, j, offset, joffset)
        else
            return dense_sub2ind( (Mh, N), i, j, offset + Mh_p_Mr * (N-1), joffset)
        end
    else #divide N
        Nh, Nr = divrem(N, 2)
        Nh_p_Nr = Nh + Nr
        if j < Nh_p_Nr
            return dense_sub2ind( (M, Nh_p_Nr), i, j, offset, joffset)
        else
            return dense_sub2ind( (M, Nh), i, j, offset + M * Nh_p_Nr, joffset + Nh_p_Nr)
        end
    end
end
function dense_sub2ind_quote(M, N, offset = 0, joffset = 1)
    minind, maxind = minmax(M, N)
    maxind*minind < abs2(cutoff) && return :( (j-$joffset ) * $M + i + $offset )
    
    if √2*minind > maxind # we divide the matrix into four blocks.
        Mh, Mr = divrem(M, 2)
        Nh, Nr = divrem(N, 2)
        Mh_p_Mr = Mh + Mr
        Nh_p_Nr = Nh + Nr
        return :(if i <= $Mh_p_Mr && j <= $Nh_p_Nr # enter block 1
                $(dense_sub2ind_quote(Mh_p_Mr, Nh_p_Nr, offset, joffset))
            elseif j <= $Nh_p_Nr # enter block 2
                $(dense_sub2ind_quote(Mh, Nh_p_Nr, offset + Mh_p_Mr * (Nh_p_Nr-1), joffset))
            elseif i <= $Mh_p_Mr # enter block 3
                $(dense_sub2ind_quote(Mh_p_Mr, Nh, offset + M * Nh_p_Nr, joffset + Nh_p_Nr))
            else # enter block 4
                $(dense_sub2ind_quote(Mh, Nh, offset + M * Nh_p_Nr + Mh_p_Mr * (Nh-1), joffset + Nh_p_Nr))
            end)
    elseif M > N # we are splitting the matrix into two blocks stacked on top of one another.
        Mh, Mr = divrem(M, 2)
        Mh_p_Mr = Mh + Mr

        return :(if i < $Mh_p_Mr
                $(dense_sub2ind_quote(Mh_p_Mr, N, offset, joffset))
            else
                $(dense_sub2ind_quote(Mh, N, offset + Mh_p_Mr * (N-1), joffset))
            end)
    else # we are splitting the matrix into two blocks, side by side.
        Nh, Nr = divrem(N, 2)
        Nh_p_Nr = Nh + Nr

        return :(if j < $Nh_p_Nr
                $(dense_sub2ind_quote(M, Nh_p_Nr, offset, joffset))
            else
                $(dense_sub2ind_quote(M, Nh, offset + M * Nh_p_Nr, joffset + Nh_p_Nr))
            end)
    end
end
@generated dense_sub2ind(::Val{M}, ::Val{N}, i, j) where {M,N} = dense_sub2ind_quote(M, N)

function triangle_sub2ind(N, i, j, offset = 0, joffset = 1)
    N < cutoff && return big_triangle(j-joffset) + i + offset
    
    Nh, Nr = divrem(N, 2)
    Nh_p_Nr = Nh + Nr
    if max(i,j) <= Nh_p_Nr # enter block 1 # triangle / symmetric block
        return triangle_sub2ind(Nh_p_Nr, i, j, offset, joffset)
    elseif i <= Nh_p_Nr # enter block 3 # square block
        return dense_sub2ind((Nh_p_Nr, Nh), i, j, offset + big_triangle(Nh_p_Nr), joffset + Nh_p_Nr)
    elseif j <= Nh_p_Nr # enter block 2 #symmetric or 0 block # Not supported; check i, j before calling.
        # $(dense_sub2ind_quote(Mh, Nh_p_Nr, offset + Mh_p_Mr * (Nh_p_Nr-1), joffset))
        throw("Indexing into this half not supported.")
    else # enter block 4 # triangle / symmetric block
        return triangle_sub2ind(Nh, i, j, offset + big_triangle(Nh_p_Nr) + Nh_p_Nr*(Nh-1), joffset + Nh_p_Nr)
    end
end


function triangle_sub2ind_quote(N, offset = 0, joffset = 1)
    N <= cutoff && return :( big_triangle(j-$joffset) + i + $offset )
    
    Nh, Nr = divrem(N, 2)
    Nh_p_Nr = Nh + Nr
    :(if max(i,j) <= $Nh_p_Nr # enter block 1 # triangle / symmetric block
            $(triangle_sub2ind_quote(Nh_p_Nr, offset, joffset))
        elseif i <= $Nh_p_Nr # enter block 3 # square block
            $(dense_sub2ind_quote(Nh_p_Nr, Nh, offset + big_triangle(Nh_p_Nr), joffset + Nh_p_Nr))
        # elseif j <= $Nh_p_Nr # enter block 2 #symmetric or 0 block # Not supported; check i, j before calling.
        #     # $(dense_sub2ind_quote(Mh, Nh_p_Nr, offset + Mh_p_Mr * (Nh_p_Nr-1), joffset))
        #     throw("Indexing into this half not supported.")
        else # enter block 4 # triangle / symmetric block
            $(triangle_sub2ind_quote(Nh, offset + big_triangle(Nh_p_Nr) + Nh_p_Nr*(Nh-1), joffset + Nh_p_Nr))
        end)
end

@generated function triangle_sub2ind(::Val{N}, i, j) where N
    quote
        i > j && throw("Indexing into this half not supported.")
        $(triangle_sub2ind_quote(N))
    end
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