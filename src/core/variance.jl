"""
Turn (n, ) to (n, 1) or a scaler to a 2D array```
"""
function matrix(inputs::Any; is_sparse=false)
    mat = hcat(Diagonal([inputs])...)
    if is_sparse
        return sparse(mat)
    else
        return mat
    end
end

"""
### Test
```jldoctest
    handle_diagonal([1, 5], 1) # error: size does not match
    handle_diagonal([1, 5], 2) # [1 0; 0 5]
    handle_diagonal(3, 1) # [[3]]
    handle_diagonal(3, 2) # [3 0; 0 3]
    handle_diagonal([3 6; 3 6], 2) # [3 6; 3 6]
```
"""
function handle_diagonal(inputs::Union{Array,Float64,Int64},
                         n_traits::Int64)

    # Cast variants of variances to a 2-D array
    # Case 1 When variances is a scaler, assign it as the diagonal of variances
    # if !isa(inputs, Array)
    if length(inputs) == 1
        inputs = diagm(fill(inputs[1], n_traits)) # [1] handle length 1 vector
    else
        inputs = matrix(inputs)
        # Case 2 When variances is a vector (ncol == 1),
        # assign it as the diagonal of a variance matrix
        if size(inputs)[2] == 1
            inputs = diagm(inputs[:, 1])
        end
    end

    if size(inputs)[2] != n_traits
        LOG("Dimensions don't match between n_traits and variances/h2", "error")
    end

    return inputs
end

function get_Vg(QTL_effects::Union{Array{Float64,2},SparseMatrixCSC},
                QTL_freq::Array{Float64,1})
    # Falconer and Mackay, 1996, Equation (8.3b)
    # 2pq
    D = diagm(2 * QTL_freq .* (1 .- QTL_freq))

    # 2pq*alpha^2
    Vg = QTL_effects'D * QTL_effects

    return Vg
end

function get_Ve(n_traits::Int64,
                Vg      ::Union{Array{Float64},Float64},
                h2      ::Union{Array{Float64},Float64}=0.5)

    h2 = handle_h2(h2, n_traits)
    Vg = handle_diagonal(Vg, n_traits)

    Ve = ((ones(n_traits) .- h2) .* diag(Vg)) ./ h2 # diagm can't handle 2x2 matrix
    Ve = n_traits == 1 ? Ve[1] : Ve

    return handle_diagonal(Ve, n_traits)
end


function sample_qtls(n_qtls, n_traits, vg)
    sampler = MvNormal([0 for _ in 1:n_traits], vg)
    qtl_effects = rand(sampler, n_qtls)
    return matrix(qtl_effects')
end

"""
### Test
```jldoctest
    handle_h2(missing, 1)    # [0.5]
    handle_h2(missing, 2)    # [0.5, 0.5]
    handle_h2(0.3, 1)        # [0.3]
    handle_h2([0.3], 1)      # [0.3]
    handle_h2([0.3, 0.2], 1) # [0.3]
    handle_h2(0.3, 2)        # [0.3, 0.3]
    handle_h2([0.3], 2)      # [0.3, 0.3]
    handle_h2([0.3, 0.2], 2) # [0.3, 0.2]
```
"""
function handle_h2(h2, n_traits)
    if ismissing(h2)
        # no h2 is provided
        h2 = [0.5 for _ in 1:n_traits]
    else
        # make sure h2 is an array when n_traits == 1
        h2 = ifelse(isa(h2, Array), h2, [h2])
        # make sure h2 vector is of length n_traits
        h2 = ifelse(length(h2) == n_traits, h2, [h2[1] for _ in 1:n_traits])
    end

    # avoid inf variance when h2 = 0
    is_zeros = h2 .== 0
    if n_traits > 1
        h2[is_zeros] .= 1e-5
    elseif is_zeros == true # single trait and vg == 0
        h2[is_zeros] = 1e-5
    end

    return h2
end

function handle_variance(vg, ve, vp, h2)
    # check inputs
    n_traits = length(h2)
    has_vg = !ismissing(vg)
    has_vp = !ismissing(vp)
    has_ve = !ismissing(ve)
    bool_inputs = [has_vg, has_vp, has_ve]

    # cases
    if sum(bool_inputs) == 0
        vg = handle_diagonal([1.0], n_traits)
        ve = infer_variances(vg, h2=h2, term_src="vg")

    elseif sum(bool_inputs) == 1
        if has_vg
            ve = infer_variances(vg, h2=h2, term_src="vg")
        elseif has_vp
            vg = infer_variances(vp, h2=h2, term_src="vp")
            ve = vp - vg
        elseif has_ve
            vg = infer_variances(ve, h2=h2, term_src="ve")
        end

    elseif sum(bool_inputs) == 2
        if has_vg && has_ve
            nothing
        elseif has_vg && has_vp
            ve = vp - vg
        elseif has_ve && has_vp
            vg = vp - ve
        end
        # estimate h2 based on the provided vg and ve
        h2 = diag(vg / (vg + ve))
        h2 = round.(h2, digits=3)
    end

    return vg, ve, vp, h2
end

function infer_variances(v_src;
                         h2::Array{Float64},
                         term_src::String)

    n_traits = length(h2)
    v_src = handle_diagonal(v_src, n_traits)

    if term_src == "vg"
        # out must be ve
        # g->e: e = (1 - h2) g / h2
        v_out = ((ones(n_traits) .- h2) .* diag(v_src)) ./ h2

    elseif term_src == "ve"
        # out must be vg
        # e->g: g = e * h2 / (1 - h2)
        v_out = (diag(v_src) .* h2) ./ (ones(n_traits) .- h2)

    elseif term_src == "vp"
        # out must be vg
        # p->g: g = p * h2
        v_out = diag(v_src) .* h2
    end

    v_out = n_traits == 1 ? v_out[1] : v_out

    return handle_diagonal(v_out, n_traits)
end

