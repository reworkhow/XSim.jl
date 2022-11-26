mutable struct Checkpoint
    # t is the number of traits
    bv_p::Array{Float64,1} #  BV: t x 1 vector
    vg_p::Array{Float64,2} # genetic variance: t by t matrix
    bv_dh::Array{Float64,1} #  BV: t x 1 vector
    vg_dh::Array{Float64,2} # genetic variance: t by t matrix
    cost_gp_per_line::Float64 # cost of genotyping per line
    cost_dh_per_line::Float64 # cost of generating DHs per line
    cost_gp::Float64 # cost of genotyping
    cost_dh::Float64 # cost of generating DHs
    # time information
    cycle::Int64 # breeding cycle
    season::Int64 # seasons in the cycle (cycle=2, season=3, actual time=2+3-1=4)
    time::Int64 # time
    ## TODO: LD DECAY O(p^2), p = 3k

    """
    c: which cycle
    s: which season
    """
    function Checkpoint(c=0, s=0)
        instance = new()
        # bv and vg
        instance.bv_p = zeros(1)
        instance.vg_p = zeros(1, 1)
        instance.bv_dh = zeros(1)
        instance.vg_dh = zeros(1, 1)
        # cost
        instance.cost_gp_per_line = 10.0
        instance.cost_dh_per_line = 10.0
        instance.cost_gp = 0.0
        instance.cost_dh = 0.0
        # time
        instance.cycle = c
        instance.season = s
        instance.time = c + s - 1
        # return
        return instance
    end
end

function update!(
    checkpoint::Checkpoint,
    cohort::Cohort,
    effects::SparseMatrixCSC;
    option::String)

    update_bv!(checkpoint, cohort; option=option)
    update_vg!(checkpoint, cohort, effects; option=option)
end

function update_bv!(
    checkpoint::Checkpoint,
    cohort::Cohort;
    option::String)

    bvs = get_BVs(cohort)
    bvs_avg = mean(bvs, dims=1)[:, 1] # (n, 1) -> (n,)
    if option == "parents"
        checkpoint.bv_p = bvs_avg
    elseif option == "products"
        checkpoint.bv_dh = bvs_avg
    end
end

function update_vg!(
    checkpoint::Checkpoint,
    cohort::Cohort,
    effects::SparseMatrixCSC;
    option::String)

    vg = get_Vg(cohort, effects)
    if option == "parents"
        checkpoint.vg_p = vg
    elseif option == "products"
        checkpoint.vg_dh = vg
    end
end

function cost!(checkpoint::Checkpoint; n_gp::Int=0, n_dh::Int=0)
    checkpoint.cost_gp += n_gp * checkpoint.cost_gp_per_line
    checkpoint.cost_dh += n_dh * checkpoint.cost_dh_per_line
end



# function update_bv!(checkpoint::Checkpoint, gs_pool::GS_pool)
#     update_bv!(checkpoint, gs_pool.cohort)
# end

# function update_vg!(checkpoint::Checkpoint, gs_pool::GS_pool)
#     update_vg!(checkpoint, gs_pool.cohort, gs_pool.effects)
# end

function update_season!(checkpoint::Checkpoint, season::Int)
    checkpoint.season = season
end


mutable struct CheckpointList

    checkpoints::Array{Checkpoint,1}

    function CheckpointList(c=0, s=0)
        checkpoints = Checkpoint[]
        instance = new(checkpoints)
        # push a starting checkpoint
        push!(instance, Checkpoint(c, s))
        # return
        return instance
    end
end

function print(ckl::CheckpointList)
    nck = size(ckl.checkpoints)[1]
    println("CheckpointList: ")
    println("(t, s, c): parents, products, costs")
    for i in 1:nck
        ck = ckl.checkpoints[i]
        t, s, c = ck.time, ck.season, ck.cycle
        bvp  = round.(ck.bv_p, digits=4)[1]
        vgp  = round.(mean(diag(ck.vg_p)), digits=4)
        bvdh = round.(ck.bv_dh, digits=4)[1]
        vgdh = round.(mean(diag(ck.vg_dh)), digits=4)
        cgp  = round.(ck.cost_gp, digits=4)
        cgd  = round.(ck.cost_dh, digits=4)
        println("($t, $s, $c): [$(bvp), $(vgp)], [$(bvdh), $(vgdh)], [$(cgp), $(cgd)]")
    end
end

Base.show(io::IO, ckl::CheckpointList) = print(ckl)



function update!(
    ckl::CheckpointList;
    n_gp::Int=0, # number of new genotyped individuals
    n_dh::Int=0, # number of new DHs
    new_cohort::Cohort=Cohort(),
    h2::Float64=0.5)

    update!(ckl, GS_pool(); n_gp=n_gp, n_dh=n_dh, new_cohort=new_cohort, h2=h2)
end

function update!(
    ckl::CheckpointList,
    gs_pool::GS_pool;
    n_gp::Int=0, # number of new genotyped individuals
    n_dh::Int=0, # number of new DHs
    new_cohort::Cohort=Cohort(),
    h2::Float64=0.5)

    # instantiate a new checkpoint
    ck = Checkpoint()
    # get the last checkpoint
    last_ck = ckl.checkpoints[end]
    # GS pool
    if isempty(gs_pool)
        # duplicate the last checkpoint
        ck = deepcopy(last_ck)
    else
        # update the pool using new cohort
        update!(gs_pool, new_cohort; h2=h2)
        # update bv and vg
        update_bv!(ck, gs_pool)
        update_vg!(ck, gs_pool)
    end
    update_cost!(ck, n_gp, n_dh)
    update_season!(ck, last_ck.season + 1)
    # add the checkpoint to the list
    push!(ckl, ck)
    # TODO: estimated bv and vg
end

# ------------------------------

import Base.length
function length(ckl::CheckpointList)
    return length(ckl.checkpoints)
end

function isempty(ckl::CheckpointList)
    return isempty(ckl.checkpoints)
end

function push!(ckl::CheckpointList, ck::Checkpoint)
    push!(ckl.checkpoints, ck)
end

function push!(ckl_ori::CheckpointList, ckl_new::CheckpointList)
    for ck in ckl_new.checkpoints
        push!(ckl_ori, ck)
    end
end

function Base.:+(ckl1::CheckpointList, ckl2::CheckpointList)
    ckl1.checkpoints = vcat(ckl1.checkpoints, ckl2.checkpoints)
    return ckl1
    # TODO
    # sum the cost and merge the time
end

