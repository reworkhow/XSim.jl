mutable struct Checkpoint
    # t is the number of traits
    bv::Array{Float64,1} # average BV: t x 1 vector
    ebv::Array{Float64,1} # average EBV: t x 1 vector
    vg::Array{Float64,2} # genetic variance: t by t matrix
    evg::Array{Float64,2} # genetic variance: t by t matrix
    cost_gp::Float64 # cost of genotyping
    cost_dh::Float64 # cost of generating DHs
    time::Float64 # number of seasons
    ## TODO: LD DECAY O(p^2), p = 3k
    cost_gp_per_line::Float64 # cost of genotyping per line
    cost_dh_per_line::Float64 # cost of generating DHs per line

    function Checkpoint()
        instance = new()
        instance.cost_gp_per_line = 10.0
        instance.cost_dh_per_line = 10.0
        return instance
    end
end

function update_bv!(checkpoint::Checkpoint, cohort::Cohort)
    bvs = get_bv(cohort)
    checkpoint.bv = mean(bvs, dims=2) # breeding value by traits
end

function update_bv!(checkpoint::Checkpoint, gs_pool::GS_pool)
    update_bv!(checkpoint, gs_pool.cohort)
end

function update_vg!(checkpoint::Checkpoint,
    cohort::Cohort, effects::Array)
    checkpoint.vg = get_vg(cohort, effects)
end

function update_vg!(checkpoint::Checkpoint, gs_pool::GS_pool)
    update_vg!(checkpoint, gs_pool.cohort, gs_pool.effects)
end

function update_cost!(checkpoint::Checkpoint, n_gp::Int, n_dh::Int)
    checkpoint.cost_gp += n_gp * checkpoint.cost_gp_per_line
    checkpoint.cost_dh += n_dh * checkpoint.cost_dh_per_line
end

function update_time!(checkpoint::Checkpoint, time::Int)
    checkpoint.time = time
end


mutable struct CheckpointList

    checkpoints::Array{Checkpoint,1}

    function CheckpointList()
        checkpoints = Checkpoint[]
        instance = new(checkpoints)
        return instance
    end
end

function update!(ckl::CheckpointList,
    gs_pool::GS_pool=nothing;
    n_gp::Int=0, # number of new genotyped individuals
    n_dh::Int=0) # number of new DHs

    # instantiate a new checkpoint
    ck = Checkpoint()
    # get the last checkpoint
    last_ck = ckl.checkpoints[end]
    # update the checkpoint
    if isnothing(gs_pool)
        ck = copy(last_ck)
    else
        update_bv!(ck, gs_pool)
        update_vg!(ck, gs_pool)
    end
    update_cost!(ck, n_gp, n_dh)
    update_time!(ck, last_ck.time + 1)
    # add the checkpoint to the list
    push!(ckl.checkpoints, ck)
    # TODO: estimated bv and vg
end

function Base.:+(ckl1::CheckpointList, ckl2::CheckpointList)
    ckl1.checkpoints = vcat(ckl1.checkpoints, ckl2.checkpoints)
    return ckl1
end