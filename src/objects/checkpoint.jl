mutable struct CheckPoint
    # n is the number of samples
    # t is the number of traits

    bv   ::Array, # breeding values: n by t array
    vg   ::Array, # genetic variance: t by t matrix
    cost ::Float64, # cost in total
    n    ::Int,    # number of individuals
    t    ::Int     # number of traits

    # Base constructor
    function CheckPoint()
        instance    = new()
        return instance
    end

    # criteria
    ## GS pool marker acc
    ## Breeding cycle time
    ## DH population performance (BV, Vg)
    ## cost per DH line
end

function update!(ck      ::CheckPoint,
                 cohort  ::Cohort,
                 effects ::Array,
                 cost    ::Float64)

    # obtain breeding values
    bvs     = cohort |> get_bv
    bvs_avg = mean(bvs, dims=2) # breeding value by traits
    ck.bv   = bvs_avg
    # obtain genetic variance
    frequency = get_MAF(cohort)
    ck.vg     = get_Vg(effects, frequency)
    # update parameters
    m, t    = size(effects)
    ck.n    = cohort.n
    ck.t    = t
    # TODO: how to estimate cost?
    ck.cost = cost
end
