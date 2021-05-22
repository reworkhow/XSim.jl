function sample_select(sires             ::Cohort,
                       dams              ::Cohort,
                       n                 ::Int64,
                       n_sires           ::Int64,
                       n_dams            ::Int64,
                       n_gen             ::Int64;
                       weights           ::Array{Float64, 1}=[1.0],
                       is_positive       ::Bool             =true)

    progeny_male, progeny_female = Cohort(), Cohort()
    pool_sires, pool_dams        = Cohort(sires.animals), Cohort(dams.animals)
    n_mates                      = round(Int, n / 2)

    for i in 1:n_gen
        # select individuals for sires and dams
        selected_sires = select(pool_sires, n_sires, weights=weights,
                                is_positive=is_positive)
        selected_dams  = select(pool_dams, n_dams, weights=weights,
                                is_positive=is_positive)
        # mate between selected individauls
        progeny_male   = random_mate(selected_sires, selected_dams, n_mates)
        progeny_female = random_mate(selected_sires, selected_dams, n_mates)
        # concat progenies into the breeding pool
        pool_sires += progeny_male
        pool_dams  += progeny_female
    end

    return progeny_male, progeny_female
end


function sample_random(sires ::Cohort,
                       dams  ::Cohort,
                       n     ::Int64,
                       n_gen ::Int64)

    progeny_male, progeny_female = Cohort(), Cohort()
    n_mate = round(Int, n / 2)

    for i in 1:n_gen
        progeny_male   = random_mate(sires, dams, n_mate)
        progeny_female = random_mate(sires, dams, n_mate)
        sires = progeny_male
        dams  = progeny_female
    end
    return progeny_male, progeny_female
end


function random_mate(sires::Cohort, dams::Cohort, n::Int64)
    return mate(sires, dams,
                n_common=n, n_pool=1, n_per_mate=1,
                replace_common=true, replace_pool=true)
end

function self_mate(cohort::Cohort, n::Int64)
    return random_mate(cohort, cohort, n)
end

function all_mate(sires::Animal, dams::Animal)
    return mate(sires, dams)
end


function embryo_transfer(dams::Animal, sires::Animal)
    return mate(dams, sires)
end

function build_demo()
    chromosome    = [1, 1, 2, 2, 2,
                     3, 3, 3, 4, 4]
    bp            = [20, 50, 10, 20, 30,
                     30, 50, 60, 20, 40]
    cM            = [25.0, 55.0, 15.0, 25.0, 35.0,
                     35.0, 55.0, 65.0, 25.0, 45.0]
    maf           = fill(0.5, 10)
    rate_mutation = 0.0
    rate_error    = 0.0
    build_genome(chromosome, bp, cM, maf, rate_mutation, rate_error)

    n_qtl = [2, 2]
    Vg    = [ 1 .5
             .5  1]
    build_phenome(n_qtl, Vg)

end