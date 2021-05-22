function build_demo()
    n_chr = 10
    n_loci_chr = 5
    n_loci = n_chr * n_loci_chr

    chromosome = [i         for i in 1:n_chr for j in 1:n_loci_chr]
    bp         = [10 * j    for i in 1:n_chr for j in 1:n_loci_chr]
    cM         = [15.0 * j  for i in 1:n_chr for j in 1:n_loci_chr]
    maf        = fill(0.5, n_loci)
    rate_mutation = 0.0
    rate_error    = 0.0
    build_genome(chromosome, bp, cM, maf, rate_mutation, rate_error)

    n_qtl = [3, 8]
    Vg    = [ 1 .5
             .5  1]
    build_phenome(n_qtl, Vg)

end

function simple_breed(sires  ::Cohort,
                      dams   ::Cohort,
                      n_gens ::Int64;
                      n_sires::Int64=sires.n,
                      n_dams ::Int64=dams.n,
                      args...)

    # LOG for parents
    sm = summary(sires + dams)
    LOG("Gen 0 -> Mean of BVs: $(sm["mu_g"]), Variance of BVs: $(sm["var_g"])")
    for i in 1:n_gens
        # Mate and select
        SILENT(true)
        males, females = mate(sires, dams; args...)
        sires          = select(males,   n_sires; args...)
        dams           = select(females, n_dams ; args...)
        SILENT(false)
        # LOG
        sm             = summary(sires + dams)
        LOG("Gen $i -> Mean of BVs: $(sm["mu_g"]), Variance of BVs: $(sm["var_g"])")
    end

    return sires, dams
end

simple_breed(cohort::Cohort, n_gens::Int64, args...) =
simple_breed(cohort, cohort, n_gens; args...)


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
