function build_demo()
    n_chr = 10
    n_loci_chr = 100
    n_loci = n_chr * n_loci_chr

    chromosome = [i         for i in 1:n_chr for j in 1:n_loci_chr]
    bp         = [10 * j    for i in 1:n_chr for j in 1:n_loci_chr]
    cM         = [1.5 * j  for i in 1:n_chr for j in 1:n_loci_chr]
    maf        = fill(0.5, n_loci)
    rate_mutation = 0.0
    rate_error    = 0.0
    build_genome(chromosome, bp, cM, maf;
                 rate_mutation=rate_mutation,
                 rate_error=rate_error)

    n_qtl = [3, 8]
    vg    = [1.0  0
             0  1.0]
    build_phenome(n_qtl; vg=vg)

end

function build_demo_small()
    n_chr = 2
    n_loci_chr = 5
    n_loci = n_chr * n_loci_chr

    chromosome = [i         for i in 1:n_chr for j in 1:n_loci_chr]
    bp         = [10 * j    for i in 1:n_chr for j in 1:n_loci_chr]
    cM         = [1.5 * j  for i in 1:n_chr for j in 1:n_loci_chr]
    maf        = fill(0.5, n_loci)
    rate_mutation = 0.0
    rate_error    = 0.0
    build_genome(chromosome, bp, cM, maf;
                 rate_mutation=rate_mutation,
                 rate_error=rate_error)

    n_qtl = [2, 4]
    vg    = [1.0  0
             0  1.0]
    build_phenome(n_qtl; vg=vg)

end

function breed(sires            ::Cohort,
               dams             ::Cohort;
               n_gens           ::Int64=1,
               n_select         ::Int64=sires.n + dams.n,
               n_select_males   ::Int64=sires.n,
               n_select_females ::Int64=dams.n,
               args...)

    # LOG for parents
    sm = summary(sires + dams)
    LOG("Gen 0 -> Mean of BVs: $(sm["mu_g"]), Variance of BVs: $(sm["var_g"])")

    for i in 1:n_gens

        # If assigned male, female ratio
        if haskey(args, :ratio_malefemale)
            # Mate
            males, females = mate(sires, dams; silent=true, args...)

            # Select on males
            if n_select_males > 0
                sires = select(males, n_select_males; silent=true, args...)
            else
                sires = males
            end

            # Select on females
            if n_select_females > 0
                dams = select(females, n_select_females; silent=true, args...)
            else
                dams = females
            end

        else
            # Mate
            progenies = mate(sires, dams; silent=true, args...)

            # Select
            if n_select > 0
                progenies = select(progenies, n_select; silent=true, args...)
            end

            # Define parents for next round
            sires, dams = progenies, progenies
        end

        # LOG
        sm = summary(sires + dams)
        LOG("Gen $i -> Mean of BVs: $(sm["mu_g"]), Variance of BVs: $(sm["var_g"])")
    end

    return sires, dams
end

breed(cohort::Cohort, n_gens::Int64, args...) =
breed(cohort, cohort, n_gens; args...)


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
                n_shared=n, n_per_shared=1, n_per_mate=1,
                replace_shared=true, replace_per_shared=true)
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
