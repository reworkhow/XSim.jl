function mate(cohort_common   ::Cohort,
              cohort_pool     ::Cohort;
              n_common        ::Int64                =cohort_common.n,
              n_pool          ::Int64                =cohort_pool.n,
              n_per_mate      ::Int64                =1,
              replace_common  ::Bool                 =false,
              replace_pool    ::Bool                 =false,
              ratio_malefemale::Union{Float64, Int64}=-1.0,
              silent          ::Bool                 =GLOBAL("silent"),
              args...)

    # Pre-allocate animals
    n_animals = n_common * n_pool * n_per_mate
    animals = Array{Animal}(undef, n_animals)

    # Sample common and pool breeds
    select_common = sample(cohort_common, n_common,          replace=replace_common)
    select_pool   = sample(cohort_pool,   n_common * n_pool, replace=replace_pool)

    # Mating
    for i in 1:n_common
        animal_common = select_common[i]

        for j in 1:n_pool
            animal_pool = select_pool[(i - 1) * n_pool + j]

            for k in 1:n_per_mate
                idx = (i - 1) * n_pool * n_per_mate + (j - 1) * n_per_mate + k
                animals[idx] = Animal(animal_common, animal_pool)
            end
        end
    end
    # print(animals)
    cohort = Cohort(animals)

    # Offspring gender ratio
    if ratio_malefemale != -1.0
        float_males = round(n_animals * ratio_malefemale / (ratio_malefemale + 1))
        n_males     = convert(Int64, float_males)
        cohort      = cohort[1:n_males], cohort[(n_males + 1):end]
    end
    # Log
    log_mate(silent, n_animals, n_common, n_pool, n_per_mate)

    return cohort
end

mate(cohort::Cohort; args...) =  mate(cohort, cohort; args...)

function log_mate(silent, n_animals, n_common, n_pool, n_per_mate)
    if !silent
        LOG("--------- Mating Summary ---------")
        LOG("Generate $n_animals individuals from $n_common shared breeds")
        LOG("Every shared breeds mates with $n_pool breeds")
        LOG("And each mating reproduces $n_per_mate progenies")
        LOG("")
        LOG("--------- Offsprings Summary ---------")
        # print(cohort)
    end
end