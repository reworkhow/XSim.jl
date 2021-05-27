function mate(cohort_shared      ::Cohort,
              cohort_per_shared  ::Cohort;
              n_shared           ::Int64=cohort_shared.n,
              n_per_shared       ::Int64=1,
              n_per_mate         ::Int64=1,
              n_pop              ::Int64=-1,
              replace_shared     ::Bool =false,
              replace_per_shared ::Bool =false,
              ratio_malefemale   ::Union{Float64, Int64}=-1.0,
              silent             ::Bool =GLOBAL("silent"),
              args...)

    # Preallocate animals, infer size
    animals, n_shared, n_per_shared, n_per_mate, n_pop =
        preallocate_animals(n_shared, n_per_shared, n_per_mate, n_pop)

    # Sample shared and per_shared breeds
    select_shared     = sample(cohort_shared, n_shared, replace=replace_shared)
    select_per_shared = sample_per_shared(cohort_per_shared,
                                          n_shared, n_per_shared,
                                          replace_per_shared)
    # Mating
    for i in 1:n_shared
        animal_shared = select_shared[i]

        for j in 1:n_per_shared
            animal_per_shared = select_per_shared[(i - 1) * n_per_shared + j]

            for k in 1:n_per_mate
                idx = (i - 1) * n_per_shared * n_per_mate + (j - 1) * n_per_mate + k
                animals[idx] = Animal(animal_shared, animal_per_shared)
            end
        end
    end

    cohort = split_by_ratio(animals, ratio_malefemale, n_pop)

    log_mate(silent, n_pop, n_shared, n_per_shared, n_per_mate)

    return cohort
end

mate(cohort::Cohort; args...) =  mate(cohort, cohort; args...)

function log_mate(silent, n_pop, n_shared, n_per_shared, n_per_mate)
    if !silent
        LOG("--------- Mating Summary ---------")
        LOG("Generate $n_pop individuals from $n_shared shared individuals")
        LOG("Every shared individuals mates with $n_per_shared individuals")
        LOG("And each mating reproduces $n_per_mate progenies")
        LOG("")
        LOG("--------- Offsprings Summary ---------")
        # print(cohort)
    end
end

function preallocate_animals(n_shared     ::Int64,
                             n_per_shared ::Int64,
                             n_per_mate   ::Int64,
                             n_pop        ::Int64)

    # If no population size assigned, the rest three are assumed to be correct
    if n_pop == -1
        n_pop = n_shared * n_per_shared * n_per_mate

    # If assigned population size, but doens't match the rest
    elseif n_pop != n_shared * n_per_shared * n_per_mate
        # If n_per_shared is not provided, and the rest are provided
        if n_per_shared == 1
            n_per_shared = trunc(Int, n_pop / n_shared / n_per_mate)

        # If n_per_mate is not provided, and the rest are provided
        elseif n_per_mate == 1
            n_per_mate   = trunc(Int, n_pop / n_shared / n_per_shared)
        end
    end

    return Array{Animal}(undef, n_pop), n_shared, n_per_shared, n_per_mate, n_pop
end

function sample_per_shared(cohort_per_shared,
                           n_shared,
                           n_per_shared,
                           replace_per_shared)

    # Handle 'more samples without replacement' error, assumed to mate all
    if (!replace_per_shared) && (n_shared * n_per_shared > cohort_per_shared.n)
       select_per_shared = Cohort(repeat(cohort_per_shared.animals, outer=n_shared))

    else
       select_per_shared = sample(cohort_per_shared, n_shared * n_per_shared,
                                  replace=replace_per_shared)
   end

   return select_per_shared
end

function split_by_ratio(animals,
                        ratio_malefemale,
                        n_pop)

    cohort = Cohort(animals)
    if ratio_malefemale != -1.0
        float_males = round(n_pop * ratio_malefemale / (ratio_malefemale + 1))
        n_males     = convert(Int64, float_males)
        return cohort[1:n_males], cohort[(n_males + 1):end]

    else
        return cohort
    end
end

