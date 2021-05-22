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

    idx = 1
    # Sample common animals
    select_common = sample(cohort_common, n_common, replace=replace_common)
    for animal_common in select_common

        # Sample pool animals
        select_pool = sample(cohort_pool, n_pool, replace=replace_pool)
        for animal_pool in select_pool

            # Generate n progenies per mate from selected parents
            for _ in 1:n_per_mate
                animals[idx] = Animal(animal_common, animal_pool)
                idx += 1
            end
        end
    end
    cohort = Cohort(animals)

    # Offspring gender ratio
    if ratio_malefemale != -1.0
        float_males = round(n_animals * ratio_malefemale / (ratio_malefemale + 1))
        n_males     = convert(Int64, float_males)
        cohort      = cohort[1:n_males], cohort[(n_males + 1):end]
    end

    # Log
    if !silent
        LOG("--------- Mating Summary ---------")
        LOG("Generate $n_animals individuals from $n_common shared breeds")
        LOG("Every shared breeds mates with $n_pool breeds")
        LOG("And each mating reproduces $n_per_mate progenies")
        LOG("")
        LOG("--------- Offsprings Summary ---------")
        # print(cohort)
    end

    return cohort
end

mate(cohort::Cohort; args...) =  mate(cohort, cohort; args...)


# Base.Iterators.Pairs{Symbol, Any, NTuple{6, Symbol},
#                      NamedTuple{(:replace_common, 
#                                  :n_per_mate, 
#                                  :replace_pool, 
#                                  :ratio_malefemale, 
#                                  :h2, :is_random),
#                                  Tuple{Bool,Int64,Bool,Int64,Array{Float64,1},Bool}}}(:replace_common => false,:n_per_mate => 2,:replace_pool => false,:ratio_malefemale => 1,:h2 => [0.8, 0.2],:is_random => false)