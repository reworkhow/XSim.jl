function mate(cohort_common ::Cohort,
              cohort_pool   ::Cohort;
              n_common      ::Int64=cohort_common.n,
              n_pool        ::Int64=cohort_pool.n,
              n_per_mate    ::Int64=1,
              replace_common::Bool=false,
              replace_pool  ::Bool=false)

    # pre-allocate animals
    n_animals = n_common * n_pool * n_per_mate
    animals = Array{Animal}(undef, n_animals)

    idx = 1
    # sample common animals
    select_common = sample(cohort_common, n_common, replace=replace_common)
    for animal_common in select_common

        # sample pool animals
        select_pool = sample(cohort_pool, n_pool, replace=replace_pool)
        for animal_pool in select_pool

            # generate n progenies per mate from selected parents
            for i in 1:n_per_mate
                animals[idx] = Animal(animal_common, animal_pool)
                idx += 1
            end
        end
    end

    return Cohort(animals)
end


