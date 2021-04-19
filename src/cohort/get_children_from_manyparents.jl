function mate(cohort_common ::Cohort,
              cohort_pool   ::Cohort;
              n_common      ::Int16=cohort_common.n,
              n_pool        ::Int16=cohort_pool.n,
              n_per_mate    ::Int16=1,
              replace_common::Bool=false,
              replace_pool  ::Bool=false)

    n_animals = n_common * n_pool * n_per_mate
    animals = Array{Animal}(undef, n_animals)
    println("Get ", n_animals, " offspring into the next generation.")

    # Dicts used to count in how many matings an indiviudal is used
    dict_dam  = Dict{Int64,Int64}()
    dict_sire = Dict{Int64,Int64}()

    idx = 1
    # select animal as common breed
    select_common = sample(cohort_common.n, n_common, replace=replace_common)
    for i in 1:n_common
        animal_common = cohort_common.animals[select_common[i]]

        # select animal from pool
        select_pool = sample(cohort_pool.n, n_pool, replace=replace_pool)
        for j in 1:n_pool
            animal_pool = cohort_pool.animals[select_pool[j]]

            # generate progenies
            for k in 1:n_per_mate
                animals[idx] = get_child(animal_common, animal_pool)
                dict_dam[animal_common.ID] = get(dict_dam, animal_common.ID, 0) + 1
                dict_sire[animal_pool.ID]  = get(dict_sire, animal_pool.ID, 0) + 1
                idx += 1
            end
        end
    end
    println("Number of fathers used: ", length(dict_dam))
    println("Number of mothers used: ", length(dict_sire))

    return Cohort(animals)
end

function sample(pool::Int, n::Int; replace::Bool=false)
    # select n samples from 1:pool numbers, return a size of n 1-d array
    if replace
        samples = rand(1:pool, n)
    else
        n = n > pool ? pool : n
        samples = shuffle(1:pool)[1:n]
    end

    return n == 1 ? samples[1] : samples
end
