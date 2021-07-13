function mate(cohort_A         ::Cohort,
              cohort_B         ::Cohort;
              nA               ::Int64=cohort_A.n,
              nB_per_A         ::Int64=1,
              n_per_mate       ::Int64=1,
              n_pop            ::Int64=-1,
              replace_A        ::Bool =false,
              replace_B        ::Bool =false,
              ratio_malefemale ::Union{Float64, Int64}=1, # 0 is no genders
              is_selfing       ::Bool =false,
              silent           ::Bool =GLOBAL("silent"),
              scheme           ::String ="none",
              args...)

    if scheme != "none"
        if scheme == "random"
            return mate(cohort_A, cohort_B)

        elseif scheme == "diallel cross"
            return mate(cohort_A, cohort_B; nA=cohort_A.n, nB_per_A=cohort_B.n)

        else
            LOG("The available scheme options are: [random, diallel cross]", "error")
        end

    else
        if is_selfing
            # Preallocate animals
            n_pop   = nA * n_per_mate
            animals = Array{Animal}(undef, n_pop)
            # Sample parent from cohort_A
            select_A = sample(cohort_A, nA, replace=replace_A)
            # Mating
            for i in 1:nA
                parent = select_A[i]
                for j in 1:n_per_mate
                    idx = (i - 1) * n_per_mate + j
                    animals[idx] = Animal(parent, parent)
                end
            end
            cohort = Cohort(animals)

        else
            # Preallocate animals, infer size
            animals, nA, nB_per_A, n_per_mate, n_pop =
                preallocate_animals(nA, nB_per_A, n_per_mate, n_pop)
            # Sample A and B breeds
            select_A = sample(  cohort_A, nA, replace=replace_A)
            select_B = sample_B(cohort_B, nA, nB_per_A, replace_B)
            # Mating
            for i in 1:nA
                animal_A = select_A[i]

                for j in 1:nB_per_A
                    animal_B = select_B[(i - 1) * nB_per_A + j]

                    for k in 1:n_per_mate
                        idx = (i - 1) * nB_per_A * n_per_mate + (j - 1) * n_per_mate + k
                        animals[idx] = Animal(animal_A, animal_B)
                    end
                end
            end
            cohort = split_by_ratio(animals, ratio_malefemale, n_pop)

        end
        log_mate(silent, n_pop, nA, nB_per_A, n_per_mate, is_selfing)
        return cohort
    end
end

mate(cohort::Cohort; args...) =  mate(cohort, cohort; args...)

function preallocate_animals(nA        ::Int64,
                             nB_per_A  ::Int64,
                             n_per_mate::Int64,
                             n_pop     ::Int64)

    # If no population size assigned, the rest three are assumed to be correct
    if n_pop == -1
        n_pop = nA * nB_per_A * n_per_mate

    # If assigned population size, but doens't match the rest
    elseif n_pop != nA * nB_per_A * n_per_mate
        # If nB_per_A is not provided, and the rest are provided
        if nB_per_A == 1
            nB_per_A = trunc(Int, n_pop / nA / n_per_mate) # ground to int

        # If n_per_mate is not provided, and the rest are provided
        elseif n_per_mate == 1
            n_per_mate   = trunc(Int, n_pop / nA / nB_per_A)
        end
    end

    return Array{Animal}(undef, n_pop), nA, nB_per_A, n_per_mate, n_pop
end

function sample_B(cohort_B,
                  nA,
                  nB_per_A,
                  replace_B)

    n_crosses = nA * nB_per_A

    if n_crosses > cohort_B.n
        if replace_B
            # > 10x10, 5x3, T, T
            select_B = sample(cohort_B, n_crosses, replace=replace_B)

        else
            # > 10x10, 5x3, T, F
            select_B = Cohort()
            for _ in 1:nA
                select_B += sample(cohort_B, nB_per_A, replace=replace_B)
            end

        end
    else
        # > 10x10, 5x2, T, T    # > 10x10, 5x2, T, F
        select_B = sample(cohort_B, n_crosses, replace=replace_B)
    end

    return select_B

    # NOTE
    # > 10x10, 5x2, T, T
    # (10, 10, T)

    # > 10x10, 5x2, T, F
    # (10, 10, F)

    # > 10x10, 5x3, T, T
    # all crosses = 15
    # (10, 15, T)

    # > 10x10, 5x3, T, F
    # all crosses = 15
    # (10, 3, F)* 5
end

function split_by_ratio(animals,
                        ratio_malefemale,
                        n_pop)

    cohort = Cohort(animals)
    # shuffle the order
    cohort = sample(cohort, cohort.n, replace=false)

    if ratio_malefemale != 0
        float_males = round(n_pop * ratio_malefemale / (ratio_malefemale + 1))
        n_males     = convert(Int64, float_males)
        return cohort[1:n_males], cohort[(n_males + 1):end]

    else
        return cohort
    end
end

function log_mate(silent, n_pop, nA, nB_per_A, n_per_mate, is_selfing)
    if !silent
        if is_selfing
            LOG("--------- Mating Summary ---------")
            LOG("Generated $n_pop individuals from $nA cohort_A individuals")
            LOG("Each selfing parent reproduces $n_per_mate progenies")
            LOG("")
            LOG("--------- Offsprings Summary ---------")

        else
            LOG("--------- Mating Summary ---------")
            LOG("Generated $n_pop individuals from $nA cohort_A individuals")
            LOG("Every cohort_A individual mates with $nB_per_A individuals from cohort_B")
            LOG("And each mating reproduces $n_per_mate progenies")
            LOG("")
            LOG("--------- Offsprings Summary ---------")
        end
        # print(cohort)
    end
end

Base.:*(x::Cohort, y::Cohort) = mate(x, y)
