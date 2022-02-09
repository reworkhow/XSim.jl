"""
# Wrap-up function: breed()

    breed(cohort_A         ::Cohort,
          cohort_B         ::Cohort;
          n_gens           ::Int64=1,
          n_select         ::Int64=cohort_A.n + cohort_B.n,
          n_select_A   ::Int64=cohort_A.n,
          n_select_B ::Int64=cohort_B.n,
          select_all_gens  ::Bool=false,
          args...)

    breed(cohort::Cohort, n_gens::Int64, args...) =
        breed(cohort, cohort, n_gens; args...)

### Arguments
Positional arguments
- `cohort_A` : A `cohort` object that is treated as common mating parents. It's asssumed to be sires/males parents.
- `cohort_B` : A `cohort` object that is a mating pool from which individuals are sampled to mate with `cohort_A`. It's assumed to be dams/female parents.

Keywords arguments
- `n_gens` : An integer specifying how many mate-select generations are carried out.
- `n_select` : Used when `ratio_malefemale` is set to `0`. `n_select` individuals are selected as parents for the next generation.
- `n_select_A` : Used when `ratio_malefemale` is not `0`. `n_select_A` will be selected as male parents for the next generation.
- `n_select_B` : Used when `ratio_malefemale` is not `0`. `n_select_B` will be selected as female parents for the next generation.
- `select_all_gens` : Default "false" and parents are not included in the next generation pool for selections. Set `select_all_gens` to "true" if the selections consider individuals from all generations.

### Outputs
By default, two `cohort` objects will be returned. The first `cohort` is assumed to be male progenies and the other `cohort` are female progenies. The size of two cohorts will folow the ratio `raiot_malefemale`. When `ratio_malefemale` is set to `0`, only one `cohort` will be returned.

### Examples
We can have `10` sires and mate each sire with `5` dams for `5` generations. In each generation, we randomly select `10` male progenies as sires and all female progenies as dams for the next generation. We can derive such breeding scheme as below:

```jldoctest
julia> args  = Dict(# mating arguments
                    :nA               => 10,
                    :nB_per_A         => 5,
                    :n_per_mate       => 2,
                    :ratio_malefemale => 1.0,
                    # selection arguments
                    :is_random        => true,
                    # breeding arguments
                    :n_gens           => 5,
                    :nA_select        => 10,
                    :select_all_gens  => true)
```
Simulate 10 sires and 50 dams as founders.
```jldoctest
julia> sires = Founders(10)
julia> dams  = Founders(50)
```
Breed cohorts based on the defined arguments.

```jldoctest
julia> sires, dams = breed(sires, dams; args...)
[ Info: Gen 0 -> Mean of BVs: [1.665 2.745], Variance of BVs: [1.008 0.479]
[ Info: Gen 1 -> Mean of BVs: [1.719 2.715], Variance of BVs: [0.96 0.546]
[ Info: Gen 2 -> Mean of BVs: [1.754 2.777], Variance of BVs: [0.936 0.724]
[ Info: Gen 3 -> Mean of BVs: [1.766 2.796], Variance of BVs: [0.927 0.775]
[ Info: Gen 4 -> Mean of BVs: [1.773 2.771], Variance of BVs: [0.951 0.781]
[ Info: Gen 5 -> Mean of BVs: [1.813 2.751], Variance of BVs: [0.952 0.804]
([ Info: Cohort (60 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.855 2.593]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.884 0.697]
, [ Info: Cohort (300 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.805 2.782]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.969 0.822]
```

The result is equivalent to the following mate-select iterations:
```jldoctest
julia> for _ in 1:5
            males, females = mate(sires, dams, args...)
            sires         += select(males, n_sires, args...)
            dams          += females
       end
```

"""
function breed(cohort_A         ::Cohort,
               cohort_B         ::Cohort;
               n_gens           ::Int64=1,
               n_select         ::Int64=cohort_A.n + cohort_B.n,
               n_select_A       ::Int64=cohort_A.n,
               n_select_B       ::Int64=cohort_B.n,
               select_all_gens  ::Bool=false,
               ratio_malefemale ::Union{Float64, Int64}=0,
               args...)

    # LOG for parents
    sm = summary(cohort_A + cohort_B)
    LOG("Gen 0 -> Mean of BVs: $(sm["mu_g"]), Variance of BVs: $(sm["var_g"])")

    for i in 1:n_gens
        # If assigned male, female ratio
        if ratio_malefemale != 0
            # Mate
            males, females = mate(cohort_A, cohort_B; silent=true, ratio_malefemale=ratio_malefemale, args...)
            print(males)
            print(females)
            # Select on males
            if n_select_A > 0
                males_s = select(males, n_select_A; silent=true, criteria="random", args...)
            else
                males_s = males
            end

            # Select on females
            if n_select_B > 0
                females_s = select(females, n_select_B; silent=true, criteria="random", args...)
            else
                females_s = females
            end

            # Define parents for next round
            if select_all_gens
                cohort_A += males_s
                cohort_B += females_s
            else
                cohort_A, cohort_B = males_s, females_s
            end

        else
            # Mate
            progenies = mate(cohort_A, cohort_B; silent=true, args...)
            # Select
            if n_select > 0
                progenies = select(progenies, n_select; silent=true, criteria="random", args...)
            end

            # Define parents for next round
            if select_all_gens
                cohort_A += progenies
                cohort_B += progenies
            else
                cohort_A = progenies
                cohort_B = progenies
            end

        end

        if length(cohort_A) == 1
            cohort_A = Cohort(cohort_A)
        end
        if length(cohort_B) == 1
            cohort_B = Cohort(cohort_B)
        end

        # LOG
        sum_cohort = cohort_A + cohort_B
        sm = summary(sum_cohort)
        LOG("Gen $i -> Mean of BVs: $(sm["mu_g"]), Variance of BVs: $(sm["var_g"])")
    end

    # Return
    if ratio_malefemale != 0
        return cohort_A, cohort_B
    else
        return cohort_A
    end
end

breed(cohort::Cohort; args...) = breed(cohort, cohort; args...)
breed(animal::Animal; args...) = breed(Cohort(animal); args...)
