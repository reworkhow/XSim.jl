"""
Find a 2-D coordinate from a sequential number
________________
|___|_0__1__2__|
| 0 | 0        |
| 1 | 1, 2     |
| 2 | 3, 4, 5  |
|--------------|

In |  Out
---|-------
0  | (0, 0)
1  | (1, 1)
2  | (1, 1)
3  | (2, 0)
________________________________
IF is_mating=true
________________
|___|_1__2__3__|
| 2 | 0        |
| 3 | 1, 2     |
| 4 | 3, 4, 5  |
|--------------|

In |  Out
---|-------
0  | (2, 1)
1  | (3, 1)
2  | (3, 2)
3  | (4, 1)
"""
function map_to_2D(x; i=1, is_mating=true)
    if x < i
        if is_mating
            return i + 1, x + 1
        else
            return i - 1, x
        end
    end
    return map_to_2D(x - i, i=(i + 1), is_mating=is_mating)
end


"""
Sample n_pair of mating parents from the input cohort
Return a 2D matrix in a dimenstion of n_pair by 2
"""
function sample_matings(cohort::Cohort, n_pair::Int)
    ids            = get_IDs(cohort)
    n_combinations = Int((cohort.n * (cohort.n - 1)) / 2)
    parents_1D     = sample(1:n_combinations, n_pair, replace=false) .- 1; # 0-index
    parents        = [map_to_2D(x) for x in parents_1D]
    p1             = [ids[p[1]] for p in parents]
    p2             = [ids[p[2]] for p in parents]
    return hcat(p1, p2)
end


"""
ids_bi is a 2D matrix in a dimension of mating events by two parental IDs
"""
function produce_F1_by_IDs(ids_bi)
    F1 = Cohort()
    for ids in eachrow(ids_bi)
        p = GET_LINES(Array(ids));
        if p.n == 1
            # when two parents are the same
            p1, p2 = p, p
        else
            p1, p2 = p
        end
        F1 += p1 * p2
    end
    return F1
end

n_bi = 100
n_rows = 20
n_seeds = 100
n_reps_AYT = 4

function burn_in_cycle(cohort_in ::Cohort;
                       n_bi      ::Int=100,
                       n_rows    ::Int=20,
                       n_seeds   ::Int=100,
                       n_reps_AYT::Int=4)

    # YEAR 1: Sample 100 mating parents for each biparental population
    ids_biparental = sample_matings(cohort_in, n_bi);

    # YEAR 1-2: Generate 100 DH F1
    F1  = produce_F1_by_IDs(ids_biparental)

    # YEAR 1-2: Derive 2000 (100 * 20 rows) DH lines
    DHs = get_DH(F1, n_rows);

    # YEAR 3: Headrows
    phenotypes = get_phenotypes(DHs, "XSim_ID";
                                sort=["y1"], h2=.1,
                                n_reps=n_seeds);
    ids_PYT    = phenotypes[1:500, :ID];
    lines_PYT  = GET_LINES(ids_PYT)

    # YEAR 4: PYT
    phenotypes    = get_phenotypes(lines_PYT, "XSim_ID";
                                   sort=["y1"], h2=.2,
                                   n_reps=n_seeds);
    ids_AYT       = phenotypes[1:50,  :ID];
    ids_gs_AYT    = phenotypes[1:20,  :ID];
    ids_gs_nonPYT = phenotypes[21:40, :ID];
    lines_AYT     = GET_LINES(ids_AYT)

    # YEAR 5: AYT
    phenotypes  = get_phenotypes(lines_AYT, "XSim_ID";
                                 sort=["y1"], h2=.5,
                                 n_reps=n_seeds * n_reps_AYT);
    ids_EYT     = phenotypes[1:10, :ID];
    ids_gs_EYT  = phenotypes[1:10, :ID];
    lines_EYT   = GET_LINES(ids_EYT);

    # next mating block
    cohort_out  = GET_LINES(ids_gs_nonPYT) +
                  GET_LINES(ids_gs_EYT) +
                  GET_LINES(ids_gs_AYT)
    return cohort_out
end

function burn_in(founders::Cohort;
                 n_cycles::Int=4)

    cohort = founders
    for i in 1:n_cycles
        cohort = burn_in_cycle(cohort)
    end
    return cohort
end


# # YEAR 6-7: EYT
# eval_EYT6 = get_phenotypes(adv_yr5, "XSim_ID"; h2=.67,
#                            n_reps=PARAM["n_seed"] * PARAM["n_reps_EYT"]);
# eval_EYT7 = get_phenotypes(adv_yr5, "XSim_ID"; h2=.67,
#                            n_reps=PARAM["n_seed"] * PARAM["n_reps_EYT"]);
# # YEAR 8: Variety
# data_var = hcat(get_IDs(adv_yr5), eval_EYT6[!, :y1] + eval_EYT7[!, :y1]) |> XSim.DataFrame
# sort!(data_var, "x2", rev=true)

# # parameters
# PARAM = Dict(
#     "n_founders" => 50,
#     "n_matings"  => 100,
#     "n_bi"       => 100,
#     "n_rows"     => 20,
#     "n_seed"     => 100,
#     "n_reps_AYT" => 4,
#     "n_reps_EYT" => 8
# )
