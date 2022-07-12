
"""
________________
|___|_1__2__3__|
| 1 | 0        |
| 2 | 1, 2     |
| 3 | 3, 4, 5  |
|--------------|

In |  Out
---|-------
0  | (1, 1)
1  | (2, 1)
2  | (2, 2)
3  | (3, 1)
"""
function map_to_half_2D(x; i=1)
    if x < i
        return i, x + 1
    end
    return map_to_half_2D(x - i, i=i + 1)
end

"""
n = 3 (n col)
________________
|___|_1__2__3__|
| 1 | 0  1  2  |
| 2 | 3  4  5  |
| 3 | 6  7  8  |
| 4 | 9 10 11  |
|--------------|

In |  Out
---|-------
0  | (1, 1)
1  | (1, 2)
2  | (1, 3)
3  | (2, 1)
"""
function map_to_full_2D(x, n; i=1)
    if x < n
        return i, x + 1 # row, col
    end
    return map_to_full_2D(x - n, n; i=i + 1)
end


"""
Sample n_pair of mating parents from the input cohort(s)
Return a 2D matrix in a dimenstion of n_pair by 2 (sire and dams)
"""
function sample_matings(cohort_a::Cohort, # rows
    cohort_b::Cohort, # columns
    n_pair::Int;
    half::Bool=false)
    ids_a = get_IDs(cohort_a)
    ids_b = get_IDs(cohort_b)
    if half
        n_combinations = Int((cohort_a.n * (cohort_b.n - 1)) / 2) + cohort_a.n
    else
        n_combinations = cohort_a.n * cohort_b.n
    end
    parents_1D = sample(1:n_combinations, n_pair, replace=false) .- 1 # 0-index

    if half
        parents = [map_to_half_2D(x) for x in parents_1D]
    else
        parents = [map_to_full_2D(x, cohort_b.n) for x in parents_1D]
    end
    p1 = [ids_a[p[1]] for p in parents]
    p2 = [ids_b[p[2]] for p in parents]
    return hcat(p1, p2)
end


"""
### Test codes
n_bi = 100
n_rows = 20
n_seeds = 100
n_reps_AYT = 4
"""
function burn_in_cycle(cohort_in::Cohort;
    n_bi::Int=100,
    n_rows::Int=20,
    n_seeds::Int=100,
    n_reps_AYT::Int=4)

    # YEAR 1: Sample 100 mating parents for each biparental population
    ids_biparental = sample_matings(cohort_in, n_bi)

    # YEAR 1-2: Generate 100 DH F1
    F1 = produce_F1_by_IDs(ids_biparental)

    # YEAR 1-2: Derive 2000 (100 * 20 rows) DH lines
    DHs = get_DH(F1, n_rows)

    # YEAR 3: Headrows
    phenotypes = get_phenotypes(DHs, "XSim_ID";
        sort=["y1"], h2=0.1,
        n_reps=n_seeds)
    ids_PYT = phenotypes[1:500, :ID]
    lines_PYT = GET_LINES(ids_PYT)

    # YEAR 4: PYT
    phenotypes = get_phenotypes(lines_PYT, "XSim_ID";
        sort=["y1"], h2=0.2,
        n_reps=n_seeds)
    ids_AYT = phenotypes[1:50, :ID]
    ids_gs_AYT = phenotypes[1:20, :ID]
    ids_gs_nonPYT = phenotypes[21:40, :ID]
    lines_AYT = GET_LINES(ids_AYT)

    # YEAR 5: AYT
    phenotypes = get_phenotypes(lines_AYT, "XSim_ID";
        sort=["y1"], h2=0.5,
        n_reps=n_seeds * n_reps_AYT)
    ids_EYT = phenotypes[1:10, :ID]
    ids_gs_EYT = phenotypes[1:10, :ID]
    lines_EYT = GET_LINES(ids_EYT)

    # next mating block
    cohort_out = GET_LINES(ids_gs_nonPYT) +
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

