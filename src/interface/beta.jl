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