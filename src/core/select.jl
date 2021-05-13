function select(cohort            ::Cohort,
                n                 ::Int;
                Ve                ::Union{Array{Float64}, Float64},
                weights           ::Array{Float64, 1}=[1.0],
                is_positive_select::Bool=true)

    traits   = get_traits(cohort, criteria)
    n_traits = size(traits, 2)
    weights  = length(weights) != n_traits ? ones(n_traits) : weights
    direction = is_positive_select * 2 - 1 # turn bool to 1 or -1
    select_index = traits * weights * direction
    return cohort[sortperm(select_index, rev=true)][1:n]
end



