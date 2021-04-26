function select(cohort::Cohort, n::Int;
                weights::Array{Float64, 1}=[1.0],
                is_positive_select::Bool=true,
                criteria="phenotypic")

    traits   = get_traits(cohort, criteria)
    n_traits = size(traits, 2)
    weights  = length(weights) != n_traits ? ones(n_traits) : weights
    direction = is_positive_select * 2 - 1 # turn bool to 1 or -1
    select_index = traits * weights * direction
    return cohort[sortperm(select_index, rev=true)][1:n]
end