function select(cohort            ::Cohort,
                n                 ::Int;
                h2                ::Union{Array{Float64}, Float64}=.5,
                Ve                ::Union{Array{Float64}, Float64}=-999.99,
                weights           ::Array{Float64, 1}=[1.0],
                is_positive_select::Bool=true)

    n_traits     = GLOBAL("n_traits")
    phenotypes   = get_phenotypes(cohort, h2=h2, Ve=Ve)
    weights      = length(weights) != n_traits ? ones(n_traits) : weights
    direction    = is_positive_select * 2 - 1 # turn bool to 1 or -1
    select_index = phenotypes * weights * direction
    return cohort[sortperm(select_index, rev=true)][1:n]
end

