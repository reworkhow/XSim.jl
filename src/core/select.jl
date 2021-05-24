function select(cohort      ::Cohort,
                n           ::Int64;
                h2          ::Union{Array{Float64}, Float64}=.5,
                Ve          ::Union{Array{Float64}, Float64}=get_Ve(GLOBAL("n_traits"), GLOBAL("Vg"), h2),
                weights     ::Array{Float64, 1}             =[1.0],
                is_positive ::Bool                          =true,
                is_random   ::Bool                          =false,
                silent      ::Bool                          =GLOBAL("silent"),
                args...)

    phenotypes, Ve   = get_phenotypes(cohort, h2=h2, Ve=Ve, return_Ve=true)

    # Skip selection
    if cohort.n == n
        idx_sel = 1:cohort.n

    # Random select
    elseif is_random
        idx_sel = sample(1:cohort.n, n, replace=false)

    # Select by phenotypes
    else
        n_traits     = GLOBAL("n_traits")
        weights      = length(weights) != n_traits ? ones(n_traits) : weights
        direction    = is_positive * 2 - 1 # turn bool to 1 or -1
        select_index = phenotypes * weights * direction
        idx_sel      = (1:cohort.n)[sortperm(select_index, rev=true)][1:n]
    end
    cohort_sel = cohort[idx_sel]

    log_select(silent, cohort, idx_sel, phenotypes, Ve, n)

    return cohort_sel

end

select(cohort::Cohort, ratio::Float64; args...) =
    select(cohort, convert(Int64, cohort.n * ratio))

function log_select(silent, cohort, idx_sel, phenotypes, Ve, n)
    if !silent
        # Compute summary info
        g_ori = get_BVs(cohort)
        g_sel = g_ori[idx_sel, :]
        p_ori = phenotypes
        p_sel = p_ori[idx_sel, :]

        # Compute stats
        p_ori_mu  = round.(XSim.mean(p_ori, dims=1), digits=3)
        p_sel_mu  = round.(XSim.mean(p_sel, dims=1), digits=3)
        p_ori_sd  = round.(XSim.std(p_ori, dims=1), digits=3)

        g_ori_mu = round.(XSim.mean(g_ori, dims=1), digits=3)
        g_sel_mu = round.(XSim.mean(g_sel, dims=1), digits=3)
        g_ori_sd = round.(XSim.std(g_ori, dims=1), digits=3)
        # g_sel_sd = round.(XSim.std(g_sel, dims=1), digits=3)

        # Compute results
        selection_differential = round.((p_sel_mu .- p_ori_mu) ./ p_ori_sd, digits=3)
        selection_response     = round.((g_sel_mu .- g_ori_mu) ./ g_ori_sd, digits=3)

        # Print
        LOG("--------- Selection Summary ---------")
        LOG("Select $n individuals out of $(cohort.n) individuals")
        LOG("Selection differential (P): $selection_differential")
        LOG("Selection response     (G): $selection_response")
        @info "" Residual_Variance=Ve
        LOG("--------- Offsprings Summary ---------")
        # print(cohort)
    end
end