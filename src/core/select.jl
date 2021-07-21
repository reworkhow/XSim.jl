function select(cohort      ::Cohort,
                n           ::Int64,
                criteria    ::String = "phenotypes";
                h2          ::Union{Array{Float64}, Float64}=GLOBAL("h2"),
                ve          ::Union{Array{Float64}, Float64}=GLOBAL("Ve"),
                weights     ::Array{Float64, 1}             =[1.0],
                return_log  ::Bool                          =false,
                is_random   ::Bool                          =false,
                silent      ::Bool                          =GLOBAL("silent"),
                args...)

    # Computation ----------------------------------------------------------
    # Phenotype
    if criteria == "phenotypes"
        phenotypes, ve = get_phenotypes(cohort, "XSim", h2=h2, ve=ve, return_ve=true)
        values_select  = phenotypes
    elseif criteria == "GBLUP"
        phenotypes, ve = get_phenotypes(cohort, "JWAS", h2=h2, ve=ve, return_ve=true)
        values_select  = GBLUP(cohort, phenotypes)
        phenotypes     = Matrix(phenotypes[:, 2:end]) # turn JWAS objects to regular dataframe
    else
        LOG("Available criteria are: ['phenotypes', 'GBLUP']", "error")
    end

    # Selection ------------------------------------------------------------
    # Skip selection
    if cohort.n <= n
        idx_sel = 1:cohort.n

    # Random select
    elseif is_random
        idx_sel = sample(1:cohort.n, n, replace=false)

    # Select by phenotypes
    else
        n_traits     = GLOBAL("n_traits")
        weights      = length(weights) != n_traits ? ones(n_traits) : weights
        select_index = values_select * weights
        idx_sel      = (1:cohort.n)[sortperm(select_index, rev=true)][1:n]
    end
    cohort_sel = cohort[idx_sel]

    # Log ------------------------------------------------------------
    if return_log
        sel_P, sel_G = log_select(silent, cohort, idx_sel, phenotypes, ve, n, true)
        return Dict("cohort"=>cohort_sel, "sel_P"=>sel_P, "sel_G"=>sel_G)
    else
        log_select(silent, cohort, idx_sel, phenotypes, ve, n, false)
        return cohort_sel
    end

end

function select(cohort  ::Cohort,
                subset_n::Array{Int64};
                args...)

    cohort_sel = select(cohort, maximum(subset_n); args...)
    return cohort_sel[subset_n]
end

select(cohort::Cohort, ratio::Float64; args...) =
    select(cohort, convert(Int64, cohort.n * ratio))

function log_select(silent, cohort, idx_sel, phenotypes, Ve, n, return_log::Bool=false)
    if !silent
        # Compute summary info
        g_ori = get_BVs(cohort)
        g_sel = g_ori[idx_sel, :]
        p_ori = phenotypes
        p_sel = p_ori[idx_sel, :]

        # Compute stats
        p_ori_mu = round.(XSim.mean(p_ori, dims=1), digits=3)
        p_sel_mu = round.(XSim.mean(p_sel, dims=1), digits=3)
        p_ori_sd = round.(XSim.std(p_ori, dims=1), digits=3)

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

        # Log
        if return_log
            return selection_differential, selection_response
        end
    end
end