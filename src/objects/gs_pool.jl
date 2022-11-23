mutable struct GS_pool

    cohort     ::Cohort
    phenotypes ::DataFrame
    effects    ::Array
    n          ::Int # number of individuals
    p          ::Int # number of traits

    # Base constructor
    function GS_pool()
        cohort      = Cohort()
        phenotypes  = DataFrame()
        effects     = Array([])
        instance    = new(cohort, phenotypes, effects, 0, 0)
        return instance
    end
end

function isempty(gs_pool::GS_pool)
    return gs_pool.n == 0
end

function update!(gs_pool::GS_pool,
    new_cohort::Cohort;
    h2::Float64=0.5)

    # update cohort and phenotypes
    add_cohort!(gs_pool, new_cohort)
    # update DH phenotypes
    dh_p = get_phenotypes(new_cohort, "JWAS", h2=h2)
    add_phenotypes!(gs_pool, dh_p)
    # update GS effects
    update_effects!(gs_pool)
end

function add_cohort!(pool::GS_pool, new_cohort::Cohort)
    pool.cohort += new_cohort
    pool.n = pool.cohort.n
end

function add_phenotypes!(pool::GS_pool, new_phenotypes::DataFrame)
    if nrow(pool.phenotypes) != 0
        pool.phenotypes = vcat(pool.phenotypes, new_phenotypes)
    else
        pool.phenotypes = new_phenotypes
    end
    pool.p = ncol(pool.phenotypes) - 1 # the 1st col is ID
end

function update_effects!(pool::GS_pool)
    gs = genetic_evaluation(pool.cohort, pool.phenotypes;
                            methods="BayesC", return_out=true)
    # get the dimension of marker effects
    effects  = gs["marker effects geno"][!, "Estimate"]
    n_traits = pool.p
    n_loci   = Int(length(effects) / n_traits)
    # reshape to n by p
    pool.effects = reshape(effects, (n_loci, n_traits))
end
