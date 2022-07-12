"""
# Selection function
    select(cohort      ::Cohort,
           n           ::Int64,
           criteria    ::Union{String, Array} = "phenotypes";
           h2          ::Union{Array{Float64}, Float64}=GLOBAL("h2"),
           ve          ::Union{Array{Float64}, Float64}=GLOBAL("Ve"),
           weights     ::Array{Float64, 1}  =[1.0],
           return_log  ::Bool               =false,
           is_random   ::Bool               =false,
           silent      ::Bool               =GLOBAL("silent")

    select(cohort::Cohort, ratio::Float64; args...)

### Arguments
Positional arguments
- `cohort` : A `cohort` from which individuals are selected.
- `n` : `n` individuals are selected.
- `ratio` : `ratio` portion of individuals are selected.
- `criteria` : `Criteria` that will be used for the selecition. Default
  "phenotypes", the options are ["phenotypes", "GBLUP", array]. If set to
  "GBLUP",  a genetic evaluation is carried out by `JWAS` and the estimated
  breeding values will be the `criteria`. It's also avaialbe to provdie
  the `criteria` (e.g., phenotypes matrix) directly for the selection.

Keyword arguments
- `h2` : The heritability `h2` of the simulated phenotypes.
- `ve` : The residual covariance `ve` of the simulated phenotypes.
- `weight` : Linear coefficients of traits for the selection. The selection is
  more sensitive to traits with greater `weight`. Negative
- `return_log` : Default `false`. Set `true` to return selection differential
  and selection response besides the selected cohort.
- `silent` : Default `false`. Set `true` to mute the log messages.

### Outputs
A selected `cohort` object will be returned. If `return_log` is set to `true`,
a `dictionary` object containing the selected cohort, selection differential,
and selection response will be returned.

─────────────────────────────────────────────────────────
### Example
#### Single Trait Selection
Set demo genome and phenome with single traits controlled by 50 QTLs.
```jldoctest
julia> build_demo()
julia> build_phenome(50)
julia> summary()
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 10
[ Info: 
[ Info: Chromosome Length (cM): 1500.0
[ Info: [150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0]
[ Info: 
[ Info: Number of Loci        : 1000
[ Info: [100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 1
[ Info: Heritability (h2)     : [0.5]
┌ Info: 
│   Genetic_Variance =
│    1×1 Array{Float64,2}:
└     1.0
┌ Info: 
│   Residual_Variance =
│    1×1 Array{Float64,2}:
└     1.0
[ Info: Number of QTLs        : [50]
```

Initialize a cohort with 100 individuals
```jldoctest
julia> cohort = Cohort(100)
```

##### Select 30 individuals
```jldoctest
# Select top 30 individuals
julia> cohort_s = select(cohort, 30)
# Equivalent
julia> cohort_s = select(cohort, 0.3)

[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [1.174]
[ Info: Selection response     (G): [0.843]
┌ Info:
│   Residual_Variance =
│    1×1 Array{Float64,2}:
└     1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [1.448]
[ Info:
[ Info: Variance of breeding values:
[ Info: [0.367]
```

##### Assign Heritability `h2` or Residual Covariance `ve`
```jldoctest
# Assign heritability
julia> progenies = select(cohort, 30, h2=0.1)

# Equivalent in the case where genetic variance `vg` is 1.0
julia> progenies = select(cohort, 30, ve=9.0)

[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [1.182]
[ Info: Selection response     (G): [0.338]
┌ Info:
│   Residual_Variance =
│    1×1 Array{Float64,2}:
└     9.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [0.956]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.643]
```

##### Negative Selection
Set `is_positive=false` to rank individuals in ascending order
```jldoctest
julia> progenies = select(cohort, 30, is_positive=false)
[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [-1.19]
[ Info: Selection response     (G): [-0.89]
┌ Info: 
│   Residual_Variance =
│    1×1 Array{Float64,2}:
└     1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.24]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.566]
```

##### Random Selection
```jldoctest
julia> progenies = select(cohort, 30, is_random=true)
[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [-0.06]
[ Info: Selection response     (G): [-0.191]
┌ Info: 
│   Residual_Variance =
│    1×1 Array{Float64,2}:
└     1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [0.441]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.946]
```

##### Selection wiht Multiple Parameters
It's possible specify multiple parameters described above in one selection.
User can either enclose parameters as keyword arguments, or pass them through a `dictionary` object.

```jldoctest
# Keyword args
julia> progenies = select(cohort, 30, h2=0.3, is_positive=false)

# Equivalent
julia> args = Dict(:h2=>0.3,
                   :is_positive=>false)
julia> progenies = select(cohort, 30; args...)

[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [-1.086]
[ Info: Selection response     (G): [-0.486]
┌ Info: 
│   Residual_Variance =
│    1×1 Array{Float64,2}:
└     2.3333333333333335
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [0.154]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.818]
```

─────────────────────────────────────────────────────────
#### Multi-Trait Selection
Set demo genome and phenome with single traits controlled by 50 QTLs.
```jldoctest
julia> build_demo()
julia> summary()
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 10
[ Info: 
[ Info: Chromosome Length (cM): 1500.0
[ Info: [150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0]
[ Info: 
[ Info: Number of Loci        : 1000
[ Info: [100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 2
[ Info: Heritability (h2)     : [0.5, 0.5]
┌ Info: 
│   Genetic_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
[ Info: Number of QTLs        : [3 8]
```

Initialize a cohort with 100 individuals
```jldoctest
julia> cohort = Cohort(100)
```

##### Assign Heritabilities for Multiple Traits
```jldoctest
julia> progenies = select(cohort, 30, h2=[0.9, 0.3])
[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [0.468 1.028]
[ Info: Selection response     (G): [0.383 0.636]
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     0.111111  0.0
└     0.0       2.33333
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.889 0.28]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.947 0.625]
```

##### Assign Trait Correlations via Residual Covariance
```jldoctest
julia> progenies = select(cohort, 30, ve=[1   0.3
                                          0.3   1])
[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [0.866 0.925]
[ Info: Selection response     (G): [0.662 0.762]
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.3
└     0.3  1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.608 0.406]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.549 0.476]
```

##### Derive Selection Index for Multiple Traits
Assigning a vector to the parameter `weights` to derive a selection index which is a linear combintation of the weights and the phenotypes. 
In this example, we demonstrate two traits with the heritability of 0.3 and 0.8, respectively.
And we can select traits with more weight on the second trait which is more heritable, and negatively select the first trait.

```jldoctest
julia> progenies = select(cohort, 30, h2=[.3, .8], weights=[-0.1, 0.9])
[ Info: --------- Selection Summary ---------
[ Info: Select 30 individuals out of 100 individuals
[ Info: Selection differential (P): [-0.318 1.027]
[ Info: Selection response     (G): [-0.233 0.869]
┌ Info:
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     2.33333  0.0
└     0.0      0.25
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [-1.508 0.513]
[ Info:
[ Info: Variance of breeding values:
[ Info: [1.053 0.458]
```
"""
function select(cohort      ::Cohort,
                n           ::Int64;
                criteria    ::Union{String, Array} = "phenotypes",
                h2          ::Union{Array{Float64}, Float64}=GLOBAL("h2"),
                ve          ::Union{Array{Float64}, Float64}=GLOBAL("Ve"),
                weights     ::Array{Float64, 1}             =[1.0],
                return_log  ::Bool                          =false,
                silent      ::Bool                          =GLOBAL("silent"),
                args...)

    # Computation ----------------------------------------------------------
    # Phenotype
    if criteria == "phenotypes"
        phenotypes, ve = get_phenotypes(cohort, "XSim", h2=h2, ve=ve,           return_ve=true)
        values_select  = phenotypes

    elseif criteria == "EBV"
        phenotypes, ve = get_phenotypes(cohort, "JWAS", h2=h2, ve=ve, return_ve=true)
        values_select  = genetic_evaluation(cohort, phenotypes)
        phenotypes     = Matrix(phenotypes[:, 2:end]) # turn JWAS objects to regular dataframe

    elseif criteria == "random"
        nothing

    elseif isa(criteria, Array)
        # use provided phenotypes
        values_select = criteria
        phenotypes = criteria

    else
        LOG("Available criteria are: ['phenotypes', 'EBV', 'random']", "error")
    end

    # Selection ------------------------------------------------------------
    # Skip selection
    if cohort.n <= n
        idx_sel = 1:cohort.n

    # Random select
    elseif criteria == "random"
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
    if criteria != "random"
        if return_log
            sel_P, sel_G = log_select(silent, cohort, idx_sel, phenotypes, ve, n, true)
            return Dict("cohort"=>cohort_sel, "sel_P"=>sel_P, "sel_G"=>sel_G)
        else
            log_select(silent, cohort, idx_sel, phenotypes, ve, n, false)
            return cohort_sel
        end
    else
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
        p_ori_mu = round.(mean(p_ori, dims=1), digits=3)
        p_sel_mu = round.(mean(p_sel, dims=1), digits=3)
        # p_ori_sd = round.(std(p_ori, dims=1), digits=3)

        g_ori_mu = round.(mean(g_ori, dims=1), digits=3)
        g_sel_mu = round.(mean(g_sel, dims=1), digits=3)
        # g_ori_sd = round.(std(g_ori, dims=1), digits=3)
        # g_sel_sd = round.(std(g_sel, dims=1), digits=3)

        # Compute results
        selection_differential = round.((p_sel_mu .- p_ori_mu), digits=3)
        selection_response     = round.((g_sel_mu .- p_ori_mu), digits=3)

        # Print
        LOG("--------- Selection Summary ---------")
        
        LOG("Residual (Co)variance")
        Base.print_matrix(stdout, Ve)
        LOG("")
        LOG("--------- Offsprings Summary ---------")

        # Log
        if return_log
            return selection_differential, selection_response
        end
    end
end