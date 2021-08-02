"""
# Genetic evaluation
    genetic_evaluation(cohort         ::Cohort,
                       phenotypes     ::DataFrame=DataFrame();
                       model_equation ::String="",
                       covariates     ::String="",
                       random_iid     ::String="",
                       random_str     ::String="",
                       methods        ::String="GBLUP",
                       add_genotypes  ::Bool=true,
                       idx_missing_p  ::Any=[],
                       return_out     ::Bool=false,
                       args...)
### Arguments
- `cohort` : The evaluated `cohort`.
- `phenotypes` : A `dataframe` wiht columns of `id`, traits, and other studied factors.
- `h2` : Default 0.5. Heritability for simulating cohort phenotypes if no `phenotype` is provided.
- `ve` : Default 1. Residual variance for simulating cohort phenotypes if no `phenotype` is provided.
- `model_equation` : Default "y ~ intercept". The equation for fitting the breeding values.
- `covariates` : Specifies terms in `phenotypes` as continuous factors. Must be included in `model_equation` as well.
- `random_iid` : Specifies terms in `phenotypes` as random effects (i.i.d.). Must be included in `model_equation` as well.
- `random_str` : Specifies terms in `phenotypes` as random effects (pedigree). Must be included in `model_equation` as well.
- `methods` : Default "GBLUP". Defines what models/methods used in the genetic evaluation.
- `add_genotypes` : Default `true`. Genotypes will be included in the equation.
- `idx_missing_p` : A vector assigning which individuals are not phenotyped.
- `return_out` : Default `false`. Set to `true` to return the complete `JWAS` outputs.

### Outputs
A `n` by `t` matrix containing breeding values will be return if `return_out = false`, where `n` is number of individuals and `t` is number of evaluated traits. A complete `JWAS` outputs will be returned if `return_out = true`.

### Example of the `phenotypes` dataframe
```jldoctest
 Row │ ID           y1              y2               factor_1  factor_2
     │ String       Float64?        Float64?         Int64     Int64
─────┼─────────────────────────────────────────────────────────────
   1 │ 1             0.88976        -0.0798048         1         1
   2 │ 2            -0.783203        0.988616          1         1
   3 │ 3             missing         missing           1         1
   4 │ 4             missing         missing           1         1
   5 │ 5             missing         missing           1         1
   6 │ 6             missing         missing           1         2
   7 │ 7             missing         missing           1         2
   8 │ 8            -1.76058         0.277289          1         2
   9 │ 9             0.938871        2.57784           1         2
  10 │ 10            0.37026         3.15993           1         2
  11 │ 11           -1.91869         0.0935064         2         3
  12 │ 12           -0.89847         1.8987            2         3
  13 │ 13            1.69663        -0.949513          2         3
  14 │ 14            3.48862         0.654378          2         3
  15 │ 15            1.39615         1.80355           2         3
  16 │ 16            2.31685         2.13446           2         4
  17 │ 17           -3.81017         0.0186156         2         4
  18 │ 18           -1.71216        -0.0976809         2         4
  19 │ 19            1.80917         1.34104           2         4
  20 │ 20           -0.504771        3.22665           2         4
```

### Examples
We will use demo `genome` and `phenome` for the example. A `cohort` with 20 individuals is simulated.

```jldoctest
julia> build_demo()
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 10
[ Info: 
[ Info: Chromosome Length (cM):
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

julia> cohort = Cohort(20)
[ Info: Cohort (20 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.862 -0.913]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.814 1.448]
```

─────────────────────────────────────────────────────────
#### Basic usage
By default, GBLUP will be performed without providing any argument.
```jldoctest
julia> out = genetic_evaluation(cohort)
20×2 Array{Float64,2}:
  42.2908     -5.89224
   6.38041    -4.5895
  19.5065    -57.4513
  24.05      -84.5086
  11.5925     45.5247
 -11.5926     -1.91879
 -29.2152    -35.5977
  -1.92563    63.8841
  19.8686     76.2722
  22.6728     45.5179
 -58.2754     42.5953
  28.556     -21.0314
  11.5918     25.1126
  -0.517653   -0.237724
 -11.4712    -47.3073
 -10.9575     -3.58754
 -11.7835     -9.62802
 -17.7895    -22.0821
 -39.5114     -3.46026
   6.53008    -1.61416
```

The breeding values `out` can be used as criteria for the selection directly.
```jldoctest
julia> select(cohort, 5, out)
[ Info: --------- Selection Summary ---------
[ Info: Select 5 individuals out of 20 individuals
[ Info: Selection differential (P): [0.783 1.192]
[ Info: Selection response     (G): [0.344 1.063]
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (5 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.552 0.367]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.181 0.539]
```

Conditions of the evaluation can be further defined by `h2` or `ve`.
```jldoctest
julia> out = genetic_evaluation(cohort, ve = [1  0.5
                                             0.5  1])
```

Set `return_out = true` to obtain the complete `JWAS` outputs.
```jldoctest
julia> out = genetic_evaluation(cohort, return_out = true)
Dict{Any,Any} with 7 entries:
  "EBV_y2"              => 20×3 DataFrame…
  "EBV_y1"              => 20×3 DataFrame…
  "heritability"        => 2×3 DataFrame…
  "location parameters" => 2×5 DataFrame…
  "residual variance"   => 4×3 DataFrame…
  "marker effects geno" => 40×5 DataFrame…
  "genetic_variance"    => 4×3 DataFrame…
```

─────────────────────────────────────────────────────────
#### Customized phenotypes and factors.
Obtain JWAS-compatible dataframe. 
```jldoctest
dt_p = get_phenotypes(cohort, "JWAS")
```

Assign un-phenotyped individuals
```jldoctest
julia> idx = 3:6
julia> allowmissing!(dt_p);
julia> dt_p[idx, 2:end] .= missing;
julia> first(dt_p, 10)
10×3 DataFrame
 Row │ ID            y1               y2             
     │ String?       Float64?         Float64?       
─────┼──────────────────────────────────────────
   1 │ 1             -0.0933375       -0.882781
   2 │ 2             -1.50748         -2.12898
   3 │ 3              missing          missing        
   4 │ 4              missing          missing        
   5 │ 5              missing          missing        
   6 │ 6              missing          missing        
   7 │ 7             -0.431361        -0.624666
   8 │ 8             -1.17867          0.415607
   9 │ 9              0.266733         1.02123
  10 │ 10            -1.15588          0.737982
```

Provide customized phenotypes for the evaluation
```jldoctest
julia> out = genetic_evaluation(dams, dt_p)
```

It's equivalent to set `idx_missing_p = 3:6` for the missing phenotypes.
```jldoctest
julia> out = genetic_evaluation(dams, dt_p, idx_missing_p = 3:6)
The marker IDs are set to 1,2,...,#markers
The individual IDs is set to 1,2,...,#observations
Genotype informatin:
#markers: 1000; #individuals: 20
The folder results is created to save results.
Checking genotypes...
Checking phenotypes...
Individual IDs (strings) are provided in the first column of the phenotypic data.
Phenotypes for 16 individuals are used in the analysis.These individual IDs are saved in the file IDs_for_individuals_with_phenotypes.txt.
```

We can simulated `factor_1` and `factor_2` as fixed and random effects, respectively. And we use both `factor_1` and `factor_2` to fit `y1`,  and `factor_1` alone to fit `y2`.

```jldoctest
julia> dt_p[:, "factor_1"] = [i for i in 1:4 for j in 1:5];
julia> dt_p[:, "factor_2"] = [i for i in 1:2 for j in 1:10];
julia> out = genetic_evaluation(cohort, dt_p,
                model_equation="y1 = intercept + factor_1 + factor_2
                                y2 = intercept + factor_1",
                random_iid="factor_2",
                return_out=true)

factor_2 is not found in model equation 2.
The marker IDs are set to 1,2,...,#markers
The individual IDs is set to 1,2,...,#observations
Genotype informatin:
#markers: 1000; #individuals: 20
The folder results is created to save results.
Checking genotypes...
Checking phenotypes...
Individual IDs (strings) are provided in the first column of the phenotypic data.
Phenotypes for 16 individuals are used in the analysis.These individual IDs are saved in the file IDs_for_individuals_with_phenotypes.txt.
Prior information for random effect variance is not provided and is generated from the data.

Pi (Π) is not provided.
Pi (Π) is generated assuming all markers have effects on all traits.

A Linear Mixed Model was build using model equations:

y1 = intercept + factor_1 + factor_2
y2 = intercept + factor_1

Model Information:

Term            C/F          F/R            nLevels
intercept       factor       fixed                1
factor_1        factor       fixed                4
factor_2        factor       random               2

MCMC Information:

chain_length                                    100
burnin                                            0
starting_value                                 true
printout_frequency                              101
output_samples_frequency                          1
constraint                                    false
missing_phenotypes                             true
update_priors_frequency                           0
seed                                          false

Hyper-parameters Information:

random effect variances (y1:factor_2):
 0.45
residual variances:           
 1.0f0  0.0f0
 0.0f0  1.0f0

Genomic Information:

complete genomic data (i.e., non-single-step analysis)

Genomic Category                               geno
Method                                        GBLUP
genetic variances (genomic):  
 1.0  0.0
 0.0  1.0
estimateScale                                 false

Degree of freedom for hyper-parameters:

residual variances:                           6.000
random effect variances:                      5.000
marker effect variances:                      6.000

The version of Julia and Platform in use:

Julia Version 1.5.4
Commit 69fcb5745b (2021-03-11 19:13 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin18.7.0)
  CPU: Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-9.0.1 (ORCJIT, skylake)
Environment:
  JULIA_EDITOR = code
  JULIA_NUM_THREADS = 

The analysis has finished. Results are saved in the returned variable and text files. MCMC samples are saved in text files.

Dict{Any,Any} with 7 entries:
  "EBV_y2"              => 20×3 DataFrame…
  "EBV_y1"              => 20×3 DataFrame…
  "heritability"        => 2×3 DataFrame…
  "location parameters" => 12×5 DataFrame…
  "residual variance"   => 4×3 DataFrame…
  "marker effects geno" => 32×5 DataFrame…
  "genetic_variance"    => 4×3 DataFrame…
```

"""
function genetic_evaluation(cohort         ::Cohort,
                            phenotypes     ::DataFrame=DataFrame();
                            ve             ::Union{Array{Float64}, Float64}=GLOBAL("Ve"),
                            model_equation ::String="",
                            covariates     ::String="",
                            random_iid     ::String="",
                            random_str     ::String="",
                            methods        ::String="GBLUP",
                            add_genotypes  ::Bool=true,
                            idx_missing_p  ::Any=[],
                            return_out     ::Bool=false,
                            args...)

    # 1. Acquire phenotypes
    if nrow(phenotypes) == 0
        phenotypes = get_phenotypes(cohort, "JWAS"; args...)
    end

    # 1.1 If missing is provided
    if length(idx_missing_p) != 0
        XSim.allowmissing!(phenotypes);
        phenotypes[idx_missing_p, 2:end] .= missing
    end

    # 2. Build Model Equations
    if model_equation == ""
        traits = names(phenotypes)[2:end]
        array_eq = ["$(trait) = intercept" for trait in traits]
        model_equation = join(array_eq, "\n")
    end

    # 3. Build model
    model = JWAS.build_model(model_equation, ve);

    # 4: Set Factors or Covariates
    if covariates != ""
        JWAS.set_covariate(model, covariates);
    end

    # 5: Set Random or Fixed Effects
    if random_iid != ""
        JWAS.set_random(model, random_iid);
    end
    if random_str != ""
        pedigree = get_pedigree(cohort, "JWAS")
        JWAS.set_random(model, random_str, pedigree);
    end

    # 6. If GBLUP, add genotypes
    if add_genotypes
        # Add genotype for GBLUP
        genotypes= get_genotypes(cohort) # 0 1 2
        JWAS.add_genotypes(model, float.(genotypes), GLOBAL("Vg"))
    end


    # 7. Run MCMC
    out = JWAS.runMCMC(model, phenotypes, methods=methods);

    # Remove outputs
    try
        rm("results", recursive=true)
        rm("IDs_for_individuals_with_genotypes.txt")
        rm("IDs_for_individuals_with_pedigree.txt")
        rm("IDs_for_individuals_with_phenotypes.txt")
    catch e
        nothing
    end

    # Outputs
    if return_out
        # return complete JWAS report
        return out
    else
        # return EBV only
        return hcat([out["EBV_$x"][:, "EBV"] for x in traits]...)
    end

    # Note
    # jwas_ped      = get_pedigree(cohort, "JWAS");
    # jwas_p        = get_phenotypes(cohort, "JWAS");
    # allowmissing!(jwas_p);
    # true_p = jwas_p.y1[1:200]
    # jwas_p.y1[1:200] .= missing;

    # genotypes     = get_genotypes(cohort, "JWAS");
    # model         = XSim.JWAS.build_model(model_equation);
    # out = XSim.JWAS.runMCMC(model, jwas_p, methods="GBLUP")

    # Not working
    # eq="y1 = intercept + genotypes"
    # genetic_evaluation(cohort, model_equation=eq)
end
