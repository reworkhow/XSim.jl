# Founders
```julia; echo=false
using CSV
using DataFrames
```

### Step 0. Load XSim and set random seed
```julia; eval=true;
using XSim
```

### Step 1. Setup genome and phenome with small examples
```julia; eval=false
XSim.build_demo_small()
```
```
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info: 
[ Info: Chromosome Length (cM): 16.0
[ Info: [8.0, 8.0]
[ Info: 
[ Info: Number of Loci        : 10
[ Info: [5, 5]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 2
┌ Info: 
│   Genetic_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
[ Info: Number of QTLs        : [2 5]
```

### Step 2. Inspect genotypes
```julia;results="hidden";echo=false
root      = dirname(dirname(pathof(XSim)))
filepath  = joinpath(root, "data", "demo_genotypes.csv")
genotypes = CSV.read(filepath, DataFrame, header=false)

```

```julia; eval=false
# XSim provide example data with XSim.data() function
genotypes = XSim.data("genotypes")
```
```julia; echo=false
genotypes
```

### Step 3. Initialize founders by files or dataframes
```julia; eval=false
cohort = Founders(genotypes)
```
```
[ Info: MAF has been updated based on provided haplotypes/genotypes
[ Info: Cohort (5 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.003 0.965]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [1.734 1.261]
```

```julia; eval=false
root      = dirname(dirname(pathof(XSim)))
filepath  = joinpath(root, "data", "demo_genotypes.csv")
cohort    = Founders(filepath)
```
```
[ Info: MAF has been updated based on provided haplotypes/genotypes
[ Info: Cohort (5 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.003 0.965]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [1.734 1.261]
```

Inspect the founders
```julia; eval=false
get_genotypes(cohort)
```
```
5×10 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 2  0  0  0  1  0  0  2  1  0
 0  0  0  1  0  0  0  2  0  1
 0  1  0  0  2  1  0  0  0  0
 1  1  0  0  2  1  0  0  2  0
 2  2  0  1  2  0  0  2  0  2
```

```julia; eval=false
get_QTLs(cohort)
```
```
5×6 Array{Int64,2}:
2  0  0  0  0  2
0  0  0  0  0  2
0  1  0  1  0  0
1  1  0  1  0  0
2  2  0  0  0  2
```

```julia; eval=false
get_BVs(cohort)
```
```
5×2 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:
0.0         3.07603
0.0         1.47043
0.00488477  0.780246
0.00488477  1.58304
0.0         4.63652
```

```julia; eval=false
get_pedigrees(cohort)
```
```
5×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
6  0  0
7  0  0
8  0  0
9  0  0
10  0  0
```

### Step 4. Mate and select for 5 generations
```julia; eval=false
args = Dict(:n_per_mate      => 4,
            :n_gens          => 5,
            :ratio_malefemale=> 2,
            :h2              => [.8, .8])
males, females = breed(cohort, cohort; args...)
```
```

[ Info: Gen 0 -> Mean of BVs: [0.002 2.309], Variance of BVs: [0.0 2.127]
[ Info: Gen 1 -> Mean of BVs: [0.0 3.48], Variance of BVs: [0.0 0.581]
[ Info: Gen 2 -> Mean of BVs: [0.0 4.026], Variance of BVs: [0.0 0.239]
[ Info: Gen 3 -> Mean of BVs: [0.0 4.177], Variance of BVs: [0.0 0.411]
[ Info: Gen 4 -> Mean of BVs: [0.0 4.402], Variance of BVs: [0.0 0.142]
[ Info: Gen 5 -> Mean of BVs: [0.0 4.48], Variance of BVs: [0.0 0.108]
```

### Step 5. Examine the results
```julia;eval=false
summary(males + females)
```
```
Dict{String,Any} with 3 entries:
  "mu_g"  => [0.0 4.48]
  "var_g" => [0.0 0.108]
  "n"     => 10
```

Compare the allele frequencies
```julia;eval=false
get_QTLs(cohort) |> get_MAF
```
```
6-element Array{Float64,1}:
0.5
0.4
0.0
0.2
0.0
0.4
```

```julia;eval=false
get_QTLs(males + females) |> get_MAF
```
```
6-element Array{Float64,1}:
0.0
0.1
0.0
0.0
0.0
0.0
```