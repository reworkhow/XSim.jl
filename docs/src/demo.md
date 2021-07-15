# Demo: Step by Step

## Step 0. Load XSim
```jldoctest
julia> using XSim
```

## Step 1. Setup genome and phenome
The demo example simulates `10` chromosomes with `100` loci each. And `2` independent traits are controlled by `3` and `8` QTLs, respectively.

```jldoctest
julia> build_genome(n_chr=10, n_marker=100)
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 10
[ Info: 
[ Info: Chromosome Length (cM): 1000.0
[ Info: [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
[ Info: 
[ Info: Number of Loci        : 1000
[ Info: [100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 

julia> build_phenome([3, 8])
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

## Step 2. Initialize Founders
Simulate `3` founder sires

```jldoctest
julia> n_sires = 3

julia> sires   = Founders(n_sires)
[ Info: Cohort (3 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [2.232 -0.609]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [2.374 0.873]
```
Simulate `20` founder dams
```jldoctest
julia> n_dams  = 20

julia> dams    = Founders(n_dams)
[ Info: Cohort (20 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [2.117 -0.369]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.519 0.732]
```

## Step 3. Mating
Here, the first `cohort` in the function `mate` is `cohort_A`, and the second one is `cohort_B`. We let each sire mate with `5` dams, and each mating event will produce `1` male and `1` female progenies when `:ratio_malefemale` is set to `1.0`.

```jldoctest
julia> args_mate = Dict(:nA           => 3,
                        :nB_per_A     => 5,
                        :n_per_mate   => 2,
                        :ratio_malefemale => 1.0)

julia> males, females = mate(sires, dams; args_mate...)
[ Info: --------- Mating Summary ---------
[ Info: Generated 30 individuals from 3 cohort_A individuals
[ Info: Every cohort_A individual mates with 5 individuals from cohort_B
[ Info: And each mating reproduces 2 progenies
[ Info: 
[ Info: --------- Offsprings Summary ---------
([ Info: Cohort (15 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [2.294 -0.592]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.846 0.33]
, [ Info: Cohort (15 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.614 -0.586]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [1.313 0.527]
)
```

Users can check how the mating was performed by calling `get_pedigree()` function. The first column is individual ID, and the second and the third correspond to the individual's sire and dam, respectively. We can use `+` to concatenate two `cohort`.

```jldoctest
julia> progenies = males + females

julia> get_pedigree(progenies)
30×3 Array{Int64,2}:
 24  2  22
 25  2  22
 26  2   7
 27  2   7
 28  2  11
 29  2  11
 30  2   6
 31  2   6
 32  2  17
 33  2  17
 34  1  23
 35  1  23
 36  1  19
 37  1  19
 38  1   4
 39  1   4
 40  1  16
 41  1  16
 42  1  13
 43  1  13
 44  3  20
 45  3  20
 46  3  18
 47  3  18
 48  3  12
 49  3  12
 50  3   9
 51  3   9
 52  3   5
 53  3   5
```

## Step 4. Selection
Next, although the heritability is set to `0.5`, we can reassign it in `:h2` with a new value (or a vector for multiple traits). The argument `:weights` allows us to weight two traits differently in the selction. This example we select `3` sires from the `15` male progenies and `10` dams from the `15` female progenies.

```jldoctest
julia> args_select = Dict(:h2     => [.8, .5],
                          :weights=> [.6, .4])

julia> sires       = select(males, 3; args_select...)
[ Info: --------- Selection Summary ---------
[ Info: Select 3 individuals out of 15 individuals
[ Info: Selection differential (P): [0.794 0.972]
[ Info: Selection response     (G): [0.68 1.078]
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     0.25  0.0
└     0.0   1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (3 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [2.591 0.597]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.972 0.095]

julia> dams        = select(females, 10; args_select...)
[ Info: --------- Selection Summary ---------
[ Info: Select 10 individuals out of 15 individuals
[ Info: Selection differential (P): [0.328 0.17]
[ Info: Selection response     (G): [0.37 0.237]
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     0.25  0.0
└     0.0   1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (10 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [2.883 -0.093]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.56 1.007]
```

## Step 5. Expand to Multiple Generations
We can expand the described `mate()` and `select()` to `:n_gens` generations. By assigning `:n_select_males` and `:n_select_females` to specify how many progenies will be passed to the next generation.

### Expand to multiple generations
```jldoctest
julia> args_breed  = Dict(:n_gens           => 5,
                          :n_select_males   => 3,
                          :n_select_females => 20)

julia> sires, dams = breed(sires, dams; args_breed..., args_mate..., args_select...)
[ Info: Gen 0 -> Mean of BVs: [2.816 0.066], Variance of BVs: [0.598 0.863]
[ Info: Gen 1 -> Mean of BVs: [3.063 0.429], Variance of BVs: [0.998 1.014]
[ Info: Gen 2 -> Mean of BVs: [3.596 0.483], Variance of BVs: [0.481 0.524]
[ Info: Gen 3 -> Mean of BVs: [4.001 0.84], Variance of BVs: [0.138 0.574]
[ Info: Gen 4 -> Mean of BVs: [4.18 0.965], Variance of BVs: [0.065 0.807]
[ Info: Gen 5 -> Mean of BVs: [4.18 1.182], Variance of BVs: [0.065 0.515]
([ Info: Cohort (3 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [4.192 1.523]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.03 0.83]
, [ Info: Cohort (15 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [4.177 1.114]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.074 0.477]
```

```jldoctest
julia> summary(sires + dams)
Dict{String,Any} with 3 entries:
  "mu_g"  => [4.18 1.182]
  "var_g" => [0.065 0.515]
  "n"     => 18
```

### Modularism
Codes below are equivalent to the `breed()` function
```jldoctest
julia>  for i in 1:5
            males, females = mate(sires, dams; args_mate...)
            sires          = select(males, 3; args_select...)
            dams           = select(females, 20; args_select...)
        end

```
```jldoctest
julia> summary(sires + dams)
Dict{String,Any} with 3 entries:
  "mu_g"  => [4.259 2.05]
  "var_g" => [0.024 0.454]
  "n"     => 18
```

## Complete Codes
```julia
# Load XSim
using XSim
import Random
Random.seed!(95616)

# Build genome and phenome
build_genome(n_chr=10, n_marker=100)
build_phenome([3, 8])

# Initialize founders
n_sires = 3
n_dams  = 20
sires   = Founders(n_sires)
dams    = Founders(n_dams)

# Define parameters
args     = Dict(# mating
                :nA               => 3,
                :nB_per_A         => 5,
                :n_per_mate       => 2,
                :ratio_malefemale => 1.0,
                # selection
                :h2               => [.8, .5],
                :weights          => [.6, .4],
                # breeding
                :n_gens           => 5,
                :n_select_males   => 3,
                :n_select_females => 20)

# Breeding program
sires, dams   = breed(sires, dams; args...)
```