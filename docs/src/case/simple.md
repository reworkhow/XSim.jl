# Simple Case
```julia; echo=false;
include("../src/XSim.jl")
using .XSim
```

### Step 0. Load XSim and set random seed
```julia; eval=false;
using XSim
import Random
Random.seed!(95616)
```

### Step 1. Setup genome and phenome information
```julia;
build_demo()
```
```
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 10
[ Info: 
[ Info: Chromosome Length (cM): 750.0
[ Info: [75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0]
[ Info: 
[ Info: Number of Loci        : 50
[ Info: [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
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
[ Info: Number of QTLs        : [3 8]
```

### Step 2. Initialize founders
```julia
n_sires = 5
sires   = Founders(n_sires)
```
```
[ Info: Cohort (5 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.177 1.121]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.867 0.91]
```
```julia
n_dams  = 3
dams    = Founders(n_dams)
```
```
[ Info: Cohort (3 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [0.464 0.181]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.899 0.376]
```

### Step 3. Mate
```julia; results="hidden"
args_mate     = Dict(:n_per_shared => n_dams,
                     :n_per_mate   => 2)
progenies     = mate(sires, dams; args_mate...)
```
```
[ Info: --------- Mating Summary ---------
[ Info: Generate 30 individuals from 5 shared breeds
[ Info: Every shared breeds mates with 3 breeds
[ Info: And each mating reproduces 2 progenies
[ Info: 
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (30 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [0.611 0.777]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.989 1.309]
```


### Step 4. Select
```julia
args_select   = Dict(:h2 => [.5, .5])
progenies_sel = select(progenies, 10; args_select...)
```
```
[ Info: --------- Selection Summary ---------
[ Info: Select 10 individuals out of 30 individuals
[ Info: Selection differential (P): [0.764 0.845]
[ Info: Selection response     (G): [0.586 0.73]
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (10 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.193 1.612]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [1.057 0.874]
```

### Step 5. Breed
##### Expand to multiple generations
```julia; results="hidden"
args_breed  = Dict(:n_gens   => 5,
                   :n_select => 10)
sires, dams = breed(sires, dams; args_breed..., args_mate..., args_select...)
```
```
[ Info: Gen 0 -> Mean of BVs: [0.909 0.769], Variance of BVs: [0.888 0.864]
[ Info: Gen 1 -> Mean of BVs: [1.363 1.307], Variance of BVs: [0.682 0.691]
[ Info: Gen 2 -> Mean of BVs: [1.987 1.792], Variance of BVs: [0.295 0.502]
[ Info: Gen 3 -> Mean of BVs: [1.906 2.732], Variance of BVs: [0.768 0.542]
[ Info: Gen 4 -> Mean of BVs: [2.665 2.807], Variance of BVs: [0.464 0.664]
[ Info: Gen 5 -> Mean of BVs: [3.091 3.05], Variance of BVs: [0.053 0.468]
```
```julia; eval=false;
summary(sires + dams)
```
```
Dict{String,Any} with 3 entries:
  "mu_g"  => [3.091 3.05]
  "var_g" => [0.053 0.468]
  "n"     => 20
```
##### Modularism of XSim
```julia; eval=false;
for i in 1:5
    progenies = mate(sires, dams; args_mate...)
    progenies = select(progenies, 10; args_select...)
    sires, dams = progenies, progenies
end
```
```julia; eval=false;
summary(sires + dams)
```
```
Dict{String,Any} with 3 entries:
  "mu_g"  => [3.091 3.05]
  "var_g" => [0.053 0.468]
  "n"     => 20
```


### Complete code
```julia; eval=false;
# Load XSim
using XSim
import Random
Random.seed!(95616)

# Build genome and phenome
build_demo()

# Initialize founders
n_sires = 5
n_dams  = 3
sires   = Founders(n_sires)
dams    = Founders(n_dams)

# Define parameters
args     = Dict(# mating
                :n_per_shared => n_dams,
                :n_per_mate   => 2,
                # selection
                :h2           => [.5, .5],
                # breeding
                :n_gens       => 5,
                :n_select     => 10)

# Breeding program
sires, dams   = breed(sires, dams; args...)
```