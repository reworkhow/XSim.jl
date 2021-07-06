# Crossbreed
10 generations of individuals. Parallel purebred populations are simulated
as well as a crossbred population. 1 small 2 large pure breeds,
and a crossbred (X) population, for example.
Breed 1 has 50 males 500 females,
Breeds 2-3 have 100 males 2000 females at G0.

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
```julia;results="hidden";eval=false;
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

### Derive small breed (A)
```julia;results="hidden";eval=false;
n_sires       = 50
dams_per_sire = 10
n_dams        = n_sires * dams_per_sire
args          = Dict(# Mating
                     :n_per_shared     => dams_per_sire,
                     :n_per_mate       => 2,
                     :ratio_malefemale => 1,
                     # Selection
                     :h2               => [.8, .2],
                     :is_random        => false,
                     # Breeding
                     :n_gens           => 10,
                     :n_select_males   => n_sires)
# Breed A
sires_A         = Founders(n_sires)
dams_A          = Founders(n_dams)
sires_A, dams_A = breed(sires_A, dams_A; args...)
```

### Derive large breeds (B and C)
```julia;results="hidden";eval=false;
# Large breeds
n_sires        = 100
dams_per_sire  = 20
n_dams         = n_sires * dams_per_sire
args[:n_per_shared]   = dams_per_sire
args[:n_select_males] = n_sires

# Breed B
sires_B         = Founders(n_sires)
dams_B          = Founders(n_dams)
sires_B, dams_B = breed(sires_B, dams_B; args...)

# Breed C
sires_C         = Founders(n_sires)
dams_C          = Founders(n_dams)
sires_C, dams_C = breed(sires_C, dams_C; args...)
```


### Rotational breeding
```julia;results="hidden";eval=false;
# Rotation parameters
args = Dict(:n_pop            => 2000,
            :n_per_mate       => 2,
            :ratio_malefemale => 1)
# Rotation (G1)
males_G1, females_G1 = mate(sires_B, dams_C; args...)

# Rotation (G2)
males_G2, females_G2 = mate(sires_A, females_G1; args...)

# Rotation (G3)
males_G3, females_G3 = mate(sires_C, females_G2; args...)
```
###