# Select

## Single Trait Selection
We can build genome from demo

```julia
# load demo genome and phenome
build_demo()
# override demo phenome with single traits controlled by 50 QTLs
build_phenome(50)
# inspect the settings
summary()
# initialize a cohort with 100 individuals
cohort = Cohort(100)
```
```
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

### Select 30 individuals
```julia
# Assign by integer
cohort_s = select(cohort, 30)
# it's equivalent to assigning the selected ratio
cohort_s = select(cohort, 0.3)
```
```
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

### Re-Assign Heritability or Residual Covariance
It's possible to re-assign heritability or residuals to
simulate different selection environment

```julia
# Assign heritability h2 = 0.1
progenies = select(cohort, 30, h2=0.1)
# Equivalent to assigning residual as 9.0
progenies = select(cohort, 30, ve=9.0)
```
```
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

### Negative Selection
Set `is_positive=false` to rank individuals in ascending order
```julia
progenies = select(cohort, 30, is_positive=false)
```
```
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

### Random Selection
```julia
progenies = select(cohort, 30, is_random=true)
```
```
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

### Selection wiht Multiple Parameters
It's possible specify multiple parameters described above in one selection.
User can either enclose parameters as keyword arguments, or pass them through a dictionary object

```julia
# Pass them as keyword arguments
progenies = select(cohort, 30, h2=0.3, is_positive=false)
# Or pass them via a dictionary
args = Dict(:h2=>0.3,
            :is_positive=>false)
progenies = select(cohort, 30; args...)
```
```
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

## Multi-Trait Selection
We can build genome from demo

```julia
# load demo genome and phenome with multiple traits
build_demo()
# inspect the settings
summary()
# initialize a cohort with 100 individuals
cohort = Cohort(100)
```
```
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

### Assign Heritabilities for Multiple Traits
```julia
progenies = select(cohort, 30, h2=[0.9, 0.3])
```
```
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

### Assign Trait Correlations via Residual Covariance
```julia
progenies = select(cohort, 30, ve=[1   0.3
                                   0.3   1])
```
```
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

### Derive Selection Index for Multiple Traits
Assigning a vector to the parameter `weights` to derive a index consisting a linear combintation of the weights and the phenotypes. 
In this example, we demonstrate two traits with the heritability of 0.3 and 0.8, respectively.
And we can select traits with more weight on the second trait which is more heritable, and negatively select the first trait.

```julia
progenies = select(cohort, 30, h2=[.3, .8], weights=[-0.1, 0.9])
```
```
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

