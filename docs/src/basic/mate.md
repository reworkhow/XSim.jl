# Mate

## Random mating (Default)
```julia
cohort_A = Founders(5)
cohort_B = Founders(10)
# Random mating scheme
args = Dict(:nA               => cohort_A.n,
            :nB_per_A         => 1,
            :replace_A        => false,
            :replace_B        => false,
            :n_per_mate       => 1,
            :ratio_malefemale => 1)
male, female = mate(cohort_A, cohort_B; args...)
# Equivalent results without any argument
male, female = mate(cohort_A, cohort_B)
# Equivalent results by specifying scheme argument
male, female = mate(cohort_A, cohort_B; scheme="random")
# Equivalent results with overloaded operator '*'
male, female = cohort_A * cohort_B
```
```
[ Info: --------- Mating Summary ---------
[ Info: Generated 5 individuals from 5 cohort_A individuals
[ Info: Every cohort_A individuals mates with 1 individuals from cohort_B
[ Info: And each mating reproduces 1 progenies
[ Info: 
[ Info: --------- Offsprings Summary ---------
([ Info: Cohort (2 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.403 0.411]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [1.434 0.004]
, [ Info: Cohort (3 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.973 0.405]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [2.174 2.119]
)
```

Inspect the pedigree to verify the mating behavior
```julia
get_pedigree(male + female)
```
```
5×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 19  1   8
 16  2   6
 17  3  10
 20  4  15
 18  5  14
```

## Another example with larger smaple 
```julia
sires = Founders(10)
dams  = Founders(200)
args = Dict(:nA               => 5,
            :nB_per_A         => 10,
            :replace_A        => false,
            :replace_B        => false,
            :n_per_mate       => 1,
            :ratio_malefemale => 1)
male, female = mate(sires, dams; args...)
```
```
[ Info: --------- Mating Summary ---------
[ Info: Generated 50 individuals from 5 cohort_A individuals
[ Info: Every cohort_A individual mates with 10 individuals from cohort_B
[ Info: And each mating reproduces 1 progenies
[ Info:
[ Info: --------- Offsprings Summary ---------
([ Info: Cohort (25 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [2.15 0.252]
[ Info:
[ Info: Variance of breeding values:
[ Info: [1.0 0.968]
, [ Info: Cohort (25 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [2.378 0.349]
[ Info:
[ Info: Variance of breeding values:
[ Info: [1.176 0.863]
)
```

```julia
sires = Founders(2)
dams  = Founders(5)
args = Dict(:nA        => sires.n,
            :nB_per_A  => dams.n,
            :replace_A => false,
            :replace_B => false,
            :n_per_mate=> 1,
            :ratio_malefemale=> 1)
male, female = mate(sires, dams; args...)
# Equivalent results by specifying scheme argument __
male, female = mate(sires, dams; scheme = "diallel cross")
```

Inspect the pedigree to verify the mating behavior
```julia
get_pedigree(male + female)
```
```
10×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 12  2  7
 10  2  6
 11  2  4
 14  1  3
 15  1  5
  9  2  5
 13  1  6
 17  1  4
  8  2  3
 16  1  7
```

## Selfing
```julia
parents = Founders(5)
args = Dict(:nA               => 3,
            :replace_A        => false,
            :n_per_mate       => 5,
            :scheme       => "selfing")
progenies = mate(parents; args...)
```
```
[ Info: --------- Mating Summary ---------
[ Info: Generated 250 individuals from 5 cohort_A individuals
[ Info: Each selfing reproduces 50 progenies
[ Info: 
[ Info: --------- Offsprings Summary ---------
[ Info: Cohort (250 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [-0.824 0.225]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [0.311 0.709]
```

Inspect the pedigree to verify the mating behavior
```julia
get_pedigree(progenies)
```
```
15×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
  6  4  4
  7  4  4
  8  4  4
  9  4  4
 10  4  4
 11  1  1
 12  1  1
 13  1  1
 14  1  1
 15  1  1
 16  5  5
 17  5  5
 18  5  5
 19  5  5
 20  5  5
```
