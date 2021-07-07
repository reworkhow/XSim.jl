# Founders

## Data
In this page, Files `demo_genotypes.csv`, `demo_haplotypes.csv` will be used in the demonstrated examples. Missing values can be represented by -1 or 9.

```
# Both files contains 5 individuals and 4 genetic markers

# demo_genotypes.csv
# rows: individuals, columns: markers
# homozygote is coded as 0 and 2, otherwise is coded as 1
2,0,0,1
0,0,1,0
0,1,0,2
1,1,0,2
2,0,2,0

# demo_haplotypes.csv
# rows: individuals, columns: alleles
# Reference allele is coded as 0, otherwise is coded as 1
1,1,0,0,0,0,1,0
0,0,0,0,1,0,0,0
0,0,0,1,0,0,1,1
1,0,1,0,0,0,1,1
1,1,0,0,1,1,0,0
```

## By Assigning a Population Size
```julia
founders = Founders(5)
```
```
[ Info: Cohort (5 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [1.265 1.697]
[ Info:
[ Info: Variance of breeding values:
[ Info: [1.6 1.4]
```

## By a genoetype or haplotypes file
```julia
# haplotypes
founder = Founders("demo_haplotypes.csv")
```
```
[ Info: MAF has been updated based on provided haplotypes/genotypes
[ Info: Cohort (5 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.391 0.922]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [2.08 0.747]
```

Or if users don't want to update the MAF by the data, `alter_maf=false` can achieve this purpose.
```julia
# genotypes
founders = Founders("demo_genotypes.csv", alter_maf=false)
```
```
[ Info: MAF has been updated based on provided haplotypes/genotypes
[ Info: Cohort (5 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [1.391 0.922]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [2.08 0.747]
```

Users can subset the population by assigning `n` and `random`. If `random=true`, `n` individuals are randomly selected from the file. Otherwise, only the first `n` individuals will be initialized.
```julia
founders = Founders("demo_genotypes.csv", n=3, random=true)
```
```
[ Info: MAF has been updated based on provided haplotypes/genotypes
[ Info: Cohort (3 individuals)
[ Info: 
[ Info: Mean of breeding values: 
[ Info: [2.108 0.404]
[ Info: 
[ Info: Variance of breeding values: 
[ Info: [2.133 0.49]
```

## Inspect the founders

### Genotypes
```julia
get_genotypes(founders)
```
```
5×4 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 0  0  1  0
 2  0  2  0
 2  0  0  1
 0  1  0  2
 1  1  0  2
```

### QTLs
```julia
get_QTLs(founders)
```
```
5×4 Array{Int64,2}:
 0  0  1  0
 2  0  2  0
 2  0  0  1
 0  1  0  2
 1  1  0  2
```

### Breeding Values
It's stored as a n by n_traits matrix
```julia
get_BVs(founders)
```
```
5×2 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:
 1.26491   0.0
 3.79473   0.0
 1.26491   1.21268
 0.0       1.69775
 0.632456  1.69775
```

### Pedigree
```julia
get_pedigree(founders)
```
```
5×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
1  0  0
2  0  0
3  0  0
4  0  0
5  0  0
```

### Minor Allele Frequencies (MAF)
In the case where we have 3 QTLs out of 4 markers, we want to compare their allel frequencies.

```julia
get_QTLs(founders) |> get_MAF
```
```
3-element Array{Float64,1}:
 0.5
 0.3
 0.5
```

```julia
get_MAF(founders)
```
```
4-element Array{Float64,1}:
 0.5
 0.2
 0.3
 0.5
```