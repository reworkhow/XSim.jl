# Cohort and Founders

!!! tip "Cohort and Founders"
    In the current version, `Cohort` works exactly the same as `Founders`.

# Initialize a cohort by population size
    Cohort(n::Int64=0)

## Arguments
- `n` : An integer to assign the population size.

## Examples
```jldoctest
julia> cohort = Cohort(5)
[ Info: Cohort (5 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [1.265 1.697]
[ Info:
[ Info: Variance of breeding values:
[ Info: [1.6 1.4]
```
──────────────────────────────────────────────────────────────
# Initialize a cohort by genotypes/haplotypes files
    Cohort(genetic_data ::Union{DataFrame, Array{Int64}}; args...)
    Cohort(filename     ::String; args...)

## Arguments
- `genetic_data` : A `dataframe`/`2D-array` that stores genotypes/haplotypes in the dimension of individuals by markers.
- `filename` : A `filepath` to a file storing genotypes/haplotypes data.
- `n` : Number of lines to be loaded from the file. The default value is `-1` and the entire file will be loaded.
- `random` : By default it's set to `true` to randomly select `n` lines (individuals) from the file to generate the cohort.
- `alter_maf` : It will update MAF based on the provided genotypes if it's set to `true` (default).

## Example of the `demo_genotypes.csv` and `demo_haplotypes.csv`
Both demo files store marker information for 5 individuals and 4 markers.
Use `DATA("demo_genotypes.csv")` to interact with demo files.
```
# demo_genotypes.csv
# rows: individuals, columns: markers
# Homozygote is coded as 0 and 2, otherwise is coded as 1
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

## Example
```jldoctest
# Load entire file
julia> cohort = Cohort("demo_haplotypes.csv")
julia> get_genotypes(cohort)
5×4 Array{Int64,2}:
 2  0  0  1
 2  0  2  0
 0  0  1  0
 1  1  0  2
 0  1  0  2

# Randomly load 3 individuals with a dataframe.
julia> data = DATA("demo_haplotypes.csv", header=false)
julia> cohort = Cohort(data, random=true, n=3)
julia> get_genotypes(cohort)
3×4 Array{Int64,2}:
 2  0  2  0
 0  1  0  2
 1  1  0  2

# Replace marker MAF by the provided file
julia> cohort = Cohort("demo_haplotypes.csv", alter_maf=true)
[ Info: MAF has been updated based on provided haplotypes/genotypes
[ Info: Cohort (5 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [1.418]
[ Info:
[ Info: Variance of breeding values:
[ Info: [2.012]
```
──────────────────────────────────────────────────────────────
# Functions that insepct `Cohort` properties:
All the listed functions can take a keyword argument `ID=true` to insert individuals' IDs as the first column.

## Genotypes
Genotype matirx in the dimension of `individuals` by `markers`
```jldoctest
julia> get_genotypes(cohort)
5×4 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 0  0  1  0
 2  0  2  0
 2  0  0  1
 0  1  0  2
 1  1  0  2
```
## QTLs
QTLs matirx in the dimension of `individuals` by `markers`
```jldoctest
julia> get_QTLs(cohort)
5×3 Array{Int64,2}:
 2  2  0
 0  0  2
 0  1  0
 1  0  2
 2  0  1
```
## Breeding values
Breeding values in the dimenstion `individuals` by `traits`
```jldoctest
julia> get_BVs(cohort)
5×2 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:
 1.26491   0.0
 3.79473   0.0
 1.26491   1.21268
 0.0       1.69775
 0.632456  1.69775
```
## Pedigree
Pedigree matrix, listed columns are in the order of individuals' ID, sire ID, and dam ID.
```jldoctest
julia> get_pedigree(cohort)
5×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
1  0  0
2  0  0
3  0  0
4  0  0
5  0  0
```

## Minor Allele Frequencies (MAF)
In the case where we have 3 QTLs out of 4 markers, we want to compare their allel frequencies.

```jldoctest
julia> get_MAF(cohort)
4-element Array{Float64,1}:
 0.5
 0.2
 0.3
 0.5
```

## Phenotypes
Simulate cohort phenotypes based on the defined `phenome`. `h2` and `ve` can be assigned specifically for this one-time simulation.
```jldoctest
julia> get_phenotypes(cohort)
5×1 Array{Float64,2}:
  1.1126064336992942
 -0.8337021175232547
 -0.363621019381922
  4.042256656472762
  1.7828800511223049
```