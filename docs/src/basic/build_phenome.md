# Build Phenome
`build_phenome` defines the genetic information including numbers of chromosomes and loci, genetic position, physical position, and minor allele frequency of each locus, mutation rates, and genotyping error rates.

## Data
In this page, you will need a file `map.csv` to complete the demonstrated examples.

```
# map.csv
id,chr,bp,cM,MAF,eff_1,eff_2
snp_1,1,1818249,50.8,0.5,0.1,0.0
snp_2,1,6557697,80.3,0.5,0.0,-0.3
snp_3,2,2298800,39.2,0.5,0.2,0.0
snp_4,2,5015698,66.3,0.5,0.0,0.5
```
```julia
using DataFrames
data = CSV.read("map.csv", DataFrame)
```
```
4×7 DataFrame
 Row │ id      chr    bp       cM       MAF      eff_1    eff_2
     │ String  Int64  Int64    Float64  Float64  Float64  Float64
─────┼────────────────────────────────────────────────────────────
   1 │ snp 1       1  1818249     50.8      0.5      0.1      0.0
   2 │ snp 2       1  6557697     80.3      0.5      0.0     -0.3
   3 │ snp_3       2  2298800     39.2      0.5      0.2      0.0
   4 │ snp 4       2  5015698     66.3      0.5      0.0      0.5
```

## Quick Start
```julia
build_phenome(10)
```
```
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
[ Info: Number of QTLs        : [10]
```

## Multi-Trait with Asigned Genetic Variance and Heritability
```julia
build_phenome([50, 30];
              vg = [1 .5
                   .5  1],
              h2 = [.3, .8])
```
```
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 2
[ Info: Heritability (h2)     : [0.3, 0.8]
┌ Info:
│   Genetic_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.5
└     0.5  1.0
┌ Info:
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     2.33333  0.0
└     0.0      0.25
[ Info: Number of QTLs        : [50 80]
```

## By a File or a Dataframe
Similar to `build_genome()`, users can build the phenome with a file through a filepath or a dataframe.

```julia
build_phenome("demo_map.csv", h2 = [0.3, 0.5])
```
```
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 2
[ Info: Heritability (h2)     : [0.3, 0.5]
┌ Info: 
│   Genetic_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     2.33333  0.0
└     0.0      1.0
[ Info: Number of QTLs        : [2 2]
```

## Summary
```julia
summary()
```
```
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info: 
[ Info: Chromosome Length (cM): 146.6
[ Info: [80.3, 66.3]
[ Info: 
[ Info: Number of Loci        : 4
[ Info: [2, 2]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 2
[ Info: Heritability (h2)     : [0.3, 0.5]
┌ Info: 
│   Genetic_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     2.33333  0.0
└     0.0      1.0
[ Info: Number of QTLs        : [2 2]
```

