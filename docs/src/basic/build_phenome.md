# Build Phenome
`build_phenome` defines the genetic information including numbers of chromosomes and loci, genetic position, physical position, and minor allele frequency of each locus, mutation rates, and genotyping error rates. 

## Data
In this page, you will need a file `map.file` to complete the demonstrated codeds.

```
# map.csv
id,chr,bp,cM,MAF,eff_1,eff_2
snp 1,1,1818249,50.8,0.5,1.5,2.8
snp 2,1,6557697,80.3,0.5,0.0,0.0
snp_3,2,2298800,39.2,0.5,0.3,0.0
snp 4,2,5015698,66.3,0.5,0.0,0.0
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
```julia
build_phenome("map.csv", h2 = 0.3)
```
```
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 2
[ Info: Heritability (h2)     : [0.3, 0.3]
┌ Info: 
│   Genetic_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     2.33333  0.0
└     0.0      2.33333
[ Info: Number of QTLs        : [2 2]
```
