# Build Genome
`build_genome` defines the genetic information including numbers of chromosomes and loci, genetic position, physical position, and minor allele frequency of each locus, mutation rates, and genotyping error rates. 

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
build_genome(n_chr    = 2,
             n_marker = 10000)
```
```
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info: 
[ Info: Chromosome Length (cM): 200.0
[ Info: [100.0, 100.0]
[ Info: 
[ Info: Number of Loci        : 20000
[ Info: [10000, 10000]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 
```

## By a File or a Dataframe
```julia
# By a filepath
build_genome("map.csv";
             rate_mutation=0.005, rate_error=0.01)
# or by a dataframe directly, they are equivalent.
using DataFrames
data = CSV.read("map.csv", DataFrame)
build_genome(data;
             rate_mutation=0.005, rate_error=0.01)
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
[ Info: Genotyping Error      : 0.01
[ Info: Mutation Rate         : 0.005
[ Info: 
```

## Load with Pre-Load Reference
When the genetic distance (Morgan) is missing (or provided but not accurate) in the map, XSim can infer it by physical positions from published reference given the specified species.

```julia
# Use cattle genome as reference
build_genome("map.csv"; species="cattle")
```
```
[ Info: Arias,J.A. et al. (2009) A high density linkage map of the bovine genome. BMC Genetics, 10, 18.
[ Info: Reference Genome : Btau 4.0
[ Info: SNP Chip         : Affymetrix GeneChip Bovine Mapping 10K SNP kit

┌ Warning: The provided genetic distances will be replaced with ones infered from preloaded linkage maps
└ @ XSim ~/Dropbox/projects/XSim/src/objects/global.jl:118
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info:
[ Info: Chromosome Length (cM): 16.8
[ Info: [15.1, 1.7]
[ Info:
[ Info: Number of Loci        : 4
[ Info: [2, 2]
[ Info:
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info:
```

## By Manual Inputs
```julia
ch  = [1,    1,     2,    2,    2]
bp  = [130,  205,   186,  503,  780]
cM  = [85.7, 149.1, 37.4, 83.6, 134.3]
maf = [0.5,  0.5,   0.5,  0.5,  0.5]
build_genome(ch, bp, cM, maf)
```
```
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info: 
[ Info: Chromosome Length (cM): 283.4
[ Info: [149.1, 134.3]
[ Info: 
[ Info: Number of Loci        : 5
[ Info: [2, 3]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 
```

