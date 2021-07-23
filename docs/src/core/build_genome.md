# Build Genome

`build_genome()` defines the genetic information including numbers of chromosomes and loci, genetic position, physical position, minor allele frequency, mutation rates, and genotyping error rates.

# Quick Start
Quick setup by assigning number of `markers` and `chromosomes`.

    build_genome(;n_marker ::Int64=-1,
                  n_chr    ::Int64=10,
                  species  ::String="none",
                  args...)

## Arguments
- `n_marker` : Number of simulated markers for each chromosome
- `n_chr`: Number of simulated chromosome with length of 100 centimorgan
- `species` : Infer genetic position (Morgan) by pre-load linkage maps, available species are: ["cattle", and "pig"]

## Examples
```jldoctest
julia> build_genome(n_chr    = 2,
                    n_marker = 10000)

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


# Define genome by a file or a DataFrame
Define genome by providing a formatted dataframe or a path to the file.

    build_genome(dt      ::DataFrame;
                 species ::String="none",
                 args...)

    build_genome(filename::String;
                 species ::String="none",
                 args...)

## Arguments
- `dt` : A `DataFrame` with required columns of `chr` and `cM`, or `chr` and `bp` if `species` is provided for the inference.
- `filename` : A filepath to the file containing genome information.
- `species` : Adjust genetic position (Morgan) by pre-load linkage maps, available species are: ["cattle", and "pig"]

## Example of the `DataFrame`
```
4×7 DataFrame
 Row │ id      chr    bp       cM       MAF      eff_1    eff_2
     │ String  Int64  Int64    Float64  Float64  Float64  Float64
─────┼────────────────────────────────────────────────────────────
   1 │ snp_1       1  1818249     50.8      0.5      0.1      0.0
   2 │ snp_2       1  6557697     80.3      0.5      0.0      0.0
   3 │ snp_3       2  2298800     39.2      0.5      0.2      0.0
   4 │ snp_4       2  5015698     66.3      0.5      0.0      0.5
```

## Examples
By a filepath
```jldoctest
julia> build_genome("path/map.csv";
                    rate_mutation=0.005, rate_error=0.01)
```

or a dataframe
```jldoctest
julia> using DataFrames
julia> data = CSV.read("path/map.csv", DataFrame)
julia> build_genome(data;
                    rate_mutation=0.005, rate_error=0.01)

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

Use cattle genome as reference to infer the genetic positions
```jldoctest
julia> build_genome("path/map.csv"; species="cattle")

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

# Explict Definition
Define genome by providing genetic information of each loci explicitly.

    build_genome(chromosome    ::Array{Int64,   1},
                 bp            ::Array{Int64,   1},
                 cM            ::Array{Float64, 1},
                 maf           ::Array{Float64, 1};
                 rate_mutation ::Float64=0.0,
                 rate_error    ::Float64=0.0)

## Arguments
- `chromosome` : Chromosome codes
- `bp` : Physical positions
- `cM` : Genetic positions
- `maf` : Minor allele frequencies
- `rate_mutation` : Mutation rate
- `rate_error` : Error rate of genotyping

## Examples
```jldoctest
julia> ch  = [1,    1,     2,    2,    2]
julia> bp  = [130,  205,   186,  503,  780]
julia> cM  = [85.7, 149.1, 37.4, 83.6, 134.3]
julia> maf = [0.5,  0.5,   0.5,  0.5,  0.5]
julia> build_genome(ch, bp, cM, maf)

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