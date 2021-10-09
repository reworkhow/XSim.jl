# Build Genome

`build_genome()` defines the genetic information including numbers of chromosomes and loci, genetic position, physical position, minor allele frequency, mutation rates, and genotyping error rates.

```@contents
Pages = ["build_genome.md"]
Depth = 4
```

## Quick Start
Quick setup by assigning number of `markers` and `chromosomes`.

    build_genome(;n_loci ::Int64=-1,
                  n_chr    ::Int64=10,
                  species  ::String="none",
                  args...)

### Arguments
- `n_loci` : Number of simulated markers for each chromosome
- `n_chr`: Number of simulated chromosome with length of 100 centimorgan
- `species` : Infer genetic position (Morgan) by pre-load linkage maps. Available species are: ['pig', 'cattle', 'maize', 'rice'].

### Examples
```jldoctest
julia> build_genome(n_chr    = 2,
                    n_loci = 10000)

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

## Build by a File or a `DataFrame`
Define genome by providing a formatted dataframe or a path to the file.

    build_genome(dt      ::DataFrame;
                 species ::String="none",
                 args...)

    build_genome(filename::String;
                 species ::String="none",
                 args...)

### Arguments
- `dt` : A `DataFrame` with required columns of `chr` and `cM`, or `chr` and `bp` if `species` is provided for the inference.
- `filename` : A filepath to the file containing genome information.
- `species` : Adjust genetic position (Morgan) by pre-load linkage maps, available species are: ["cattle", and "pig"]

### Example of the `DataFrame`
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

### Examples
The preloaded file can be located through:
```jldoctest
# Filepath to the preloaded map file
julia> filepath = PATH("map")
```

Build genome by a filepath
```jldoctest
julia> build_genome(filepath;
                    rate_mutation=0.005, rate_error=0.01)
```

or by a dataframe
```jldoctest
julia> using DataFrames
julia> data = CSV.read(filepath, DataFrame)
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
julia> build_genome(filepath; species="cattle")

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
## Build by Pre-Loaded Genome
Define genome by pre-loaded genome from reference.

    build_genome(;species :: String="none")

### Arguments
- `species` : Adjust genetic position (Morgan) by pre-load linkage maps, available species are: ["cattle", and "pig"]

### Examples
```jldoctest
julia> build_genome(species="pig")
[ Info: Tortereau,F. et al. (2012) A high density recombination map of the pig reveals a correlation between sex-specific recombination and GC content. BMC Genomics, 13, 586.
[ Info: Reference Genome : Sscrofa 10.2
[ Info: SNP Chip         : PorcineSNP60 BeadChip
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 18
[ Info: 
[ Info: Chromosome Length (cM):
[ Info: [98.0, 94.7, 96.8, 92.2, 89.3, 124.1, 112.8, 94.8, 95.4, 84.4, 64.6, 77.0, 97.6, 106.6, 93.5, 66.5, 53.9, 49.1]
[ Info: 
[ Info: Number of Loci        : 45292
[ Info: [6580, 2356, 1938, 3682, 2217, 1766, 3489, 2100, 2538, 1281, 1805, 1072, 3529, 4053, 2612, 1513, 1646, 1115]
[ Info: 
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info: 
```

## Build by Explicit Definition
Define genome by providing genetic information of each loci explicitly.

    build_genome(chromosome    ::Array{Int64,   1},
                 bp            ::Array{Int64,   1},
                 cM            ::Array{Float64, 1},
                 maf           ::Array{Float64, 1};
                 rate_mutation ::Float64=0.0,
                 rate_error    ::Float64=0.0)

### Arguments
- `chromosome` : Chromosome codes
- `bp` : Physical positions
- `cM` : Genetic positions
- `maf` : Minor allele frequencies
- `rate_mutation` : Mutation rate
- `rate_error` : Error rate of genotyping

### Examples
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
```
