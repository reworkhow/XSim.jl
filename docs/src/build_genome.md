# Build Genome

## Data

```
# map.csv
id,chr,bp,cM,MAF,eff_1,eff_2
snp 1,1,1818249,50.8,0.5,1.5,2.8
snp 2,1,6557697,80.3,0.5,0.0,0.0
snp_3,2,2298800,39.2,0.5,0.0,0.0
snp 4,2,5015698,66.3,0.5,0.0,0.0
```

## Quick start
```julia
build_genome(n_chr    = 2,
             n_marker = 10000)
```

## By a file
```julia
build_genome("map.csv")
```

## By a dataframe
```julia
using DataFrames
data = CSV.read("map.csv", DataFrame)
build_genome(data)
```

## Load with pre-load reference
```julia
build_genome("map.csv"; species="cattle")
build_genome(data;      species="cattle")
```


## By a manual inputs
```julia
ch  = [1,    1,     2,    2,    2]
bp  = [130,  205,   186,  503,  780]
cM  = [85.7, 149.1, 37.4, 83.6, 134.3]
maf = [0.5,  0.5,   0.5,  0.5,  0.5]
build_genome(ch  = [1,    1,     2,    2,    2],
             bp  = [130,  205,   186,  503,  780],
             cM  = [85.7, 149.1, 37.4, 83.6, 134.3],
             maf = [0.5,  0.5,   0.5,  0.5,  0.5])
```

