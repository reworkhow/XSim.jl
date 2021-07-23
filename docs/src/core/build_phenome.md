# Build Phenome
`build_phenome` defines the phenomics information including `QTL effects` and `heritability`.

## Quick Start
Quick setup by assigning number of `QTL`.

    build_phenome(n_qtls ::Union{Array{Int64, 1}, Int64};
                  args...)

### Arguments
- `n_qtls` : Number of simulated `QTLs`. It can be an array of integers for multiple traits.
- `vg` : Genetic (co)variances of `QTLs`.
- `h2` : Heritability of simulated traits. This will define the residual (co)variances.

### Examples
Single trait
```jldoctest
julia> build_phenome(10)
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

Multi-trait with additional information
```jldoctest
julia> build_phenome([10, 15];
                     vg = [1 .5
                          .5  1],
                     h2 = [.3, .8])
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
[ Info: Number of QTLs        : [10 25]
```

## Define phenome by a file or a DataFrame
Define genome by providing a formatted dataframe or a path to the file.

    build_phenome(dt        ::DataFrame; args...)
    build_phenome(filename  ::String;    args...)

### Arguments
- `dt` : A `DataFrame` with required columns of `eff_` prefixed specifying marker effects.
- `filename` : A filepath to the file containing phenome information.

### Example of the `DataFrame`
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

### Examples
By a filepath
```jldoctest
julia> build_phenome("path/map.csv", h2 = [0.3, 0.5])
```
or a dataframe
```jldoctest
julia> using DataFrames
julia> data = CSV.read("path/map.csv", DataFrame)
julia> build_phenome(data, h2 = [0.3, 0.5])

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


## Define phenome by a matrix of QTL effects

    build_phenome(QTL_effects ::Union{Array{Float64}, SparseMatrixCSC}; args...)


### Arguments

- `QTL_effects` : A matrix storing marker effects with the dimension of individuals by markers.

### Examples

```jldoctest
julia> effects = [0.1  0.0
                  0.0 -0.3
                  0.2  0.0
                  0.0  0.5]
julia> build_phenome(effects)
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
[ Info: Number of QTLs        : [2 2]
```
It's also possible to add additional information such as heritability.

```jldoctest
julia> build_phenome(effects, h2=[0.1, 0.8])
[ Info: --------- Phenome Summary ---------
[ Info: Number of Traits      : 2
[ Info: Heritability (h2)     : [0.1, 0.8]
┌ Info: 
│   Genetic_Variance =
│    2×2 Array{Float64,2}:
│     1.0  0.0
└     0.0  1.0
┌ Info: 
│   Residual_Variance =
│    2×2 Array{Float64,2}:
│     9.0  0.0
└     0.0  0.25
[ Info: Number of QTLs        : [2 2]
```