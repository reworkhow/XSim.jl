"""
# Quick Start
Quick setup by assigning number of `QTL`.

    build_phenome(n_qtls ::Int64;
                  args...)

## Arguments
- `n_qtls` : Number of simulated the `QTLs`.
- `vg` : Genetic (co)variances of the `QTLs`.
- `ve` : Residual (co)variances of the `QTLs`.
- `h2` : Heritability of simulated traits. This will define the residual (co)variances.

## Examples
Single trait
```jldoctest
julia> build_phenome(10)
--------- Phenome Summary ---------
Number of Traits      : 1
Heritability (h2)     : [0.5]
Number of QTLs        : [10;;]
Genetic (Co)variance
 1.0
Residual (Co)variance
 1.0
QTL Effects (Only first 30 QTLs are shown)
 -1.112
 -0.724
 -0.29
 -0.737
 -0.216
 -1.267
 -1.566
  1.603
 -0.537
 -0.156
```

Multi-trait with additional information
```jldoctest

julia> build_phenome(5;
                     n_traits = 2,
                     vg = [1.0 0.5
                           0.5 1.0],
                     h2 = [.1, .8])
--------- Phenome Summary ---------
Number of Traits      : 2
Heritability (h2)     : [0.1, 0.8]
Number of QTLs        : [5 5]
Genetic (Co)variance
 1.0  0.5
 0.5  1.0
Residual (Co)variance
 9.0  0.0
 0.0  0.25
QTL Effects (Only first 30 QTLs are shown)
  0.255   1.552
  0.595   0.648
 -0.146  -1.169
 -0.35   -0.659
 -0.342   0.064
```
────────────────────────────────────────────────────────────────
# Define phenome by a file or a DataFrame
Define genome by providing a formatted dataframe or a path to the file.

    build_phenome(dt       ::DataFrame; args...)
    build_phenome(filename ::String;    args...)

## Arguments
- `dt` : A `DataFrame` with required columns of `eff_` prefixed specifying marker effects.
- `filename` : A filepath to the file containing phenome information.

## Example of the `DataFrame`
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

## Examples
By a filepath
```jldoctest
julia> path = PATH("map")
julia> build_phenome(path, h2 = [0.3, 0.5])
```
or a dataframe
```jldoctest
julia> using DataFrames
julia> data = CSV.read(path, DataFrame)
julia> build_phenome(data, h2 = [0.3, 0.5])

Vg is not provided and is set to the covariance of the QTL effects
Ve is not provided and is set based on the Vg and h2
--------- Phenome Summary ---------
Number of Traits      : 2
Heritability (h2)     : [0.3, 0.5]
Number of QTLs        : [2 2]
Genetic (Co)variance
  0.009  -0.005
 -0.005   0.11
Residual (Co)variance
 0.021  0.0
 0.0    0.11
QTL Effects (Only first 30 QTLs are shown)
 0.1   0.0
 0.0  -0.3
 0.2   0.0
 0.0   0.5
```

────────────────────────────────────────────────────────────────
# Explicitly define phenome by providing a matrix of all loci effects.

    build_phenome(effects ::Union{Array{Float64}, SparseMatrixCSC}; args...)


## Arguments

- `effects` : A matrix storing marker effects with the dimension of individuals by markers.

## Examples

```jldoctest
julia> effects = [0.1  0.0
                  0.0 -0.3
                  0.2  0.0
                  0.0  0.5]
julia> build_phenome(effects)
Vg is not provided and is set to the covariance of the QTL effects
Heritability is not provided and is set to 0.5
Ve is not provided and is set based on the Vg and h2
--------- Phenome Summary ---------
Number of Traits      : 2
Heritability (h2)     : [0.5, 0.5]
Number of QTLs        : [2 2]
Genetic (Co)variance
  0.009  -0.005
 -0.005   0.11
Residual (Co)variance
 0.009  0.0
 0.0    0.11
QTL Effects (Only first 30 QTLs are shown)
 0.1   0.0
 0.0  -0.3
 0.2   0.0
 0.0   0.5
```

It's also possible to add additional information such as heritability.

```jldoctest
julia> build_phenome(effects,
                     vg=[1.0 0.5
                         0.5 1.0])
Heritability is not provided and is set to 0.5
Ve is not provided and is set based on the Vg and h2
--------- Phenome Summary ---------
Number of Traits      : 2
Heritability (h2)     : [0.5, 0.5]
Number of QTLs        : [2 2]
Genetic (Co)variance
 1.0  0.5
 0.5  1.0
Residual (Co)variance
 1.0  0.0
 0.0  1.0
QTL Effects (Only first 30 QTLs are shown)
 0.1   0.0
 0.0  -0.3
 0.2   0.0
 0.0   0.5
```
"""
function build_phenome(
    n_qtls::Int64;
    n_traits::Int64=1,
    vg::Union{Missing,Float64,Array{Float64}}=missing,
    ve::Union{Missing,Float64,Array{Float64}}=missing,
    vp::Union{Missing,Float64,Array{Float64}}=missing,
    h2::Union{Missing,Float64,Array{Float64}}=missing)

    if n_traits == 1
        LOG("If it is for a multi-trait simulation, make sure `n_traits` is set correctly.", "warn")
    end

    # collect users' inputs
    h2 = handle_h2(h2, n_traits)
    vg, ve, vp, h2 = handle_variance(vg, ve, vp, h2)

    # sample QTL effects using MVN
    QTL_effects = sample_qtls(n_qtls, n_traits, vg)

    # all loci matrix
    n_loci = GLOBAL("n_loci")
    effects = spzeros(n_loci, n_traits)
    idx_qtls = sample(1:n_loci, n_qtls, replace=false)
    effects[idx_qtls, :] = QTL_effects

    # build phenome
    build_phenome(effects, vg=vg, ve=ve, h2=h2)
end

"""
Explicitly define phenome by providing a matrix of all loci effects.
"""
function build_phenome(
    effects::Union{Array{Float64},SparseMatrixCSC};
    vg::Union{Missing,Float64,Array{Float64}}=missing,
    ve::Union{Missing,Float64,Array{Float64}}=missing,
    h2::Union{Missing,Float64,Array{Float64}}=missing)

    n_loci, n_traits = size(effects)
    h2 = handle_h2(h2, n_traits)

    if ismissing(vg)
        vg = cov(effects)
        LOG("Vg is not provided and is set to the covariance of the QTL effects")
    end
    if ismissing(h2)
        h2 = ones(n_traits) * 0.5
        LOG("Heritability is not provided and is set to 0.5")
    end
    if ismissing(ve)
        ve = get_Ve(n_traits, vg, h2)
        LOG("Ve is not provided and is set based on the Vg and h2")
    end

    SET("n_traits", n_traits)
    SET("effects", SparseArrays.SparseMatrixCSC(effects))
    SET("Vg", round.(vg, digits=3))
    SET("Ve", round.(ve, digits=3))
    SET("h2", round.(h2, digits=3))

    # Summary
    summary_phenome()
end

```
Load dataframe to define effects
```
function build_phenome(
    dt::DataFrame; args...)

    QTL_effects = from_dt_to_eff(dt)
    build_phenome(QTL_effects; args...)
end

```
Load file to define effects
```
function build_phenome(
    filename::String; args...)

    build_phenome(CSV.read(filename, DataFrame); args...)
end

function summary_phenome()
    n_traits = GLOBAL("n_traits")
    n_qtls = sum(GLOBAL("effects") .!= 0, dims=1)
    Vg = GLOBAL("Vg")
    Ve = GLOBAL("Ve")
    h2 = GLOBAL("h2")

    if !GLOBAL("silent")
        LOG("--------- Phenome Summary ---------")
        LOG("Number of Traits      : $n_traits")
        LOG("Heritability (h2)     : $h2")
        LOG("Number of QTLs        : $n_qtls")
        LOG("Genetic (Co)variance")
        Base.print_matrix(stdout, Vg)
        LOG("")
        LOG("Residual (Co)variance")
        Base.print_matrix(stdout, Ve)
        LOG("")
        LOG("QTL Effects (Only first 30 QTLs are shown)")
        effects = round.(GLOBAL("effects_QTLs") |> Matrix, digits=3)
        n_qtls = size(effects, 1) # n_rows
        if n_qtls >= 30
            Base.print_matrix(stdout, effects[begin:30, :])
        else
            Base.print_matrix(stdout, effects)
        end
        LOG("")
    end
end

function from_dt_to_eff(dt::DataFrame)
    columns = names(dt)
    idx_eff = [occursin("eff", s) for s in columns]
    return Matrix(dt[:, idx_eff])
end




