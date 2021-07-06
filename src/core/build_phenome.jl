function build_phenome(QTL_effects  ::Union{Array{Float64}, SparseMatrixCSC};
                       vg           ::Union{Array{Float64}, Float64} =1.0,
                       h2           ::Union{Array{Float64 }, Float64}=0.5)

    # Assign number of traits
    SET("n_traits", size(QTL_effects)[2])

    # Assign Vg
    vg = handle_diagonal(vg, GLOBAL("n_traits"))
    vg = Symmetric(vg)
    SET("Vg"      , Array(vg))

    # Assign QTL effects
    SET("effects" , scale_effects(matrix(QTL_effects),
                                  GLOBAL("maf"),
                                  GLOBAL("Vg"),
                                  is_sparse=true))

    # Assign heritability
    SET("Ve"      , get_Ve(GLOBAL("n_traits"), GLOBAL("Vg"), h2))

    # Summary
    summary_phenome()
end


```
A quick start by assigning number of qtls,
vg and h2 can be optional to provide
```
function build_phenome(n_qtls       ::Union{Array{Int64, 1}, Int64};
                       vg           ::Union{Array{Float64 }, Float64}=1.0,
                       args...)

    # Handle different length of n_qtls and Vg
    n_traits = maximum([size(n_qtls)..., size(vg)...])

    # Instantiate QTL effects
    n_loci      = GLOBAL("n_loci")
    QTL_effects = spzeros(n_loci, n_traits)

    # Assign QTL effects
    if (n_traits > 1) & (length(n_qtls) == 1)
        # When n_qtls is a scaler, assign same number of QTLs for all traits
        idx_qtl = sample(1:n_loci, n_qtls, replace=false)
        for i in 1:n_traits
            QTL_effects[idx_qtl, i] = randn(n_qtls[i])
        end

    else
        # When n_qtls is a vector, assign different QTL locations for multiple traits
        for i in 1:n_traits
            idx_qtl = sample(1:n_loci, n_qtls, replace=false)
            QTL_effects[idx_qtl, i] = randn(n_qtls[i])
        end

    end

    # build_genome
    build_phenome(QTL_effects, vg=vg, args...)
end

```
Load dataframe to define effects
```
function build_phenome(dt            ::DataFrame; args...)
    QTL_effects = from_dt_to_eff(dt)
    build_phenome(QTL_effects; args...)
end

```
Load file to define effects
```
function build_phenome(filename     ::String; args...)
     build_phenome(
        CSV.read(filename, DataFrame);
        args...)
end




function summary_phenome()
    n_traits = GLOBAL("n_traits")
    n_qtls   = sum(GLOBAL("effects") .!=0, dims=1)
    Vg       = GLOBAL("Vg")
    Ve       = GLOBAL("Ve")

    if GLOBAL("n_traits") == 1
        h2 = Vg / (Vg + Ve)
    else
        h2 = diag(Vg ./ (Vg + Ve))
    end

    if !GLOBAL("silent")
        LOG("--------- Phenome Summary ---------")
        LOG("Number of Traits      : $n_traits")
        LOG("Heritability (h2)     : $h2")
        @info "" Genetic_Variance=Vg
        @info "" Residual_Variance=Ve
        LOG("Number of QTLs        : $n_qtls")
    end
end

function from_dt_to_eff(dt::DataFrame)
    columns = names(dt)
    idx_eff = [occursin("eff", s) for s in columns]
    return Matrix(dt[:, idx_eff])
end




