function build_phenome(QTL_effects  ::Union{Array{Float64}, SparseMatrixCSC},
                       Vg           ::Union{Array{Float64}, Float64}=1.0)

    SET("n_traits", size(QTL_effects)[2])

    Vg = handle_diagonal(Vg, GLOBAL("n_traits"))
    Vg = Symmetric(Vg)
    SET("Vg"      , Array(Vg))

    SET("effects" , scale_effects(matrix(QTL_effects),
                                  GLOBAL("maf"),
                                  GLOBAL("Vg"),
                                  is_sparse=true))
    summary_phenome()
end

function build_phenome(n_qtls       ::Union{Array{Int64, 1}, Int64},
                       Vg           ::Union{Array{Float64 }, Float64}=1.0)

    # Could be n_traits = length(n_qtls)
    n_traits = maximum([size(n_qtls)..., size(Vg)...])

    # When n_qtls is a scaler, assign same number of QTLs for all traits
    if (n_traits > 1) & (length(n_qtls) == 1)
        n_qtls = fill(n_qtls, n_traits)
    end

    # Simulate QTL effects
    n_loci      = GLOBAL("n_loci")
    QTL_effects = spzeros(n_loci, n_traits)
    for i in 1:n_traits
        idx_qtl                 = sample(1:n_loci, n_qtls[i], replace=false)
        QTL_effects[idx_qtl, i] = randn(n_qtls[i])
        println(QTL_effects[idx_qtl, i])
    end

    # build_genome
    build_phenome(QTL_effects, Vg)
end


function summary_phenome()
    n_traits = GLOBAL("n_traits")
    n_qtls   = sum(GLOBAL("effects") .!=0, dims=1)
    Vg       = GLOBAL("Vg")

    if !GLOBAL("silent")
        LOG("--------- Phenome Summary ---------")
        LOG("Number of Traits      : $n_traits")
        @info "" Genetic_Variance=Vg
        LOG("Number of QTLs        : $n_qtls")
    end
end




