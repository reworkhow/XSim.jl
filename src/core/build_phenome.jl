function build_phenome(QTL_effects  ::Union{Array{Float64}, SparseMatrixCSC},
                       Vg           ::Union{Array{Float64}, Float64}=1.0)

    SET("n_traits", size(QTL_effects)[2])
    SET("Vg"      , handle_diagonal(Vg,
                                    GLOBAL("n_traits")))

    SET("effects" , scale_effects(matrix(QTL_effects),
                                  GLOBAL("maf"),
                                  GLOBAL("Vg"),
                                  is_sparse=true))
    summary_phenome()
end

function build_phenome(n_qtls       ::Union{Array{Int64, 1}, Int64},
                       Vg           ::Union{Array{Float64 }, Float64}=1.0)

    n_traits = length(n_qtls)
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
    end

    # build_genome
    build_phenome(QTL_effects, Vg)
end


function summary_phenome()
    println("--------- Phenome Summary ---------")
    println("Number of traits      : ", GLOBAL("n_traits"))
    print(  "Genetic variance      : ")
    display(GLOBAL("Vg"))
    print(  "Number of QTLs        : ")
    display(sum(GLOBAL("effects") .!=0, dims=1))
end


