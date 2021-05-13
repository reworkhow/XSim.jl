function build_phenome(QTL_effects  ::Array{Float64},
                       Vg           ::Union{Array{Float64}, Float64}=matrix(1))

    SET("n_traits", size(QTL_effects)(2))
    SET("Vg"      , handle_diagonal(Vg,
                                     GLOBAL("n_traits")))
    SET("effects" , scale_effects(matrix(QTL_effects),
                                  GLOBAL("maf"),
                                  GLOBAL("Vg")))
end

function build_phenome(n_traits     ::Int64,
                       n_qtls       ::Union{Array{Int64, 1}, Int64},
                       Vg           ::Union{Array{Float64 }, Float64}=matrix(1))

    # When n_qtls is a scaler, assign same number of QTLs for all traits
    if n_traits > 1 & length(n_qtls) == 1
        n_qtls = fill(n_qtls, n_triats)
    end

    # Simulate QTL effects
    n_loci = GLOBAL("n_loci")
    QTL_effects = zeros(n_loci, n_traits)
    for i in 1:n_traits
        idx_qtl = sample(1:n_loci, n_qtls[i], replace=false)
        QTL_effects[idx_qtl, i] = randn(n_qtls[i])
    end

    # Build
    build_phenome(QTL_effects, Vg)
end


