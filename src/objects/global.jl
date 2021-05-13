mutable struct Global
    # Length = number of loci
    chromosome      ::Array{Int64,   1}
    bp              ::Array{Int64,   1}
    cM              ::Array{Float64, 1}
    maf             ::Array{Float64, 1}
    effects         ::Array{Float64, 2} # Second dim for traits

    # Length = number of chromosome
    n_loci_chr      ::Array{Int64,   1}
    length_chr      ::Array{Float64, 1} # It's cM

    # Singular
    n_loci          ::Int64
    n_chr           ::Int64
    n_traits        ::Int64
    rate_mutation   ::Float64
    rate_error      ::Float64 # Genotyping error
    Vg              ::Array{Float64, 2}

    # Counter
    founders        ::Array{Animal, 1}
    count_hap       ::Int64
    count_id        ::Int64

    # Constructor
    Global() = new(Array{Int64  }(undef, 0),
                   Array{Int64  }(undef, 0),
                   Array{Float64}(undef, 0),
                   Array{Float64}(undef, 0),
                   Array{Float64}(undef, 0, 0),
                   Array{Int64  }(undef, 0),
                   Array{Float64}(undef, 0),
                   0, 0,
                   0.0, 0.0,
                   Array{Float64}(undef, 0, 0),
                   Array{Animal }(undef, 0),
                   1, 1)
end


function CLEAR()
    global Global = GLobal()
end

function SET(key   ::Any,
             value ::Any)

    setfield!(Global, Symbol(key), value)

    if key == "chromosome"
        SET("n_loci_chr", [count(==(c), value) for c in unique(value)])
        SET("n_loc"     , length(       value))
        SET("n_chr"     , length(unique(value)))

    elseif key == "cM"
        chrs = GLOBAL("chromosome")
        SET("length_chr", [max.(value[chrs.==c]...) for c in unique(chrs)])

    end
end


function GLOBAL(option    ::String;
                chromosome::Int64=-1,
                locus     ::Int64=-1)

    if option == "n_loci" & chromosome != -1
        return getfield(Global, Symbol("n_loci_chr"))[chromosome]

    elseif option == "length_chr" & chromosome != -1
        return getfield(Global, Symbol("length_chr"))[chromosome]

    elseif chromosome != -1 & locus != -1
        return get_loci(chromosome, locus, option)

    elseif chromosome != -1
        return get_loci(chromosome, option)

    else
        return getfield(Global, Symbol(option))
    end
end


function add_count_ID!(;by::Int64=1)
    Global.count_id += by
end

function add_count_haplotype!(;by::Int64=1)
    Global.count_hap += by
end

function add_founder!(animal::Animal)
    push!(Global.founders, animal)
end

```Return info of specific loci```
function get_loci(chromosome::Int64, loci::Int64, option::String="bp")
    return get_loci(chromosome, option)[loci]
end

```Return info of  all loci on the chromosome```
function get_loci(chromosome::Int64, option::String="bp")
    idx_starts = findfirst(Global.chromosome .== chromosome)
    idx_ends   = findlast(Global.chromosome .== chromosome)

    return GLOBAL(option)[idx_starts:idx_ends]
end

```Turn 1-D n-size vector, or a scaler to 2-D vector with dimension of n by 1```
function matrix(inputs::Any)
    return hcat(Diagonal([inputs])...)
end

function handle_diagonal(inputs ::Union{Array{Float64}, Float64},
                         n_traits  ::Int64)

    # Cast variants of variances to a 2-D array
    if length(inputs) == 1
        # When variances is a scaler, assign it as the diagonal of variances
        inputs = diagm(fill(inputs, n_traits))
    else
        inputs = matrix(inputs)
        if size(inputs)[2] == 1
            # When variances is a vector, assign it as the diagonal of variances
            inputs = diagm(inputs[:, 1])
        end
    end

    if size(inputs)[2] != n_traits
        error("Dimensions don't match between n_traits and variances/h2")
    end

    return inputs
end

function get_Vg(QTL_effects ::Array{Float64, 2},
                QTL_freq    ::Array{Float64, 1})

    # 2pq
    D = diagm(2 * QTL_freq .* (1 .- QTL_freq))

    # 2pq*alpha^2
    Vg = QTL_effects'D * QTL_effects

    return Vg
end

function scale_effects(QTL_effects ::Array{Float64, 2},
                       QTL_freq    ::Array{Float64, 1},
                       Vg_goal     ::Array{Float64, 2})

    # Compute Vg for input QTL_effects
    Vg_ori    = get_Vg(QTL_effects, QTL_freq)

    # Decompose original variance
    Vg_ori_U  = cholesky(Vg_ori).U'
    Vg_ori_Ui = inv(Vg_ori_U)

    # Decompose goal variance
    Vg_goal_U = cholesky(Vg_goal).U

    # m by t = m by t * t by t
    QTL_effects_scaled = QTL_effects * Vg_ori_Ui'Vg_goal_U

    return QTL_effects_scaled
end

