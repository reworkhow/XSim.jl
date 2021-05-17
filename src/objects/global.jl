mutable struct GB
    # Length = number of loci
    chromosome      ::Array{Int64,   1}
    bp              ::Array{Int64,   1}
    cM              ::Array{Float32, 1}
    maf             ::Array{Float32, 1}
    effects         ::Array{Float32, 2} # Second dim for traits
    effects_QTLs    ::Array{Float32, 2} # Second dim for traits
    is_QTLs         ::BitArray{      1} # Second dim for traits

    # Length = number of chromosome
    n_loci_chr      ::Array{Int64,   1}
    length_chr      ::Array{Float32, 1} # Unit is Morgan (100 cM)
    idx_chr         ::Array{Int64,   2} # chr by [start, end]

    # Scaler
    n_loci          ::Int64
    n_chr           ::Int64
    n_traits        ::Int64
    rate_mutation   ::Float64
    rate_error      ::Float64 # Genotyping error
    Vg              ::Union{Array{Float32, 2}, Array{Float64, 2}}

    # Counter
    founders        ::Array{Animal, 1}
    count_hap       ::Int64
    count_id        ::Int64

    # Constructor
    GB() = new(Array{Int64   }(undef, 0),
               Array{Int64   }(undef, 0),
               Array{Float32 }(undef, 0),
               Array{Float32 }(undef, 0),
               Array{Float32 }(undef, 0, 0),
               Array{Float32 }(undef, 0, 0),
               Array{BitArray}(undef, 0),
               Array{Int64   }(undef, 0),
               Array{Float32 }(undef, 0),
               Array{Int64   }(undef, 0, 0),
               0, 0, 0,
               0.0, 0.0,
               Array{Float32 }(undef, 0, 0),
               Array{Animal  }(undef, 0),
               1, 1)
end

function CLEAR()
    global gb = GB()
end

function SET(key   ::Any,
             value ::Any)

    setfield!(gb, Symbol(key), value)

    if key == "chromosome"
        SET("n_loci_chr"  , [count(==(c), value) for c in unique(value)])
        SET("n_loci"      , length(       value))
        SET("n_chr"       , length(unique(value)))
        # Index chromosome position
        idx_each_chr = [value .== c for c in 1:GLOBAL("n_chr")]
        SET("idx_chr"     , hcat(findfirst.(idx_each_chr), findlast.(idx_each_chr)))

    elseif key == "cM"
        chrs = GLOBAL("chromosome")
        SET("length_chr"  , [round(max(value[chrs.==c]...)/100, digits=2) for c in unique(chrs)])

    elseif key == "effects"
        SET("is_QTLs", sum(GLOBAL("effects"), dims=2)[:, 1] .!= 0)
        SET("effects_QTLs", GLOBAL("effects")[GLOBAL("is_QTLs"), :])
    end
end

function GLOBAL(option    ::String;
                chromosome::Int64=-1,
                locus     ::Int64=-1)

    if option == "n_loci" && chromosome != -1
        return getfield(gb, Symbol("n_loci_chr"))[chromosome]

    elseif option == "length_chr" && chromosome != -1
        return getfield(gb, Symbol("length_chr"))[chromosome]

    elseif chromosome != -1 && locus != -1
        return get_loci(chromosome, locus, option)

    elseif chromosome != -1
        return get_loci(chromosome, option)

    else
        return getfield(gb, Symbol(option))
    end
end


function add_count_ID!(;by::Int64=1)
    gb.count_id += by
end

function add_count_haplotype!(;by::Int64=1)
    gb.count_hap += by
end

function add_founder!(animal::Animal)
    push!(gb.founders, animal)
end

```Return info of specific loci```
function get_loci(chromosome::Int64, loci::Int64, option::String="bp")
    return get_loci(chromosome, option)[loci]
end

```Return info of  all loci on the chromosome```
function get_loci(chromosome::Int64, option::String="bp")
    idx_starts = GLOBAL("idx_chr")[chromosome, 1]
    idx_ends   = GLOBAL("idx_chr")[chromosome, 2]
    return GLOBAL(option)[idx_starts:idx_ends]
end

```Turn 1-D n-size vector, or a scaler to 2-D vector with dimension of n by 1```
function matrix(inputs::Any)
    mat = hcat(Diagonal([inputs])...)
    return convert(Array{Float32}, mat)
end

function handle_diagonal(inputs ::Union{Array{Float64}, Float64},
                         n_traits  ::Int64)

    # Cast variants of variances to a 2-D array
    if length(inputs) == 1 && !isa(inputs, Array)
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

function get_Vg(QTL_effects ::Array{Float32, 2},
                QTL_freq    ::Array{Float32, 1})

    # 2pq
    D = diagm(2 * QTL_freq .* (1 .- QTL_freq))

    # 2pq*alpha^2
    Vg = QTL_effects'D * QTL_effects

    return Vg
end

function scale_effects(QTL_effects ::Array{Float32, 2},
                       QTL_freq    ::Array{Float32, 1},
                       Vg_goal     ::Union{Array{Float32, 2}, Array{Float64, 2}})

    QTL_effects = convert(Array{Float32}, QTL_effects)
    QTL_freq    = convert(Array{Float32}, QTL_freq)

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

