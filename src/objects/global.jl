mutable struct GB
    # Length = number of loci
    chromosome::Array{Int64,1}
    bp::Array{Int64,1}
    cM::Array{Float64,1}
    maf::Array{Float64,1}
    effects::SparseMatrixCSC  # loci by traits
    effects_QTLs::SparseMatrixCSC  # qtls by traits
    is_QTLs::BitArray{1}

    # Length = number of chromosome
    n_loci_chr::Array{Int64,1}
    length_chr::Array{Float64,1} # Unit is Morgan (100 cM)
    idx_chr::Array{Int64,2} # chr by [start, end]

    # Scaler
    n_loci::Int64
    n_chr::Int64
    n_traits::Int64
    rate_mutation::Float64
    rate_error::Float64 # Genotyping error
    Vg::Array{Float64,2}
    Ve::Array{Float64,2}
    h2::Array{Float64,1}

    # Counter
    animals::Array{Animal,1}
    founders::Array{Animal,1}
    count_hap::Int64
    count_id::Int64

    # Others
    silent::Bool

    # Constructor
    GB() = new(Array{Int64}(undef, 0),
        Array{Int64}(undef, 0),
        Array{Float64}(undef, 0),
        Array{Float64}(undef, 0),
        spzeros(0),
        spzeros(0),
        Array{BitArray}(undef, 0),
        Array{Int64}(undef, 0),
        Array{Float64}(undef, 0),
        Array{Int64}(undef, 0, 0),
        0, 0, 0,
        0.0, 0.0,
        Array{Float64}(undef, 0, 0),
        Array{Float64}(undef, 0, 0),
        Array{Float64}(undef, 0),
        Array{Animal}(undef, 0),
        Array{Animal}(undef, 0),
        1, 1, false)
end

Base.show(io::IO, gb::GB) = ""


function CLEAR()
    global gb = GB()
    # LOG("XSim has been reset")
end

function SET(key::Any,
    value::Any)

    setfield!(gb, Symbol(key), value)

    if key == "chromosome"
        SET("n_loci_chr", [count(==(c), value) for c in unique(value)])
        SET("n_loci", length(value))
        SET("n_chr", length(unique(value)))
        # Index chromosome position
        idx_each_chr = [value .== c for c in 1:GLOBAL("n_chr")]
        SET("idx_chr", hcat(findfirst.(idx_each_chr), findlast.(idx_each_chr)))

    elseif key == "cM"
        chrs = GLOBAL("chromosome")
        SET("length_chr", [round(max(value[chrs.==c]...) / 100, digits=3) for c in unique(chrs)])

    elseif key == "effects"
        SET("is_QTLs", sum(GLOBAL("effects"), dims=2)[:, 1] .!= 0)
        SET("effects_QTLs", GLOBAL("effects")[GLOBAL("is_QTLs"), :])
    end
end

function GLOBAL(option::String="";
    chromosome::Int64=-1,
    locus::Int64=-1)

    if option == "n_loci" && chromosome != -1
        return getfield(gb, Symbol("n_loci_chr"))[chromosome]

    elseif option == "length_chr" && chromosome != -1
        return getfield(gb, Symbol("length_chr"))[chromosome]

    elseif chromosome != -1 && locus != -1
        return get_loci(chromosome, locus, option)

    elseif chromosome != -1
        return get_loci(chromosome, option)

    elseif option == "effects_QTLs"
        return Array(getfield(gb, Symbol(option)))

    elseif option == ""
        LOG("Available options are: ['chromosome', 'bp', 'cM', 'maf',
                               'effects', 'effects_QTLs', 'is_QTLs', 'animals',
                               'n_loci_chr', 'length_chr', 'idx_chr', 'n_loci',
                               'n_chr', 'n_traits', 'rate_mutation', 'rate_error',
                               'Vg', 'Ve', 'h2', 'error']", "error")

    else
        return getfield(gb, Symbol(option))
    end
end


function LOG(msg::String="",
    option::String="info";
    silent::Bool=GLOBAL("silent"))

    if !silent
        signiture = ""
        if option == "info"
            # @info "$signiture$msg"
            println("$signiture$msg")

        elseif option == "warn"
            @warn "$signiture$msg"

        elseif option == "error"
            error("$signiture$msg")

        end
    end
end

function SILENT(is_on::Bool=false)
    SET("silent", is_on)
    status = is_on ? "ON" : "OFF"
    # @info "The silent mode is $status"
end


```Return info of specific loci```
function get_loci(chromosome::Int64, loci::Int64, option::String="bp")
    return get_loci(chromosome, option)[loci]
end

```Return info of all loci on the chromosome```
function get_loci(chromosome::Int64, option::String="bp")
    idx_starts = GLOBAL("idx_chr")[chromosome, 1]
    idx_ends = GLOBAL("idx_chr")[chromosome, 2]
    return GLOBAL(option)[idx_starts:idx_ends]
end

function add_count_ID!(; by::Int64=1)
    gb.count_id += by
end

function add_count_haplotype!(; by::Int64=1)
    gb.count_hap += by
end

function add_animal!(animal::Animal)
    push!(gb.animals, animal)
end

function add_founders!(animal::Animal)
    push!(gb.founders, animal)
end

function GET_LINES(ids::Array)
    try
        ids = parse.(Int, ids)
    catch
        # Will be failed if ids are integers already
        nothing
    end
    ANIMALS = GLOBAL("animals")
    return Cohort([animal for animal in ANIMALS if animal.ID in ids])
end

function IS_EXIST(id::Int)
    return length(GET_LINES([id])) != 0
end


function scale_effects(QTL_effects::Union{Array{Float64,2},SparseMatrixCSC},
    QTL_freq::Array{Float64,1},
    Vg_goal::Array;
    is_sparse::Bool=false)

    # Compute Vg for input QTL_effects
    Vg_ori = get_Vg(QTL_effects, QTL_freq)

    # Decompose original variance
    Vg_ori_U = cholesky(Vg_ori).U'
    Vg_ori_Ui = inv(Vg_ori_U)

    # Decompose goal variance
    Vg_goal_U = cholesky(Vg_goal).U

    # (m by t) = (m x t) * (t x t)
    # Qs = Q * (v_goal / v_ori)
    QTL_effects_scaled = QTL_effects * Vg_ori_Ui'Vg_goal_U

    return is_sparse ? sparse(QTL_effects_scaled) : QTL_effects_scaled
end

function get_MAF(array::Array, is_haplotype=false)
    if is_haplotype
        array = hap_to_geno(array)
    end
    freq = sum(array, dims=1) / (2 * size(array, 1))
    maf = min.(freq, 1 .- freq)
    return round.(vcat(maf...), digits=3)
end

"""
Convert haplotypes (n by 2p) to genotypes matrix (n by p)
n: number of individuals
p: number of loci
"""
function hap_to_geno(array::Array)
    n_loci = Int(size(array)[2] / 2)
    for i in 1:n_loci
        i1 = (i * 2) - 1
        i2 = (i * 2)
        sub = array[:, i1:i2]
        array[:, i] = sum(sub, dims=2)
    end
    return array[:, 1:n_loci]
end

function uni_01(arr::Array)
    num = arr .- min(arr...)
    det = max(arr...) - min(arr...)
    return num / det
end


# function subset_dict(dict::Dict, subsets::Array)
#     k = collect(keys(dict))
#     v = collect(values(dict))

#     idx_subset = findall(in(subsets), k)
#     return Dict(k[i] => v[i] for i in idx_subset)
# end

# function make_silent(obj::Any)
#     args = convert(Dict, obj)
#     args["silent"] = true
#     return args
# end

