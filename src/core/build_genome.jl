"""
# Quick Start
Quick setup by assigning number of `markers` and `chromosomes`.

    build_genome(;n_loci ::Int64=-1,
                  n_chr    ::Int64=10,
                  species  ::String="none",
                  args...)

## Arguments
- `n_loci` : Number of simulated markers for each chromosome
- `n_chr`: Number of simulated chromosome with length of 100 centimorgan
- `species` : Infer genetic position (Morgan) by pre-load linkage maps, available species are: ["cattle", and "pig"]

## Examples
```jldoctest
julia> build_genome(n_chr  = 2,
                    n_loci = 10000)

[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info:
[ Info: Chromosome Length (cM): 200.0
[ Info: [100.0, 100.0]
[ Info:
[ Info: Number of Loci        : 20000
[ Info: [10000, 10000]
[ Info:
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info:
```


# Define genome by a file or a DataFrame
Define genome by providing a formatted dataframe or a path to the file.

    build_genome(dt      ::DataFrame;
                 species ::String="none",
                 args...)

    build_genome(filename::String;
                 species ::String="none",
                 args...)

## Arguments
- `dt` : A `DataFrame` with required columns of `chr` and `cM`, or `chr` and `bp` if `species` is provided for the inference.
- `filename` : A filepath to the file containing genome information.
- `species` : Adjust genetic position (Morgan) by pre-load linkage maps, available species are: ["cattle", and "pig"]

## Example of the `DataFrame`
```
4×7 DataFrame
 Row │ id      chr    bp       cM       MAF      eff_1    eff_2
     │ String  Int64  Int64    Float64  Float64  Float64  Float64
─────┼────────────────────────────────────────────────────────────
   1 │ snp_1       1  1818249     50.8      0.5      0.1      0.0
   2 │ snp_2       1  6557697     80.3      0.5      0.0      0.0
   3 │ snp_3       2  2298800     39.2      0.5      0.2      0.0
   4 │ snp_4       2  5015698     66.3      0.5      0.0      0.5
```

## Examples
By a filepath
```jldoctest
julia> build_genome("path/map.csv";
                    rate_mutation=0.005, rate_error=0.01)
```

or a dataframe
```jldoctest
julia> using DataFrames
julia> data = CSV.read("path/map.csv", DataFrame)
julia> build_genome(data;
                    rate_mutation=0.005, rate_error=0.01)

[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info:
[ Info: Chromosome Length (cM): 146.6
[ Info: [80.3, 66.3]
[ Info:
[ Info: Number of Loci        : 4
[ Info: [2, 2]
[ Info:
[ Info: Genotyping Error      : 0.01
[ Info: Mutation Rate         : 0.005
[ Info:
```

Use cattle genome as reference to infer the genetic positions
```jldoctest
julia> build_genome("path/map.csv"; species="cattle")

[ Info: Arias,J.A. et al. (2009) A high density linkage map of the bovine genome. BMC Genetics, 10, 18.
[ Info: Reference Genome : Btau 4.0
[ Info: SNP Chip         : Affymetrix GeneChip Bovine Mapping 10K SNP kit

┌ Warning: The provided genetic distances will be replaced with ones infered from preloaded linkage maps
└ @ XSim ~/Dropbox/projects/XSim/src/objects/global.jl:118
[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info:
[ Info: Chromosome Length (cM): 16.8
[ Info: [15.1, 1.7]
[ Info:
[ Info: Number of Loci        : 4
[ Info: [2, 2]
[ Info:
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info:

```

# Explict Definition
Define genome by providing genetic information of each loci explicitly.

    build_genome(chromosome    ::Array{Int64,   1},
                 bp            ::Array{Int64,   1},
                 cM            ::Array{Float64, 1},
                 maf           ::Array{Float64, 1};
                 rate_mutation ::Float64=0.0,
                 rate_error    ::Float64=0.0)

## Arguments
- `chromosome` : Chromosome codes
- `bp` : Physical positions
- `cM` : Genetic positions
- `maf` : Minor allele frequencies
- `rate_mutation` : Mutation rate
- `rate_error` : Error rate of genotyping

## Examples
```jldoctest
julia> ch  = [1,    1,     2,    2,    2]
julia> bp  = [130,  205,   186,  503,  780]
julia> cM  = [85.7, 149.1, 37.4, 83.6, 134.3]
julia> maf = [0.5,  0.5,   0.5,  0.5,  0.5]
julia> build_genome(ch, bp, cM, maf)

[ Info: --------- Genome Summary ---------
[ Info: Number of Chromosome  : 2
[ Info:
[ Info: Chromosome Length (cM): 283.4
[ Info: [149.1, 134.3]
[ Info:
[ Info: Number of Loci        : 5
[ Info: [2, 3]
[ Info:
[ Info: Genotyping Error      : 0.0
[ Info: Mutation Rate         : 0.0
[ Info:
```
"""
function build_genome(chromosome      ::Array{Int64,   1},
                      bp              ::Array{Int64,   1},
                      cM              ::Array{Float64, 1},
                      maf             ::Array{Float64, 1};
                      rate_mutation   ::Float64=0.0,
                      rate_error      ::Float64=0.0)

    is_silent = GLOBAL("silent")
    CLEAR()
    SET("chromosome"   , chromosome)
    SET("bp"           , bp)
    SET("cM"           , cM)
    SET("maf"          , maf)
    SET("rate_mutation", rate_mutation)
    SET("rate_error"   , rate_error)
    SET("silent"       , is_silent)

    summary_genome()
    build_phenome(1) # in some cases, users don't need simulated phenotypes
end

function build_genome(;# use ref species
                       species :: String="none",
                       # quick start
                       n_loci  :: Int64=-1,
                       n_chr   :: Int64=10,
                       args...)

    if (n_loci != -1)
        # quick start
        n_row = n_chr * n_loci
        # chr
        chr  = repeat([1:n_chr;], inner=n_loci)
        # maf
        dist = Normal(0, .05)
        maf  = .5 .- abs.(rand(dist, n_row))
        # cM
        cM   = vcat([Base.sort(uni_01(rand(dist, n_loci))) .* 100 for _ in 1:n_chr]...)
        # bp
        bp  = fill(0, n_row)
        # build genome
        build_genome(chr,
                     bp,
                     cM,
                     maf,
                     args...)

    elseif (species != "none")
        # get ref
        build_genome(load_ref(species);
                     args...)

    else
        LOG("The usage is not valid", "error")
    end

end

function build_genome(dt      :: DataFrame;
                      species :: String="none",
                      args...)

    # load reference if provdied
    ref = load_ref(species)

    # check columns
    columns = names(dt)
    has_chr, has_bp, has_cM, has_maf = in(columns).(["chr", "bp", "cM", "maf"])

    # infer cM
    if all([has_chr, has_cM])
        if ismissing(ref)
            # pass
        else
            if has_bp
                # want to recalculate cM based on bp and ref
                add_cM_by_ref!(dt, ref)
                LOG("The provided genetic distances will be replaced with ones infered from preloaded linkage maps", "warn")
            else
                # no bp but provide ref and cM, don't know what to do
                LOG("Missing required column 'bp'", "error")
            end
        end

    elseif all([has_chr, !has_cM, has_bp])
        if ismissing(ref)
            # infer cM by bp linearly
            add_cM_by_bp!(dt)
        else
            # infer cM by reference and bp
            add_cM_by_ref!(dt, ref)
        end

    else
        LOG("Missing required columns", "error")
    end

    # bp and MAF
    dt.maf = has_maf ? dt.maf : fill(0.5, nrow(dt))
    dt.bp  = has_bp  ? dt.bp  : fill(0,   nrow(dt))

    # build genome
    build_genome(dt.chr,
                 dt.bp,
                 dt.cM,
                 dt.maf;
                 args...)
end

function build_genome(filename::String; args...)
    build_genome(
        CSV.read(filename, DataFrame);
        args...)
end

function summary_genome()
    n_chr      = GLOBAL("n_chr")
    length_cM  = round.(GLOBAL("length_chr") .* 100, digits=1)
    length_all = sum(length_cM)
    n_loci     = GLOBAL("n_loci")
    n_loci_chr = GLOBAL("n_loci_chr")
    rate_error = GLOBAL("rate_error")
    rate_mut   = GLOBAL("rate_mutation")

    if !GLOBAL("silent")
        LOG("--------- Genome Summary ---------")
        LOG("Number of Chromosome  : $n_chr")
        LOG()
        LOG("Chromosome Length (cM):")
        LOG("$length_cM")
        LOG()
        LOG("Number of Loci        : $n_loci")
        LOG("$n_loci_chr")
        LOG()
        LOG("Genotyping Error      : $rate_error")
        LOG("Mutation Rate         : $rate_mut")
        LOG()
    end
end

# PRIVATE --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
function load_ref(species::String)
    if species == "pig"
        ref = DATA("genome_pig", header=true)
        LOG("Tortereau,F. et al. (2012) A high density recombination map of the pig reveals a correlation between sex-specific recombination and GC content. BMC Genomics, 13, 586.")
        LOG("Reference Genome : Sscrofa 10.2")
        LOG("SNP Chip         : PorcineSNP60 BeadChip")

    elseif species == "cattle"
        ref = DATA("genome_cattle", header=true)
        LOG("Arias,J.A. et al. (2009) A high density linkage map of the bovine genome. BMC Genetics, 10, 18.")
        LOG("Reference Genome : Btau 4.0")
        LOG("SNP Chip         : Affymetrix GeneChip Bovine Mapping 10K SNP kit")

    elseif species == "maize"
        ref = DATA("genome_maize", header=true)
        LOG("Portwood,J.L. et al. (2019) MaizeGDB 2018: the maize multi-genome genetics and genomics database. Nucleic Acids Research 47, D1146–D1154")
        LOG("Reference Genome : B73 v3 and v4")
        LOG("SNP Chip         : IBM MaizeSNP50")

    elseif species == "rice"
        ref = DATA("genome_rice", header=true)
        LOG("Kurata, N., and Yamazaki, Y. (2006). Oryzabase. An integrated biological and genome information database for rice. Plant Physiol 140, 12–17.")

    elseif species == "none"
        ref = missing

    else
        LOG("Assigned species not found, available options are: ['pig', 'cattle', 'maize', 'rice']", "error")

    end

    return ref
end

function add_cM_by_bp!(dt::DataFrame)
    dt[:, "cM"] .= 0.0
    for chr in unique(dt.chr)
        dt[dt.chr .== chr, "cM"] .= uni_01(dt[dt.chr .== chr, "bp"]) .* 100
    end
end

function add_cM_by_ref!(dt::DataFrame, ref::DataFrame)
    dt[:, "cM"] .= 0.0
    for chr in unique(dt.chr)
        # fetch info from user and reference
        bp_user = dt[dt.chr .== chr, "bp"]
        bp_ref  = ref[ref.chr .== chr, "bp"]
        cM_ref  = ref[ref.chr .== chr, "cM"]
        # find the nearest marker from reference
        idx_new_cM = [argmin(abs.(bp_ref .- bp_user[i])) for i in 1:length(bp_user)]
        dt[dt.chr .== chr, "cM"] = cM_ref[idx_new_cM]
    end
end





# """
#   build_genome(n_chr,len_chr,n_loci,gene_frequency,map_position,qtl_index,qtl_effect,mutation_rate)

# * n_chr::Int64
#   * number of chromosomes
# * len_chr::Array{Float64,1}
#   * length of each chromosome
# * n_loci::Array{Int64,1}
#   * number of loci for each chromosome
# * gene_frequency::Array{Array{Float64,1},1}
#   * gene frequency for each locus on each chromosome
# * map_position::Array{Array{Float64,1},1}
#   * map position for each locus on each chromosome
# * qtl_index::Array{Array{Int64,1},1}
#   * Index for each QTL on each chromosome
# * qtl_effect::Array{Array{Float64,1},1}
#   * Effect of QTL for each QTL on each chromosome
# * mutation_rate::Float64
# """
# function build_genome(n_chr          ::Int64,
#                       len_chr        ::Array{Float64, 1},
#                       n_loci         ::Array{Int64,   1},

#                       gene_frequency ::Array{Array{Float64, 1}, 1},
#                       map_position   ::Array{Array{Float64, 1}, 1},

#                       qtl_index      ::Array{Array{Int64,   1}, 1}, # store which loci is qtl
#                       qtl_effect     ::Array{Array{Float64, 2}, 1}, #
#                       n_trait        ::Int64=1,
#                       G0=[],
#                       mutation_rate  ::Float64=0.0,
#                       rate_error=0.0)

#     numQTLOnChr       = Array{Int64}(undef, 0)
#     QTL_index         = Array{Int64}(undef, 0)
#     QTLEffectsMat     = Array{Float64,2}(undef, 0, n_trait)
#     gene_frequencyQTL = Array{Array{Float64,1},1}(undef, 0)
#     chrs              = Array{ChromosomeInfo}(undef, 0)

#     startlocus = 0 #locus index on whole genome

#     for j in 1:n_chr
#         locus_array = Array{LocusInfo}(undef, n_loci[j])

#         # define loci in each chromosome
#         for i in 1:n_loci[j]
#             pos = map_position[j][i]
#             locus_array[i] = LocusInfo(pos, gene_frequency[j][i], false)
#         end

#         for i in qtl_index[j]
#             locus_array[i].is_qtl = true
#         end

#         chromosome = ChromosomeInfo(len_chr[j],n_loci[j],map_position[j],locus_array)
#         push!(chrs, chromosome)


#         if size(G0, 1) > 0
#             nQTLPerChr    = length(qtl_index[j])
#             numQTLOnChr   = push!(numQTLOnChr,nQTLPerChr)
#             push!(gene_frequencyQTL,(gene_frequency[j])[qtl_index[j]])
#         else
#             QTLEffectsMat = vcat(QTLEffectsMat, qtl_effect[j])
#         end
#         QTL_index = vcat(QTL_index, qtl_index[j] .+ startlocus)
#         startlocus += n_loci[j]
#     end

#     if size(G0, 1) > 0
#         qtl_effect, QTLEffectsMat = transformEffects(numQTLOnChr, qtl_effect, gene_frequencyQTL, G0)
#         println("G after transformation = ", 0.5 * QTLEffectsMat'QTLEffectsMat)
#     end
#     G = GenomeInfo(chrs, mutation_rate, rate_error, QTL_index, QTLEffectsMat)

#     # Init common
#     GLOBAL.founders=Array{Animal}(undef, 0)
#     GLOBAL.G = G
#     GLOBAL.countId = 1
#     GLOBAL.countChromosome = 1;
# end


# function build_genome(n_chr::Int64,
#                       len_chr::Float64,
#                       n_loci_each_chromosome::Int64,
#                       qtl_each_chromosome::Int64;
#                       n_trait::Int64=1,
#                       G0 = [],
#                       mutation_rate::Float64=0.0)
#     n_loci         = fill(n_loci_each_chromosome,n_chr)
#     chrLength      = fill(len_chr,n_chr)

#     cstart         = len_chr/(n_loci_each_chromosome)/2
#     cend           = len_chr-cstart
#     map_position   = fill(collect(range(cstart,stop=cend,length=n_loci_each_chromosome)),n_chr)
#     gene_frequency = fill(fill(0.5,n_loci_each_chromosome),n_chr)
#     qtl_index      = [sample(1:n_loci_each_chromosome, qtl_each_chromosome, replace=false, ordered=true) for i in 1:n_chr]
#     qtl_effects    = fill(fill(0.0,qtl_each_chromosome,n_trait),n_chr)
#     for i in 1:n_chr
#         qtl_effects[i] = randn(qtl_each_chromosome,n_trait)
#     end

#     build_genome(n_chr,chrLength,n_loci,gene_frequency,map_position,
#          qtl_index,qtl_effects,n_trait,G0,mutation_rate)
# end

# function build_genome(n_chr::Int64,
#                       len_chr::Float64,
#                       n_loci::Int64,
#                       gene_frequency::Array{Float64,1},
#                       map_position::Array{Float64,1},
#                       mutation_rate::Float64=0.0,
#                       qtl_index::Array{Int64,1}=Array{Int64,1}(undef, 0),
#                       qtl_effects::Array{Float64,2}=Array{Float64,2}(undef, 0, 0),
#                       G0=[],
#                       rate_error=0.0)

#     n_trait           = size(qtl_effects,2)
#     n_loci            = fill(n_loci, n_chr)
#     len_chr           = fill(len_chr,n_chr)
#     gene_frequencyQTL = fill(gene_frequency[qtl_index],n_chr)
#     gene_frequency    = fill(gene_frequency,n_chr)
#     map_position      = fill(map_position,n_chr)
#     qtl_index         = fill(qtl_index,n_chr)
#     qtl_effects       = fill(qtl_effects,n_chr)

#     build_genome(n_chr, len_chr, n_loci,
#                 gene_frequency, map_position,
#                 qtl_index, qtl_effects, n_trait, G0,
#                 mutation_rate, rate_error)
# end





