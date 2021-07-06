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
end


function build_genome(;# use ref species
                       species :: String="none",
                       # quick start
                       n_chr   :: Int64=-1,
                       n_marker:: Int64=-1,
                       args...)

    if (n_chr != -1 && n_marker != -1)
        # quick start
        n_row = n_chr * n_marker
        # chr
        chr  = repeat([1:n_chr;], inner=n_marker)
        # maf
        dist = Normal(0, .05)
        maf  = .5 .- abs.(rand(dist, n_row))
        # cM
        cM   = vcat([sort(uni_01(rand(dist, n_marker))) for _ in 1:n_chr]...)
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
        LOG("Chromosome Length (cM): $length_all")
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
        ref = DATA("genome_pig.csv")
        LOG("Tortereau,F. et al. (2012) A high density recombination map of the pig reveals a correlation between sex-specific recombination and GC content. BMC Genomics, 13, 586.")
        LOG("Reference Genome : Sscrofa 10.2")
        LOG("SNP Chip         : PorcineSNP60 BeadChip")

    elseif species == "cattle"
        ref = DATA("genome_cattle.csv")
        LOG("Arias,J.A. et al. (2009) A high density linkage map of the bovine genome. BMC Genetics, 10, 18.")
        LOG("Reference Genome : Btau 4.0")
        LOG("SNP Chip         : Affymetrix GeneChip Bovine Mapping 10K SNP kit")

    elseif species == "none"
        ref = missing

    else
        LOG("Assigned species not found, available options are: ['pig', 'cattle']", "error")

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





