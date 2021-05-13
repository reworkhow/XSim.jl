function build_genome(chromosome      ::Array{Int64,   1},
                      bp              ::Array{Int64,   1},
                      cM              ::Array{Float64, 1},
                      maf             ::Array{Float64, 1},
                      rate_mutation   ::Float64,
                      rate_error      ::Float64)

    CLEAR()
    SET("chromosome"   , chromosome)
    SET("bp"           , bp)
    SET("cM"           , cM)
    SET("maf"          , maf)
    SET("rate_mutation", rate_mutation)
    SET("rate_error"   , rate_error)
end


function build_genome(dt              ::DataFrame;
                      rate_mutation   ::Float64=0.0,
                      rate_error      ::Float64=0.0)

    if !all(in(names(dt)).(["chr", "bp", "cM"]))
        error("Missing required columns")
    end

    maf = all(in(names(dt)).(["maf"])) ? dt.maf : fill(0.5, GLOBAL("n_loci"))

    build_genome(dt.chr,
                 dt.bp,
                 dt.cM,
                 maf,
                 rate_mutation,
                 rate_error)

end

function build_genome(filename        ::String;
                      rate_mutation   ::Float64=0.0,
                      rate_error      ::Float64=0.0)

    build_genome(
        CSV.read(filename, DataFrame),
        rate_mutation=rate_mutation,
        rate_error   =rate_error)
end

function build_genome(;
                      species         ::String,
                      rate_mutation   ::Float64=0.0,
                      rate_error      ::Float64=0.0)

    if species == "Pig"
        build_genome(
            "../data/genome_pig.csv",
            rate_mutation=rate_mutation,
            rate_error   =rate_error)

    elseif species == "Cattle"
        build_genome(
            "../data/genome_cattle.csv",
            rate_mutation=rate_mutation,
            rate_error   =rate_error)
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





