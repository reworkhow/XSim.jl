mutable struct LocusInfo
    map_pos     ::Float64
    allele_freq ::Float64
    QTL         ::Bool
end

mutable struct ChromosomeInfo
    chrLength   ::Float64
    numLoci     ::Int64
    mapPos      ::Array{Float64,   1}
    loci        ::Array{LocusInfo, 1}
end

mutable struct GenomeInfo
    chr                 ::Array{ChromosomeInfo, 1}
    numChrom            ::Int64
    mutRate             ::Float64
    genotypeErrorRate   ::Float64
    qtl_index           ::Array{Int64,   1}
    qtl_effects         ::Array{Float64, 2} # second dim is for traits

    GenomeInfo() = new(Array{ChromosomeInfo}(undef, 0),
                       0, 0.0, 0.0,
                       Array{Int64, 1}(undef, 0),
                       Array{Float64, 2}(undef, 0, 0))
    GenomeInfo(
        chr                 ::Array{ChromosomeInfo, 1},
        numChrom            ::Int64,
        mutRate             ::Float64,
        genotypeErrorRate   ::Float64,
        qtl_index           ::Array{Int64, 1},
        qtl_effects         ::Array{Float64, 2}) =
    new(chr, numChrom, mutRate, genotypeErrorRate, qtl_index, qtl_effects)

end

mutable struct GLOBALS
    founders        ::Array{Animal, 1}
    G               ::GenomeInfo
    countChromosome ::Int64
    countId         ::Int64
    LRes            ::Array{Float64, 2} # Cholesky of residual covMat put
                                        # here instead of a scalar for varRes
                                        # as in single trait version
    GLOBALS() = new(Array{Animal}(undef, 0),
                    GenomeInfo(),
                    0,
                    0,
                    Array{Float64,2}(undef, 0, 0))
end


function clear_globals()
    global GLOBAL = GLOBALS()
end

export clear_globals

function setResidualVariance(LRes)
    common.LRes = LRes
end


function set_residual_var(LRes)
    GLOBAL.LRes = LRes
end

function set_num_loci(my::ChromosomeInfo, n::Int64)
    my.numLoci=n;
end

function set_num_chrom(my::GenomeInfo, n::Int64)
    my.numChrom=n; nothing
end

function get_num_loci(my::ChromosomeInfo)
    return my.numLoci
end

function get_num_alleles(my::LocusInfo)
    return length(my.allele_freq)
end

function get_num_chrom(my)
    return my.numChrom
end

function make_map_pos(my::GenomeInfo)
    for i in 1:my.numChrom
        mk_mappos_from_locus_info(my.chr[i])
    end
    mapPosDone = true
end

function mk_mappos_from_locus_info(my::ChromosomeInfo)
    my.mapPos.resize(my.chrLength)
    for i in 1:my.chrLength
        my.mapPos[i] = my.locus[i].map_pos
    end
end
