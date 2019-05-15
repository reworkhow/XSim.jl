Types and Methods for Information on the Genome

###Types
mutable struct LocusInfo
    map_pos::Float64
    allele_freq::Array
    QTL::Bool
    #QTL_effect::Float64 deprecated, not used
end

mutable struct ChromosomeInfo
    chrLength::Float64
    numLoci::Int64
    mapPos::Array{Float64,1}
    loci::Array{LocusInfo,1}
end

mutable struct GenomeInfo
    chr::Array{ChromosomeInfo,1}
    numChrom::Int64
    mutRate::Float64
    genotypeErrorRate::Float64
    qtl_index::Array{Int64,1}
    qtl_effects::Array{Float64,2} # second dim is for traits
end


###Methods

##set genome information: number of loci,chromosomes
function set_num_loci(my::ChromosomeInfo,n::Int64)
    my.numLoci=n; nothing
end

function set_num_chrom(my::GenomeInfo,n::Int64)
    my.numChrom=n; nothing
end

##get genome information: number of loci, chromosomes
function get_num_loci(my::ChromosomeInfo)
    return my.numLoci
end

function get_num_alleles(my::LocusInfo)
    return length(my.allele_freq)
end

function get_num_chrom(my)
    return my.numChrom
end

##set map positions information
function make_map_pos(my::GenomeInfo)
    for i in 1:my.numChrom
        mk_mappos_from_locus_info(my.chr[i])
    end
    mapPosDone = true
end

function mk_mappos_from_locus_info(my::ChromosomeInfo)
    my.mapPos.resize(my.chrLength)
    for i in 1:my.chrLength
        my.mapPos[i]=my.locus[i].map_pos
    end
end
