module XSim

using Distributions
using DataFrames
using JWAS

tempPos=Array(Float64,100000)
tempOri=Array(Int64,  100000)

include("genome/genome.jl")
include("cohort/cohort.jl")
include("output/output.jl")
include("global/global.jl")
include("global/setParms.jl")
include("mating/mating.jl")
include("deprecated.jl")


# initialize genome
function init(numChr::Int64,numLoci::Int64,chrLength::Float64,geneFreq::Array{Float64,1},
        mapPos::Array{Float64,1},qtl_marker::Array{Bool,1},qtl_effect::Array{Float64,1},mutRate::Float64,genotypeErrorRate=0.0,myCommon=common) #assume same chromosomes

    #create genome
    locus_array = Array(LocusInfo,numLoci)
    QTL_index = Array(Int64,0)
    QTL_effect = Array(Float64,0)
    chr = Array(ChromosomeInfo,0)

    for j in 1:numChr
      for i in 1:numLoci
          if mapPos[i]>=chrLength
           error("Map position is not on the chromosome (map position >= chromosome length)")
          end

          locus_index    = numLoci*(j-1)+i
          locus_array[i] = LocusInfo(mapPos[i],[geneFreq[i],1-geneFreq[i]],qtl_marker[locus_index],qtl_effect[locus_index])
          if qtl_marker[locus_index]
            push!(QTL_index,locus_index) #make an array of QTL index for whole Genome
            push!(QTL_effect,locus_index)  #make an array of QTL effects for whole Genome
          end
      end
      chromosome = ChromosomeInfo(chrLength,numLoci,mapPos,locus_array)
      push!(chr,chromosome)
    end
    chromosome = ChromosomeInfo(chrLength,numLoci,mapPos,locus_array)
    chr = fill(chromosome,numChr)
    G = GenomeInfo(chr,numChr,mutRate,genotypeErrorRate,QTL_index,QTL_effect)

    # Init common
    myCommon.founders=Array(Animal,0)
    myCommon.G = G
    myCommon.countId = 1
    myCommon.countChromosome = 1; nothing
end


function init(numChr::Int64,
              numLoci::Array{Int64,1},
              chrLength::Array{Float64,1},
              geneFreq,mapPos,
              qtl_index,qtl_effect,
              mutRate::Float64=0.0,genotypeErrorRate=0.0,myCommon=common)

    QTL_index =[] #for whole genome
    QTL_effect=[] #for whole genome
    startlocus= 0 #locus index on whole genome

    chr = Array(ChromosomeInfo,0)
    for j in 1:numChr
      locus_array = Array(LocusInfo,numLoci[j])
      for i in 1:numLoci[j]
        if mapPos[j][i]>=chrLength[j]
          error("Map posion is not on the chromosome (map position >= chromosome length)")
        end
        pos = mapPos[j][i]
        locus_array[i] = LocusInfo(pos,[geneFreq[j][i],1-geneFreq[j][i]],false,0.0)
      end

      whichqtl=1
      for i in qtl_index[j]
        locus_array[i].QTL       =true
        locus_array[i].QTL_effect=qtl_effect[j][whichqtl]
        whichqtl += 1
      end

      chromosome = ChromosomeInfo(chrLength[j],numLoci[j],mapPos[j],locus_array)
      push!(chr,chromosome)

      QTL_index =vcat(QTL_index,qtl_index[j]+startlocus)
      QTL_effect=vcat(QTL_effect,qtl_effect[j])
      startlocus += numLoci[j]
    end
    G = GenomeInfo(chr,numChr,mutRate,genotypeErrorRate,QTL_index,QTL_effect)

    # Init common
    myCommon.founders=Array(Animal,0)
    myCommon.G = G
    myCommon.countId = 1
    myCommon.countChromosome = 1; nothing
end

"""
  build_genome(nChromosome,chromosome_length,nLoci,gene_frequency,map_position,qtl_index,qtl_effect,mutation_rate)

* nChromosome::Int64
  * number of chromosomes
* chromosome_length::Array{Float64,1}
  * length of each chromosome
* nLoci::Array{Int64,1}
  * number of loci for each chromosome
* gene_frequency::Array{Array{Float64,1},1}
  * gene frequency for each locus on each chromosome
* map_position::Array{Array{Float64,1},1}
  * map position for each locus on each chromosome
* qtl_index::Array{Array{Int64,1},1}
  * Index for each QTL on each chromosome
* qtl_effect::Array{Array{Float64,1},1}
  * Effect of QTL for each QTL on each chromosome
* mutation_rate::Float64
"""
function build_genome(nChromosome::Int64,
                      chromosome_length::Array{Float64,1},
                      nLoci::Array{Int64,1},
                      gene_frequency::Array{Array{Float64,1},1},
                      map_position::Array{Array{Float64,1},1},
                      qtl_index::Array{Array{Int64,1},1},
                      qtl_effect::Array{Array{Float64,1},1},
                      mutation_rate::Float64=0.0)
    init(nChromosome,nLoci,chromosome_length,gene_frequency,map_position,
         qtl_index,qtl_effect,mutation_rate)
end

function build_genome(nChromosome::Int64,
                      chromosome_length::Float64,
                      nLoci_each_chromosome::Int64,
                      qtl_each_chromosome::Int64,
                      mutation_rate::Float64=0.0)
    nLoci          = fill(nLoci_each_chromosome,nChromosome)
    chrLength      = fill(chromosome_length,nChromosome)
    map_position   = fill(collect(linspace(0.005,chromosome_length-0.005,nLoci_each_chromosome)),nChromosome)
    gene_frequency = fill(fill(0.5,nLoci_each_chromosome),nChromosome)
    qtl_index      = [sample(1:nLoci_each_chromosome,qtl_each_chromosome) for i in 1:nChromosome]
    qtl_effect     = fill(fill(0.0,qtl_each_chromosome),nChromosome)
    for i in 1:nChromosome
        qtl_effect[i] = randn(qtl_each_chromosome)
    end
    init(nChromosome,nLoci,chrLength,gene_frequency,map_position,
         qtl_index,qtl_effect,mutation_rate)
end


function init(numChr::Int64,numLoci::Int64,chrLength::Float64,geneFreq::Array{Float64,1},
        mapPos::Array{Float64,1},mutRate::Float64,genotypeErrorRate=0.0,myCommon=common)
    qtl_marker = fill(false,numLoci)
    qtl_effect = fill(0.0,numLoci)
    init(numChr,numLoci,chrLength,geneFreq, mapPos,qtl_marker,qtl_effect,mutRate,genotypeErrorRate,myCommon)
end

export build_genome
export sampleFounders,sampleRan,sampleSel,samplePed,concatCohorts,cohortSubset
export getOurGenotypes,getOurPhenVals,getOurGenVals
export outputPedigree,outputGenData,outputHapData,outputGenData,outputCatData
export getIDs,getPedigree
export startrPop #deprecated

end # module
