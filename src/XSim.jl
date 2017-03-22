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


function init(numChr::Int64,numLoci::Array{Int64,1},chrLength::Array{Float64,1},
              geneFreq,mapPos,qtl_marker, #Array{Array{T,N},1}
              qtl_effect,mutRate::Float64,
              genotypeErrorRate=0.0,myCommon=common)

    #create genome
    QTL_index = Array(Int64,0)    #for whole genome
    QTL_effect = Array(Float64,0) #for whole genome
    chr = Array(ChromosomeInfo,0)

    locus_array = Array(LocusInfo,0)
    whichlocus=0
    for j in 1:numChr
        locus_array = Array(LocusInfo,numLoci[j])
        for i in 1:numLoci[j]
          if mapPos[j][i]>=chrLength[j]
            error("Map posion is not on the chromosome (map position >= chromosome length)")
          end
          pos = mapPos[j][i]
          locus_array[i] = LocusInfo(pos,[geneFreq[j][i],1-geneFreq[j][i]],
            qtl_marker[j][i],qtl_effect[j][i])
          if qtl_marker[j][i]
            push!(QTL_index,whichlocus+i)       #make an array of QTL index for whole Genome
            push!(QTL_effect,qtl_effect[j][i])  #make an array of QTL effects for whole Genome
          end
        end
        whichlocus = whichlocus + numLoci[j]
        chromosome = ChromosomeInfo(chrLength[j],numLoci[j],mapPos[j],locus_array)
        push!(chr,chromosome)
    end
    G = GenomeInfo(chr,numChr,mutRate,genotypeErrorRate,QTL_index,QTL_effect)

    # Init common
    myCommon.founders=Array(Animal,0)
    myCommon.G = G
    myCommon.countId = 1
    myCommon.countChromosome = 1; nothing
end

function init(numChr::Int64,numLoci::Int64,chrLength::Float64,geneFreq::Array{Float64,1},
        mapPos::Array{Float64,1},mutRate::Float64,genotypeErrorRate=0.0,myCommon=common)
    qtl_marker = fill(false,numLoci)
    qtl_effect = fill(0.0,numLoci)
    init(numChr,numLoci,chrLength,geneFreq, mapPos,qtl_marker,qtl_effect,mutRate,genotypeErrorRate,myCommon)
end

export sampleFounders,sampleRan,sampleSel,samplePed,concatCohorts,cohortSubset
export getOurGenotypes,getOurPhenVals,getOurGenVals
export outputPedigree,outputGenData,outputHapData,outputGenData,outputCatData
export getIDs,getPedigree
export startrPop #deprecated

end # module
