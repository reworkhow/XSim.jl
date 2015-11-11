module XSim

using Distributions

tempPos=Array(Float64,100000)
tempOri=Array(Int64,  100000)

include("genome/genome.jl")
include("cohort/cohort.jl")
include("output/output.jl")
include("global/global.jl")
include("mating/mating.jl")
include("deprecated.jl")

# initialize genome
function init(numChr::Int64,numLoci::Int64,chrLength::Float64,geneFreq::Array{Float64,1},
        mapPos::Array{Float64,1},qtl_marker::Array{Bool,1},qtl_effect::Array{Float64,1},mutRate::Float64,genotypeErrorRate=0.0,myCommon=common)

    #create genome
    locus_array = Array(LocusInfo,numLoci)
    QTL_index = Array(Int64,0)
    QTL_effect = Array(Float64,0)
    chr = Array(ChromosomeInfo,0)

    for j in 1:numChr
      for i in 1:numLoci
          if mapPos[i]>=chrLength
           error("Map posion is not on the chromosome (map position >= chromosome length)")
          end

          locus_array[i] = LocusInfo(mapPos[i],[geneFreq[i],1-geneFreq[i]],qtl_marker[i],qtl_effect[i])
          if qtl_marker[i]
            push!(QTL_index,numLoci*(j-1)+i) #make an array of QTL index for whole Genome
            push!(QTL_effect,qtl_effect[i])  #make an array of QTL effects for whole Genome
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
