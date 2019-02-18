module XSim

using Distributions
using DataFrames
using CSV
using JWAS
using Printf
using LinearAlgebra

tempPos=Array{Float64}(undef,100000)
tempOri=Array{Int64}(undef,100000)
tempMut=Array{Float64}(undef,100000)

include("genome/genome.jl")
include("cohort/cohort.jl")
include("output/output.jl")
include("global/global.jl")
include("global/setParms.jl")
include("mating/mating.jl")
include("deprecated.jl")

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
                      mutation_rate::Float64=0.0,
                      genotypeErrorRate=0.0,myCommon=common)
    QTL_index  = Array{Int64}(undef, 0)  #for whole genome
    QTL_effect = Array{Float64}(undef, 0)#for whole genome
    chr        = Array{ChromosomeInfo}(undef, 0)#for whole genome

    startlocus= 0 #locus index on whole genome

    for j in 1:nChromosome
      locus_array = Array{LocusInfo}(undef, nLoci[j])
      for i in 1:nLoci[j]
        if map_position[j][i]>=chromosome_length[j]
          error("Map posion is not on the chromosome (map position >= chromosome length)")
        end
        pos = map_position[j][i]
        locus_array[i] = LocusInfo(pos,[gene_frequency[j][i],1-gene_frequency[j][i]],false,0.0)
      end

      whichqtl=1
      for i in qtl_index[j]
        locus_array[i].QTL       =true
        locus_array[i].QTL_effect=qtl_effect[j][whichqtl]
        whichqtl += 1
      end

      chromosome = ChromosomeInfo(chromosome_length[j],nLoci[j],map_position[j],locus_array)
      push!(chr,chromosome)

      QTL_index =vcat(QTL_index,qtl_index[j] .+ startlocus)
      QTL_effect=vcat(QTL_effect,qtl_effect[j])
      startlocus += nLoci[j]
    end
    G = GenomeInfo(chr,nChromosome,mutation_rate,genotypeErrorRate,QTL_index,QTL_effect)

    # Init common
    myCommon.founders=Array{Animal}(undef, 0)
    myCommon.G = G
    myCommon.countId = 1
    myCommon.countChromosome = 1; nothing
end

function build_genome(nChromosome::Int64,
                      chromosome_length::Float64,
                      nLoci_each_chromosome::Int64,
                      qtl_each_chromosome::Int64,
                      mutation_rate::Float64=0.0)

    nLoci          = fill(nLoci_each_chromosome,nChromosome)
    chrLength      = fill(chromosome_length,nChromosome)

    cstart         = chromosome_length/(nLoci_each_chromosome)/2
    cend           = chromosome_length-cstart
    map_position   = fill(collect(range(cstart,stop=cend,length=nLoci_each_chromosome)),nChromosome)
    gene_frequency = fill(fill(0.5,nLoci_each_chromosome),nChromosome)
    qtl_index      = [sample(1:nLoci_each_chromosome,qtl_each_chromosome,replace=false,ordered=true) for i in 1:nChromosome]
    qtl_effect     = fill(fill(0.0,qtl_each_chromosome),nChromosome)
    for i in 1:nChromosome
        qtl_effect[i] = randn(qtl_each_chromosome)
    end
    build_genome(nChromosome,chrLength,nLoci,gene_frequency,map_position,
         qtl_index,qtl_effect,mutation_rate)
end

function build_genome(nChromosome::Int64,
                      chromosome_length::Float64,
                      nLoci::Int64,
                      gene_frequency::Array{Float64,1},
                      map_position::Array{Float64,1},
                      mutation_rate::Float64=0.0,
                      qtl_index::Array{Int64,1}=Array{Int64,1}(undef, 0),
                      qtl_effect::Array{Float64,1}=Array{Float64,1}(undef, 0),
                      genotypeErrorRate=0.0,myCommon=common)

    nLoci             = fill(nLoci, nChromosome)
    chromosome_length = fill(chromosome_length,nChromosome)
    gene_frequency    = fill(gene_frequency,nChromosome)
    map_position      = fill(map_position,nChromosome)
    qtl_index         = fill(qtl_index,nChromosome)
    qtl_effect        = fill(qtl_effect,nChromosome)
    build_genome(nChromosome,chromosome_length,nLoci,gene_frequency,map_position,
                qtl_index,qtl_effect,mutation_rate)
end



export build_genome
export sampleFounders,sampleRan,sampleSel,samplePed,concatCohorts,cohortSubset,sampleBLUPSel,sampleDHOffspringFrom,sampleOneDHOffspringFrom
export getOurGenotypes,getOurPhenVals,getOurGenVals
export outputPedigree,outputGenData,outputHapData,outputGenData,outputCatData
export getIDs,getPedigree
export recode
export startrPop #deprecated

end # module
