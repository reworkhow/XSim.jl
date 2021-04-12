module XSim

using Distributions
using DataFrames
using CSV
using JWAS
using Printf
using LinearAlgebra

tempPos = Array{Float64}(undef, 100000)
tempOri = Array{Int64}(undef, 100000)
tempMut = Array{Float64}(undef, 100000)

"""
base type for genotype and Haplotype storage.
Shouldn't be exported but needs to be defined. Used throughout the included src files.
Originally an Int64
"""
const AlleleIndexType = Int8

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
                      chromosome_length::Array{Float64, 1},
                      nLoci::Array{Int64, 1},
                      gene_frequency::Array{Array{Float64, 1}, 1},
                      map_position::Array{Array{Float64,1 }, 1},
                      qtl_index::Array{Array{Int64, 1}, 1},
                      qtl_effect::Array{Array{Float64, 2}, 1},
                      nTraits::Int64=1,
                      G0=[],
                      mutation_rate::Float64=0.0,
                      genotypeErrorRate=0.0,
                      myCommon=common)
                      println("G0 = ", G0)

    numQTLOnChr       = Array{Int64}(undef, 0)  #for whole genome
    QTL_index         = Array{Int64}(undef, 0)  #for whole genome
    QTLEffectsMat     = Array{Float64,2}(undef, 0, nTraits) #for whole genome
    gene_frequencyQTL = Array{Array{Float64,1},1}(undef, 0)  #for whole genome
    chr               = Array{ChromosomeInfo}(undef, 0) #for whole genome

    startlocus = 0 #locus index on whole genome

    for j in 1:nChromosome
      locus_array = Array{LocusInfo}(undef, nLoci[j])
      for i in 1:nLoci[j]
        if map_position[j][i] >= chromosome_length[j]
          error("Map posion is not on the chromosome (map position >= chromosome length)")
        end
        pos = map_position[j][i]
        locus_array[i] = LocusInfo(pos, [gene_frequency[j][i], 1 - gene_frequency[j][i]], false)
      end

      whichqtl = 1
      for i in qtl_index[j]
        locus_array[i].QTL = true
        #locus_array[i].QTL_effect=qtl_effect[j][whichqtl] # we think this is not used; G.qtl_effects[][] is used instead
        whichqtl += 1
      end

      chromosome = ChromosomeInfo(chromosome_length[j],nLoci[j],map_position[j],locus_array)
      push!(chr,chromosome)

      if size(G0,1) > 0
        nQTLPerChr    = length(qtl_index[j])
        numQTLOnChr   = push!(numQTLOnChr,nQTLPerChr)
        push!(gene_frequencyQTL,(gene_frequency[j])[qtl_index[j]])
      else
          QTLEffectsMat = vcat(QTLEffectsMat,qtl_effect[j])
      end
      QTL_index =vcat(QTL_index,qtl_index[j] .+ startlocus)
      startlocus += nLoci[j]
    end
    if size(G0,1)>0
        qtl_effect, QTLEffectsMat = transformEffects(numQTLOnChr, qtl_effect, gene_frequencyQTL, G0)
        println("G after transformation = ",0.5*QTLEffectsMat'QTLEffectsMat)
    end
    G = GenomeInfo(chr,nChromosome,mutation_rate,genotypeErrorRate,QTL_index,QTLEffectsMat)

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
                      nTraits::Int64=1,
                      G0 = [],
                      mutation_rate::Float64=0.0)
    nLoci          = fill(nLoci_each_chromosome,nChromosome)
    chrLength      = fill(chromosome_length,nChromosome)

    cstart         = chromosome_length/(nLoci_each_chromosome)/2
    cend           = chromosome_length-cstart
    map_position   = fill(collect(range(cstart,stop=cend,length=nLoci_each_chromosome)),nChromosome)
    gene_frequency = fill(fill(0.5,nLoci_each_chromosome),nChromosome)
    qtl_index      = [sample(1:nLoci_each_chromosome,qtl_each_chromosome,replace=false,ordered=true) for i in 1:nChromosome]
    qtl_effects    = fill(fill(0.0,qtl_each_chromosome,nTraits),nChromosome)
    for i in 1:nChromosome
        qtl_effects[i] = randn(qtl_each_chromosome,nTraits)
    end

    build_genome(nChromosome,chrLength,nLoci,gene_frequency,map_position,
         qtl_index,qtl_effects,nTraits,G0,mutation_rate)
end

function build_genome(nChromosome::Int64,
                      chromosome_length::Float64,
                      nLoci::Int64,
                      gene_frequency::Array{Float64,1},
                      map_position::Array{Float64,1},
                      mutation_rate::Float64=0.0,
                      qtl_index::Array{Int64,1}=Array{Int64,1}(undef, 0),
                      qtl_effects::Array{Float64,2}=Array{Float64,2}(undef, 0, 0),
                      G0=[],
                      genotypeErrorRate=0.0,
                      myCommon=common)

    nTraits           = size(qtl_effects,2)
    nLoci             = fill(nLoci, nChromosome)
    chromosome_length = fill(chromosome_length,nChromosome)
    gene_frequencyQTL = fill(gene_frequency[qtl_index],nChromosome)
    gene_frequency    = fill(gene_frequency,nChromosome)
    map_position      = fill(map_position,nChromosome)
    qtl_index         = fill(qtl_index,nChromosome)
    qtl_effects       = fill(qtl_effects,nChromosome)

    build_genome(nChromosome, chromosome_length, nLoci,
                gene_frequency, map_position,
                qtl_index, qtl_effects, nTraits, G0,
                mutation_rate, genotypeErrorRate)
end

function transformEffects(numQTLOnChr, qtlEffects, geneFreqQTL, G0)

    nTraits = size(G0, 1)
    qtlEffectsMat = Array{Float64, 2}(undef, 0, nTraits)
    for i in qtlEffects # qtlEffects is a vector of length 'numChr', each element of the vector is a matrix,
                        # with dimension numQTL x numTraits with effects for the QTL on that chrommosome
        qtlEffectsMat = [qtlEffectsMat; i] # we are concatenating those matrices into a numQTL x numTraits matrix
    end

    geneFreqQTLVec = []
    for i in geneFreqQTL # we are doing the same for the  genefrequencies of the QTL, but here we do not need a
        geneFreqQTLVec = [geneFreqQTLVec; i] # separate column per trait, as gene frequency does not depend on trait
    end

    d = 2 * geneFreqQTLVec .* (1 .- geneFreqQTLVec)
    D = diagm(0=>d) # creating diagonal mat with d on diagonal
    V = qtlEffectsMat'D * qtlEffectsMat
    L = cholesky(V).U'
    Li = inv(L)
    U = cholesky(G0).U

    A = qtlEffectsMat * Li'U

    AoM = Array{Array{Float64, 2}, 1}(undef, 0)
    numChr = length(numQTLOnChr)
    let
        k = 0
        for i = 1:numChr
            chrMat = Array{Float64, 2}(undef, numQTLOnChr[i], nTraits)
            for j = 1:numQTLOnChr[i]
                chrMat[j, :] = A[k + j, :]
            end
            push!(AoM, chrMat)
            k += numQTLOnChr[i]
        end
    end

    return AoM, A
end

export build_genome, transformEffects
export sampleFounders, sampleRan, sampleSel, samplePed,
       concatCohorts, cohortSubset,
       sampleBLUPSel, sampleDHOffspringFrom, sampleOneDHOffspringFrom
export getOurGenotypes, getOurPhenVals, getOurGenVals
export outputPedigree, outputGenData, outputHapData,
       outputGenData, outputCatData
export getIDs, getPedigree
export recode
export startrPop #deprecated

end # module
