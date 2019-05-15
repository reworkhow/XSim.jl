#Types and Methods for Simulating Genotypes of Animals

mutable struct Chromosome
    haplotype::Array{Int64,1}
    ori::Array{Int64,1}
    pos::Array{Float64,1}
    mut::Array{Float64,1}
end

mutable struct Animal
    genomePat::Array{Chromosome,1}
    genomeMat::Array{Chromosome,1}
    breedComp::Array{Float64,1}        ##################### 2015.6.4
    myID::Int64
    sireID::Int64
    damID::Int64
    phenVal::Array{Float64,1}
    genVal::Array{Float64,1}
    ebv::Float64
end

mutable struct Cohort
    animalCohort::Array{Animal,1}
    npMatrix::Array{Int64,2}
end

function Animal(mySire::Int64, myDam::Int64)
    my=Animal(Array{Chromosome}(undef,common.G.numChrom),Array{Chromosome}(undef,common.G.numChrom), Array{Float64}(undef,0),0,0,0,
              Array{Float64,1}(undef,0),Array{Float64,1}(undef,0),-9999.0)
    my.sireID = mySire
    my.damID  = myDam
    my.myID   = common.countId
    common.countId += 1
    for i in 1:common.G.numChrom
        my.genomePat[i]=Chromosome(Array{Int64}(undef,0),Array{Int64}(undef,0),Array{Float64}(undef,0),Array{Float64}(undef,0))
        my.genomeMat[i]=Chromosome(Array{Int64}(undef,0),Array{Int64}(undef,0),Array{Float64}(undef,0),Array{Float64}(undef,0))
    end
    return my
end

include("founders.jl")
include("nonfounders.jl")
include("checkCohort.jl")
