# Type for storing Globals
mutable struct CommonToAnimals
    founders::Array{Animal,1}
    G::GenomeInfo
    countChromosome::Int64
    countId::Int64
    LRes::Array{Float64,2} #Cholesky of residual covMat put here instead of a scalar for varRes as in single trait version
end

# Make object for storing globals
G = GenomeInfo(Array{ChromosomeInfo}(undef,0),0,0.0,0.0,[],Array{Float64,2}(undef,0,0))
common = CommonToAnimals(Array{Animal}(undef,0),G,0,0,Array{Float64,2}(undef,0,0))
