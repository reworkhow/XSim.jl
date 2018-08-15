# Type for storing Globals
mutable struct CommonToAnimals
    founders::Array{Animal,1}
    G::GenomeInfo
    countChromosome::Int64
    countId::Int64
    varRes::Float64
end

# Make object for storing globals
G = GenomeInfo(Array{ChromosomeInfo}(undef,0),0,0.0,0.0,[],[])
common = CommonToAnimals(Array{Animal}(undef,0),G,0,0,0.0)
