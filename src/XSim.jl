module XSim
    using  Distributions
    using  DataFrames
    using  CSV
    using  Random
    using  JWAS
    using  Printf
    using  LinearAlgebra

    import Base.length
    import Base.getindex
    import Base.print
    import StatsBase.sample

    tempPos = Array{Float64}(undef, 100000)
    tempOri = Array{Int64}(undef, 100000)
    tempMut = Array{Float64}(undef, 100000)

    """
    base type for genotype and Haplotype storage.
    Shouldn't be exported but needs to be defined. Used throughout the included src files.
    Originally an Int64
    """
    const AlleleIndexType = Int64

    # Objects
    include("objects/chromosome.jl")
    include("objects/animal.jl")
    include("objects/cohort.jl")
    include("objects/global.jl")
    # Core functions
    include("core/build.jl")
    include("core/build_genome.jl")
    include("core/build_phenome.jl")
    include("core/genome.jl")
    include("core/mate.jl")
    include("core/select.jl")
    # Interface
    include("interface/interface.jl")

    export Chromosome, Animal, Cohort
    export get_BVs, get_phenotypes, get_genotypes,
           get_QTLs,
           get_IDs, get_pedigree, get_DH
    export CELAR, SET, GLOBAL
    export build, build_genome, build_phenome
    export summary, summary_genome, summary_phenome
    export mate, select
    export sample_select, sample_random,
           random_mate, self_mate, all_mate, embryo_transfer
    export Global
end


# Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
# Variant identifier
# Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
# Base-pair coordinate