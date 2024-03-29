module XSim
using Base
using CSV
using Distributions
using DataFrames
using CSV
using Random
using JWAS
using Printf
using LinearAlgebra
using SparseArrays
using StatsBase

import Base.length
import Base.getindex
import Base.print
import Base.sort
import Base.push!
import Base.pop!
import Base.popfirst!
import StatsBase.sample
import StatsBase.mean
import StatsBase.var

const AlleleIndexType = Int64

# Objects
include("objects/chromosome.jl")
include("objects/animal.jl")
include("objects/cohort.jl")
include("objects/global.jl")
include("objects/data.jl")
include("objects/gs_pool.jl")
include("objects/checkpoint.jl")
# Core functions
include("core/build.jl")
include("core/build_genome.jl")
include("core/build_phenome.jl")
include("core/genome.jl")
include("core/mate.jl")
include("core/select.jl")
include("core/genetic_evaluation.jl")
include("core/variance.jl")
# Interface
include("interface/breed.jl")
include("interface/build_demo.jl")
include("interface/others.jl")
include("interface/beta.jl")

# Initialize global
gb = GB()

# Base
export sort, sample, length
# Core struct
export Chromosome, Animal, Cohort, Founders
# Breeding structs
export GS_pool, Checkpoint, CheckpointList, update!, cost!
# Getters
export get_BVs, get_EBVs, get_phenotypes, get_genotypes,
    get_QTLs,
    get_IDs, get_pedigree, get_DH
# Cohort
export center_BV!
export get_Vg, get_MAF, scale_effects
export genetic_evaluation, GBLUP
# Global
export CLEAR, SET, GLOBAL, LOG, SILENT, DATA, PATH,
    GLobal, gb, GET_LINES
# Variance
export handle_h2, handle_variance, sample_qtls
# Build
export build, build_genome, build_phenome, build_demo, build_demo_small
# Summary
export summary, summary_genome, summary_phenome
# Core
export mate, select
# Interface
export breed, sample_select, sample_random,
    random_mate, self_mate, all_mate, embryo_transfer,
    save_map, sample_matings
XSim.LOG("--- --- --- XSim v2.1.4 --- --- --- ")
XSim.LOG("Last update: Nov 23, 2022")
XSim.LOG("Developers: James Chen <niche@vt.edu>
            Hao Cheng <qtlcheng@ucdavis.edu>,
            Rohan Fernando <rohanluigi@icloud.com>,
            Mike Keehan <m.a.keehan@gmail.com>,
            Chris Stricker <stricker@genetics-network.ch>")
XSim.LOG("User manual: https://reworkhow.github.io/XSim.jl")
end

# Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
# Variant identifier
# Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
# Base-pair coordinate
