# using DataFrames: DataAPI
using XSim

dt = DATA("demo_map.csv")

build_

DATA("demo_map.csv")
### conversion  of bp to cM based on reference
build_genome(dt)
build_phenome(dt, h2=[.3, .5])


hap = DATA("demo_haplotypes.csv", header=false)



summary()
XSim.diag(GLOBAL("Vg") ./ (GLOBAL("Ve") + GLOBAL("Vg")))

f = [3]
f[[false]] = 0

build_genome(species="cattle")
build_demo()

XSim.gb.chromosome

n_chr = 2
n_marker = 3
n_row = n_chr * n_marker

# chr
chr = repeat([1:n_chr;], inner=n_marker)

# maf
dist = XSim.Normal(0, .05)
maf = .5 .- abs.(XSim.rand(dist, n_row))

# cM
cM = vcat([sort(XSim.uni_01(XSim.rand(dist, n_marker))) for _ in 1:n_chr]...)


build_genome(n_chr  = 3, n_marker= 5)


rep(34, 6)


XSim.Normal(0, 1)
fill(1:5)


#
ref = XSim.DATA("genome_cattle.csv")
dt = XSim.DATA("demo_map.csv")

columns = names(dt)

has_chr, has_bp, has_cM, has_maf = in(columns).(["chr", "bp", "cM", "maf"])

species = "cattle"
ref = load_ref(species)


# Infer cM
if all([has_chr, has_cM])
    if ismissing(ref)
        # pass
    else
        if has_bp
            dt = add_cM_by_ref!(dt, ref)
            LOG("The provided genetic distances will be replaced with ones infered from preloaded linkage maps", "warn")
        else
            LOG("Missing required column 'bp'", "error")
        end
    end

elseif all([has_chr, !has_cM, has_bp])
    if ismissing(ref)
        dt = add_cM_by_bp!(dt)
    else
        dt = add_cM_by_ref!(dt, ref)
    end

else
    LOG("Missing required columns", "error")
end

# bp and MAF
dt.maf = has_maf ? dt.maf : fill(0.5, nrow(dt))
dt.bp  = has_bp  ? dt.bp  : fill(0,   nrow(dt))

# build genome
build_genome(dt.chr,
             dt.bp,
             dt.cM,
             dt.maf;
             args...)









# if bp or cM
build_genome("file.csv")


###
jwas_P   = get_phenotypes(cohort) |> XSim.DataFrame


using XSim
using Lazy
#  ====== test dict arguments ====== ====== ====== ====== ====== ====== ======
# https://en.wikibooks.org/wiki/Introducing_Julia/Dictionaries_and_sets
ddd = Dict(:a => 3, :b => 4, :c => 5)
Pair(:a=>3)

k = collect(keys(ddd))
v = collect(values(ddd))

idx_supported = findall(in([:b, :c]), collect(keys(ddd)))
Dict(k[i] => v[i] for i in idx_supported)

subset_dict(ddd, [:b, :c])

ddd[idx_supported]
pars= Iterators.Pairs(:a=>1, :b=>2)


#  ====== WHY BVs are mostly negative?? ====== ====== ====== ======
build_demo()
GLOBAL("effects_QTLs")
Cohort(30)


20 5*4*1*1


a = [1, 2, 3, 4]

repeat(a, outer=5)

# ======  ====== ====== ====== ====== ====== ====== ======








GLOBAL("silent")
SILENT(true)
males
nothing
# use dict to return summary
# get_MAF(sires)


get_BVs(sires)
get_BVs(dams)

get_genotypes(dams)
get_QTLs(dams)
get_QTLs(sires)
GLOBAL("effects_QTLs")



scores = [32,40, 50, 20, 10]
id = [1,2,3,4,5]
(1:5)[sortperm(scores, rev=true)][1:3]

a = [1.403 1.377]
b = [1.106 1.004]
a./b

a .- b


original = get_BVs(f1_dams)
sel = original[1:5, :]

bv = round.(XSim.mean(original, dims=1), digits=3)
bv_sel = round.(XSim.mean(sel, dims=1), digits=3)

var_g      = round.(XSim.var(bvs,  dims=1), digits=3)

bv
sum_f1 = summary(f1_dams)
sum_f2 = summary(f1_dams[1:5])


(sum_f2["mu_g"] - sum_f1["mu_g"]) / sum(sum_f1["mu_g"])



ratio = 1.0
n_animals = length(c)


ratio = m / f
n = m + f

m = n - f
isa(1/5, Float64)
ratio = 1
n_males = convert(Int64, round(n_animals * ratio / (ratio + 1)))
n_female = n_animals - n_males

c[1:n_males]
c[(n_males+1):end]


#  ====== parallelism ====== ====== ====== ====== ====== ====== ======
Threads.nthreads()
# julia --threads 4
run(`pwd`)
ENV
ENV["JULIA_NUM_THREADS"]="4"
ENV["SHELL"]

# @time Threads.@threads for i in 1:1000000
#            randn(10)
#        end
# @time for i in 1:1000000
#            randn(10)
#       end
# julia> @time self_mate(founders, 10000)
#   1.988459 seconds (5.06 M allocations: 185.074 MiB, 31.52% gc time)


#  ====== manual case ====== ====== ====== ====== ====== ====== ======
build_demo()
founders = Founders(100)
summary(founders)

mate(founders, n=3)
mate(founders, founders, 3)

#  ====== reference case ====== ====== ====== ====== ====== ====== ======

build_genome(species="cattle")


SILENT(true)

x = Cohort()

n_qtl = [2, 2]
Vg    = [ 1 .6
         .5  1]
params = [n_qtl, Vg]
build_phenome(params...)


build_phenome(3, [3, 6.4])
build_phenome([3, 2], 5.0)

function test(;a::Int64, b::Int64)
    return a+b
end

dd = Dict(:a=>3, :b=>4)
test(;dd...)

test(a=3, b=5)

[size(n_qtl)..., size(Vg)...]


h2 = [0.8, 0.9]
weights = [1.0, 0.0]
n = 100
n_sel = 10
SILENT(false)

founders = Founders(50)
@> f1 = self_mate(founders, n) select(n_sel, h2=h2, weights=weights)
@> f2 = self_mate(f1, n)       select(n_sel, h2=h2, weights=weights)
@> f3 = self_mate(f2, n)       select(n_sel, h2=h2, weights=weights)

summary(founders)["Mu_g"]
summary(f1)["Mu_g"]
summary(f2)["Mu_g"]
summary(f3)["Mu_g"]

founders
f1

#  ====== Genotype ====== ====== ====== ====== ====== ====== ======
g = get_genotypes(founders)

i1 = XSim.CSV.read("test.csv", DataFrame, header=false, missingstrings=["-1", "9"])
i1

co = Cohort(i1)
get_genotypes(co)
idx_row = 3:96
idx_col = 50:69



# Validation
sum(get_genotypes(co)[idx_row, idx_col])
sum(sum.(eachrow(i1[idx_row, idx_col])))

genotypes = Array(i1)
freq = sum(genotypes, dims=1) / (2 * size(genotypes, 1))
maf = min.(freq, 1 .- freq)



# [0 1 1]2/6 = 0.33




length(co[1].genome_sire[2].haplotype)
co[1].genome_sire[2].haplotype
co[1].genome_sire[1].haplotype
summary_genome()

#  ====== Test chr ====== ====== ====== ====== ====== ====== ======
chr = [1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4]
idx_each_chr = [chr .== c for c in unique(chr)]
hcat(findfirst.(idx_each_chr), findlast.(idx_each_chr))

####
using  SparseArrays

spzeros(3, 3)
#  ====== Test function  ====== ====== ====== ====== ====== ====== ======
sp = sparse([1, 3], [1, 2], [5, 3])
XSim.matrix(sp)
mat = [3 5; 1 2]
sparse(sp)

#  ====== Old XSim ====== ====== ====== ====== ====== ====== ======

# using XSim, CSV, DataFrames

# dt = CSV.read("genome_pig.csv", DataFrame)

# numChr = length(unique(dt.chr))
# chrLength = [round(max(dt.cM[dt.chr.==c]...)/100, digits=2) + .01 for c in unique(dt.chr)]
# numLoci = [count(==(c), dt.chr) for c in unique(dt.chr)]
# nTraits = 2
# numQTLOnChr = [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,1,1]
# numQTL = 2
# qtlIndex = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[], [500],[500]]
# geneFreq = [fill(0.5, n) for n in numLoci]
# mPos = [dt.cM[dt.chr .== chr]./100 for chr in unique(dt.chr)]
# G0 = [1. 0.5;0.5 2.]

# qtlEffects = Array{Array{Float64,2},1}(undef,0)
# push!(qtlEffects,randn(1,nTraits))
# push!(qtlEffects,randn(1,nTraits))

# build_genome(numChr,
#             chrLength,
#             numLoci,
#             geneFreq,
#             mPos,
#             qtlIndex,
#             qtlEffects,
#             nTraits,
#             G0)

# @time founders = sampleFounders(50)
# #  ====== Old XSim ====== ====== ====== ====== ====== ====== ======


using Revise
# includet("test/james_tests.jl")
includet("test/runtests.jl")




using CSV, DataFrames, LinearAlgebra, StatsBase

dt = CSV.read("data/genome_cattle.csv", DataFrame)
# dt = CSV.read("data/genome_pig.csv", DataFrame)
rand(XSim.Bernoulli(.5))
randn!(XSim.Bernoulli(.5))

# Load genome
try
    SET("chromosome", dt.chr)
    SET("bp"        , dt.bp)
    SET("cM"        , dt.cM)
catch
    error("Missing required columns")
end

# Define MAF
if all(in(names(dt)).(["maf"]))
    SET("maf", dt.maf)
else
    SET("maf", fill(0.5, GLOBAL("n_loci")))
end



function fa(a::Int64; b::Int64)
    return a + b*2
end

function fa(a::Int64; c::Union{Array{Int64, 1}, Int64})
    return a + c*3
end

fa(1, b=3)
fa(1, c=3)



Vg = [1 0.5
      0.5 1]
h2 = [.5, .8]

Vg = matrix(.7)
h2 = .5

n_traits = size(Vg)[2]

ve = ((ones(n_traits) .- h2) .* diag(Vg)) ./ h2
handle_diagonal(ve[1], n_traits)


inputs = ve
# Formate variances to 2-D matrix
inputs = matrix(inputs)

# Cast variants of variances to a 2-D array
if length(inputs) == 1
    # When variances is a scaler, assign it as the diagonal of variances
    inputs = diagm(fill(inputs, n_traits))
else
    inputs = matrix(inputs)
    if size(inputs)[2] == 1
        # When variances is a vector, assign it as the diagonal of variances
        inputs = diagm(inputs[:, 1])
end

if size(inputs)[2] != n_traits
    error("Dimensions don't match between n_traits and variances/h2")
end

    return inputs






vg = h2*vg  + h2 * ve
(1 - h2)Vg / h2 = ve

effects = [0.3 0.5
            0 0.2
         -1.3 0.4]
freq = [0.5, 0.5, 0.5]

effects_sc = scale_effects(effects, freq, Vg)
get_Vg(effects_sc, freq)

Vg = matrix(3.0)
effects = matrix([0.3, 0.5, 0])
freq = [0.5, 0.5, 0.5]

effects_sc = scale_effects(effects, freq, Vg)

get_Vg(effects_sc, freq)



vec1 = [1,2,3]
vec2 = 34

isa(vec1, Array)


freq = [1,2,3,3,3,3,2,5,3,4]
sample(1:10, 11, replace=true)


function TU(mat ::Union{Array{Int64, 1}, Int64})
    print(mat)
end


TU([3, 6])

function ab(mat ::Array{Int64, N} where N = [1, 2])
    return mat
end

vfs([1 2 4
     1 2 6])


x = matrix([3, 5])
diagm(x)



isa([3.0], Array{Float64, 1})

fla(Float64)


mat = [1 2
       3.3 4]

       mat
diag(1)
isa(mat, Array)
isa(3, Array)

c("Gender", "Age") %in% names(df)

function funB(;
    a, b, c, d=3)
    return a+b+c+d
end

a = Array{Float64, 2}([1.0], 1, 1)
a[1, 1] = 3
(, 1, 1)

hcat(Diagonal([1])...)

funB(a       =3,
    b = 5, 
    c =     3)


compute_rec!(dt)
length(mat)
Diagonal()

mutable struct test
    a
    b
end

T = test(1, 3)
setfield!(T, Symbol("a"), 99)


# ```Append columns to dt```
# function compute_rec_prob(chromosomes::Array{Int64,   1},
#                           cM         ::Array{Flaot64, 1})
#     # --- Compute recombination rate
#     # Get diff cM
#     cM_diff = vcat(0, diff(cM))
#     n_loci = GLOBAL("n_loci")
#     for i in 1:n_loci
#         # diff(cM) < 0 when across chromosome
#         if cM_diff[i] < 0
#             cM_diff[i] = cM[i]
#         end

#         # diff(cM) == 0 when markers in the same LD block
#         if i < length(cM_diff) && cM_diff[i + 1] == 0
#             cM_diff[i + 1] = cM_diff[i]
#         end
#     end

#     # --- Compute standardized probability
#     probs = zeros(n_loci)
#     for chr in unique(chromosomes)
#         idx = chromosomes .== chr
#         probs[idx] = round.(    cM_diff[idx]/
#                             sum(cM_diff[idx]),
#                             digits=5)
#     end
#     return probs
# end


