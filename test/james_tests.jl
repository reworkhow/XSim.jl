# XSim
using XSim

build_genome(n_chr=2, n_loci=100)

build_demo()


build_phenome([500, 15];
            vg = [1 .5
                .5  1],
            h2 = [.3, .8])

build_phenome([10];
            vg = [1],
            h2 = [.3])
build_phenome(10;
            vg = 1,
            h2 = .3)
build_phenome([10];
            vg = 1,
            h2 = [.3])
build_phenome([10];
            vg = [1],
            h2 = .3)
build_phenome(10;
            vg = [1.0],
            h2 = .3)

CLEAR()

# Define a genome containing 4 loci
build_genome(n_loci=4, n_chr=1)
# Simulate a trait controlled by 2 loci
build_phenome(2, h2=.5)

# Load demo pedigree data frame
ped = DATA("pedigree")

# Create cohort
# Option 1: with simulated genotypes
cohort = Cohort(3)

# Option 2: with real genotypes
genotypes = DATA("genotypes")
cohort = Cohort(genotypes)
cohort|>get_genotypes

# Generate offsprings with the defiend pedigree
offsprings = mate(ped)

# Examine the genotypes
offsprings |> get_genotypes






build_demo()
ped


genotypes = get_genotypes(cohort)
n_loci  = size(genotypes, 2)
dt = hcat(get_IDs(cohort), genotypes) |> XSim.DataFrame
XSim.rename!(dt, vcat("ID", ["m$i" for i in 1:n_loci]))

vcat(["df", "63"], "h")

size(genotypes, 2)


using DataFrames, CSV




hcat(map_g, map_p)

DataFrames.columns(map_p)
get_columns(map_p)
map_p.columns




length(map_p)
["trait_$i" for i in ]

 @printf("trait_%d", 1)

cohort = Founders(20)

using Test
@test 3==3



QTL_effects = [1.0 .5
               0   1.0
               0   1.0
               1.0   0]
build_genome(n_chr=1, n_loci=4)
vg = [1 .5; .5 1]
build_phenome(QTL_effects, ve=vg, vp=[2 .5; .5 2])

args = Dict(:criteria => "EBV",
            :method => "GBLUP")


select(cohort, 10; args...)


round.(QTL_effects, )

x = []
x[1]

out = genetic_evaluation(cohort)

QTL_effects
# Initialize founders
n_sires = 3
n_dams  = 20
cohort_A   = Founders(n_sires)
cohort_B    = Founders(n_dams)

build_phenome(3, Vg=3, xx=6)

n_traits = 2
vg = [3.3 5;1 3]
XSim.handle_diagonal(3.3, n_traits)

|> XSim.Symmetric |> Array

n_traits[1]

XSim.matrix([3.5])[:, 1] |> XSim.diagm


function test(;a=missing, b=3)
    if ismissing(a)
        a = 4
    end
    return a + b
end

test(b=6)


n_traits = 2
h2 = .2
h2 = XSim.handle_h2(h2, n_traits)

v_src = [3]
v_src = [1.0 0.0;0.0 1.0]
term_src = "ve"
term_out = "vp"
v_src =  XSim.handle_diagonal(v_src, n_traits)

p = g + e
h2 = (p-e) / p  = g / (g + e)

XSim.get_Ve(n_traits, [1.0 0.0;0.0 1.0], .5)
XSim.get_Ve(n_traits, [1.0], .5)

XSim.diag(v_src)

h2 (g+e) = g

e = (g - gh2) / h2 = (1 - h2) * g

h2 * e = (1 - h2) * g = g - gh2

# g->e: e = (1 - h2) g / h2
v_out = ((ones(n_traits) .- h2) .* XSim.diag(v_src)) ./ h2
# e->g: g = e * h2 / (1 - h2)
v_out = (XSim.diag(v_src) .* h2) ./ (ones(n_traits) .- h2)
# p->g: g = p * h2
v_out = XSim.diag(v_src) .* h2


function infer_variances(v_src,
                         n_traits::Int64;
                         h2,
                         term_src::String,
                         term_out::String="optional")

    h2    = handle_h2(h2, n_traits)
    v_src = handle_diagonal(v_src, n_traits)

    if term_src == "vg"
        # out must be ve
        # g->e: e = (1 - h2) g / h2
        v_out = ((ones(n_traits) .- h2) .* diag(v_src)) ./ h2

    elseif term_src == "ve"
        # out must be vg
        # e->g: g = e * h2 / (1 - h2)
        v_out = (diag(v_src) .* h2) ./ (ones(n_traits) .- h2)

    elseif term_src == "vp"
        # out must be vg
        # p->g: g = p * h2
        v_out = diag(v_src) .* h2
    end

    v_out = n_traits == 1 ? v_out[1] : v_out

    return handle_diagonal(v_out, n_traits)
end

build_phenome(QTL_effects, )

XSim.handle_diagonal(1.0, 1)

h2 = g / (g + e)
(g + e) * h2 = g
h2g + h2e = g
(h2 - 1) g = -h2e

g = h2e / (1-h2)
p

if term_src == "vg"
    if    term_out == "ve"
        
    elseif term_out == "vp"

    end

    Ve = ((ones(n_traits) .- h2) .* diag(Vg)) ./ h2
Ve = n_traits == 1 ? Ve[1] : Ve

return handle_diagonal(Ve, n_traits)

p[2]



x = missing
ismissing()

XSim.get_Ve(3)

# vg, ve, vp, h2
A = [:x=>"a", :y=>"c"]
p = Pair(:vg, 9)

[:vg, :v34] in [:hg, :vg]
any([false, true])

- provide 1 arguments
h2
    traits are independent
    ve, vg = 1

vp or ve or vg
    h2 = 0.5

- provide 2 arguments
h2 and vp or ve or vg

vg





build_genome(n_chr=1, n_loci=400)
QTL_effects = [1.0 .5
               0   1.0
               0   1.0
               1.0   0]

Base.print_matrix(stdout, effects[begin:30])

build_phenome([5, 3])

GLOBAL("effects_QTLs")


args = Dict(:nA               => 5,
            :nB_per_A         => 10,
            :replace_A        => true,
            :replace_B        => false,
            :n_per_mate       => 1)

mate1 = cohort_A * cohort_B
mate2 = mate(cohort_A, cohort_B; args...)

mate1|> get_pedigree



XSim.@printf("%-30s %20.3f\n","residual variances:", QTL_effects)

(MCMCinfo.categorical_trait ? 1.0 : mme.R))

QTL_effects = [1.0 .5
               0   1.0
               0   1.0
               1.0   0]

build_phenome([3, 3])

Base.print_matrix(stdout, QTL_effects)

cohort = Cohort(3)
bvs = cohort |> get_BVs
adj = XSim.mean(bvs, dims=1)
bvs_0 = bvs .- adj

for i in 1:cohort.n
    cohort[i].val_g = bvs_0[i, :]
end


# Define parameters
args     = Dict(# mating
                :nA               => 3,
                :nB_per_A         => 5,
                :n_per_mate       => 2,
                :ratio_malefemale => 1.0,
                # selection
                :h2               => [.8, .5],
                :weights          => [.6, .4],
                # breeding
                :n_gens           => 5,
                :   => 3,
                :n_select_B => 20)

# Breeding program
sires, dams   = breed(sires, dams; args...)



dams = Founders(20)
out = genetic_evaluation(cohort, methods="BayesB", return_out=true) 

dt_p = get_phenotypes(dams, "JWAS")
idx = 3:6
XSim.allowmissing!(dt_p);
dt_p[idx, 2:end] .= missing

out = genetic_evaluation(dams, dt_p)

dt_p[:, "factor_1"] = [i for i in 1:4 for j in 1:5];
dt_p[:, "factor_2"] = [i for i in 1:2 for j in 1:10];

dt_p

out = genetic_evaluation(dams, dt_p,
                        model_equation="y1 = intercept + factor_1 + factor_2
                                        y2 = intercept + factor_1",
                        random_iid="factor_1", return_out=true)

20×2 Array{Float64,2}:
  28.1642    -74.8493
 -33.3485      9.32548
   4.23973    47.1613
   2.93291   -21.9782
  19.1834    -29.6798
  30.0874      0.426136
 -31.6773    -28.4247
  -5.30071    75.2142
  79.847      43.5332
  10.1532     78.9969
 -38.9385    -19.928
 -24.5641     24.0025
 -24.0252   -105.518
  54.0432      1.85657
  17.2012     15.3954
  40.7006     -2.13721
 -73.1188    -50.8392
 -36.7247    -46.4853
   6.45923    10.6944
 -25.3142     73.2333




file = DATA("demo_haplotypes.csv", header=false)
cohort = Cohort(file)


x = XSim.DataFrame()
XSim.nrow(x)

dtg = get_genotypes(cohort)
dtid = get_IDs(cohort)
hcat(dtid, dtg)


out["EBV"]

# XSim
using XSim
build_demo_small()
cohort = Cohort(100)
select(cohort, 3, criteria="GBLUP")


dams = Founders(20)
phenotypes, ve = get_phenotypes(cohort, "JWAS", return_ve=true)
phenotypes, ve = get_phenotypes(cohort,  return_ve=true)

traits = names(phenotypes)[2:end]

hcat([out["EBV_$x"][:, "EBV"] for x in traits]...)
traits
c = 5
"xged_$c"

out = genetic_evaluation(cohort)

# Working case
model_equation = "y1 = intercept + genotypes"
jwas_p         = get_phenotypes(cohort, "JWAS");
genotypes      = get_genotypes(cohort, "JWAS");
model          = XSim.JWAS.build_model(model_equation);
out = XSim.JWAS.runMCMC(model, jwas_p, methods="GBLUP")
​
# Not working case
function genetic_evaluation(cohort; model_equation)
    jwas_p        = get_phenotypes(cohort, "JWAS");
    genotypes     = get_genotypes(cohort, "JWAS");
    model         = XSim.JWAS.build_model(model_equation);
    return XSim.JWAS.runMCMC(model, jwas_p, methods="GBLUP")
end

eq = "y1 = intercept"

genotype = get_genotypes(cohort)
float.(genotype)
out = XSim.genetic_evaluation(cohort, model_equation=eq)



args_mate = Dict(:nB_per_A     => 5,
                 :n_per_mate   => 2)
progenies = mate(sires, dams; args_mate...)




select(cohort, 30)


using Statistics
cor(jwas_p[201:end, 2], out["EBV_y1"][201:end, 2])^2
cor(true_p, out["EBV_y1"][1:200, 2])^2
cor(vcat(true_p, jwas_p[201:end, 2]), out["EBV_y1"][:, 2])^2




n = cohort.n
n_traits = GLOBAL("n_traits")
ve_u = cholesky(GLOBAL("Ve")).U

hcat([ve_u * randn(n_traits) for _ in 1:n]...)' |> Array


g = get_genotypes(cohort)
get_BVs(cohort) 

hcat([XSim.cholesky(GLOBAL("Ve")).U * randn(2) for _ in 1:10]...)' |> Array



dt = DATA("demo_map.csv")

### conversion  of bp to cM based on reference
build_genome(dt)
build_phenome(dt, h2=[.3, .5])


hap = DATA("demo_haplotypes.csv", header=false)

x = Cohort()
for _ in 1:5
    x += XSim.sample(founders, 3) 
end

summary()
XSim.diag(GLOBAL("Vg") ./ (GLOBAL("Ve") + GLOBAL("Vg")))

f = [3]
f[[false]] = 0

build_genome(species="cattle")
build_demo()

XSim.gb.chromosome

n_chr = 2
n_loci = 3
n_row = n_chr * n_loci

# chr
chr = repeat([1:n_chr;], inner=n_loci)

# maf
dist = XSim.Normal(0, .05)
maf = .5 .- abs.(XSim.rand(dist, n_row))

# cM
cM = vcat([sort(XSim.uni_01(XSim.rand(dist, n_loci))) for _ in 1:n_chr]...)


build_genome(n_chr  = 3, n_loci= 5)


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


