"""
Author: James Chen, Virginia Tech <niche@vt.edu>
Date:   2022-12-22
Description: This script is presented for demonstrating the XSim.jl package in the 179th AGRIX Forum 百博智慧大講堂
"""

# solve a linear equation
y = [5 22 6]'
x = [0 0 1 0
    0 0 0 1
    1 1 0 0]

x \ y

# libraries
using Distributions
using Plots

# beta distribution with parameters (m - k + 1, k + 1)
m = 100;
k = 20;
dist_1 = Beta(m - k + 1, k + 1)

m = 100;
k = 80;
dist_2 = Beta(m - k + 1, k + 1)

# plotting
x = 0:0.01:1 # from 0 to 1 with step size 0.01
plot(x, pdf.(dist_1, x), color="red", label="Beta(81, 21) - 80% zero-effect")
plot!(x, pdf.(dist_2, x), color="blue", label="Beta(21, 81) - 20% zero-effect")
plot!(title="PDF of Beta Distribution",
    xlabel="\\pi", ylabel="pdf\\( \\pi \\)", legend=:top)


# XSim start ------------------------------------
using XSim

# build genome
build_genome(n_chr=2, n_loci=10000)
map = DATA("map")
build_genome(map)
build_genome(species="cattle")

# build phenome
build_phenome(10)
build_phenome(10; h2=0.3)
build_phenome(10;
    h2=0.8,
    n_traits=2,
    vg=[1 0.2
        0.2 1])
build_phenome(map)

# cohort ------------------------------------
cohort = Cohort(100)
get_genotypes(cohort)
get_QTLs(cohort)
get_pedigree(cohort, sort=true)
get_BVs(cohort)
get_phenotypes(cohort, h2=0.3)
get_phenotypes(cohort, h2=0.9)

# mating ------------------------------------
cohort_sire = Cohort(10)
cohort_dam = Cohort(100)

get_pedigree(offspring, sort=true)
get_pedigree(cohort_sire + offspring, sort=true)

offspring_1 = mate(cohort_sire, cohort_dam)
offspring_5 = mate(cohort_sire, cohort_dam, nB_per_A=5)
get_pedigree(cohort_sire * cohort_dam, sort=true)
get_pedigree(offspring_1 + cohort_sire, sort=true)

# selection ------------------------------------
select(offspring_5, 10; h2=0.3)
select(offspring_5, 10; h2=0.9)

select(offspring_5, 10; h2=0.3, criteria="random")
select(offspring_5, 10; h2=0.3, criteria="phenotypes")
select(offspring_5, 10; h2=0.3, criteria="EBV")
# select by provided phenotypes
phenotype = get_phenotypes(offspring_5, h2=0.3)
select(offspring_5, 10; criteria=phenotype)

args = Dict(
    :weights => [0.5, 0.8],
    :criteria => "EBV",
    :ve => [1 0.3
        0.3 1])
select(offspring_5, 10; args...)


# JWAS ------------------------------------
cohort_gs = cohort_sire + offspring_5
ped_gs = get_pedigree(cohort_gs)
dt_p = get_phenotypes(cohort_gs, "JWAS")
dt_p[:, "sire"] = sire_id(ped_gs)

# sire model
out = genetic_evaluation(cohort_gs, dt_p,
    model_equation="y1 = intercept + sire
                    y2 = intercept",
    random_str="sire",
    return_out=true)

# cow breeding ------------------------------------
using DataFrames
using XSim

# define background
build_genome(species="cattle")
build_phenome(50; h2=0.3)

# cohorts
bulls = Cohort(50)
cows = Cohort(1000)

# avoid lethal alleles
idx_lethal = [30, 1345, 2030]
genotypes = get_genotypes(bulls)[:, idx_lethal]
idx_select = sum(genotypes .== 2, dims=2) .== 0
idx_select = vec(idx_select)
bulls = bulls[idx_select]

# mating
daughters = mate(bulls, cows; nB_per_A=10)
get_pedigree(daughters, sort=true)

# generate phenotypes
pedigree = get_pedigree(daughters)
milk_yields = get_phenotypes(daughters, "JWAS")
milk_yields[:, "sire"] = sire_id(pedigree)

# GS sire model
cohort_gs = bulls + daughters
out = genetic_evaluation(cohort_gs, milk_yields,
    model_equation="y1 = intercept + sire",
    methods="BayesC",
    random_str="sire",
    return_out=true)

# inspect the results (EBV, PEV:prediction error variances)
out["EBV_y1"][1:bulls.n, :]

# maize breeding NAM ------------------------------------
CLEAR()

# define genetic background
build_genome(species="maize")
build_phenome(100; h2=0.3)

# founders
diverse_parents = Founders(25)
common_parents = Founders(1)
diverse_parents |> get_genotypes

# derive F1s
F1 = Cohort()
for parent in diverse_parents
    F1 += common_parents * parent
end
F1 |> get_pedigree
F1 |> get_genotypes

# Each family produce 200 progenies by selfing
args = Dict(
    # mating
    :n_per_mate => 10,
    :scheme => "selfing",
    # selection
    :criteria => "phenotypes",
    # breed
    :n_gens => 4,
    # single-seed decent
    :n_select => 1)

NAM = Cohort()
for family in F1
    F2 = mate(family; n_per_mate=200, scheme="selfing")
    for seed in F2
        NAM += breed(seed; args...)
    end
end

# visualization
using MultivariateStats
using Plots
# run PCA on g_nam matrix
g_nam = get_genotypes(NAM)
pca = fit(PCA, g_nam', maxoutdim=3)
# NAM genotypes
g_nam = get_genotypes(NAM)
PCs_nam = predict(pca, g_nam')
# parental genotypes
g_p = get_genotypes(diverse_parents)
PCs_p   = predict(pca, g_p')
# common parental genotypes
g_c = get_genotypes(common_parents)
PCs_c = predict(pca, g_c')

# plot 1:
plot(size=(1000, 1000))
scatter!(PCs_nam[1, :], PCs_nam[2, :], label="NAM")
scatter!(PCs_p[1, :], PCs_p[2, :], color=:red, label="diverse parents",
    markersize = 10, markerstrokewidth=1, alpha=0.7)
scatter!(PCs_c[1, :], PCs_c[2, :], color=:yellow, label="common parents",
    markersize=10, markerstrokewidth=1, alpha=0.7)

# plot 2:
n_nam = 200
n_p = 25
colors = palette(:Dark2_8);
plot(size=(1000, 1000), legend=false)
theme(:orange)
for i in 1:n_p
    c = i % length(colors) + 1
    idx = (i - 1) * n_nam + 1 : i * n_nam
    # offspring
    scatter!(PCs_nam[1, idx], PCs_nam[2, idx], color=colors[c],
        label="NAM_$i",
        markersize=5, markerstrokewidth=0, linewidth=.1, alpha=.5)
end
for i in 1:n_p
    c = i % length(colors) + 1
    # parents
    scatter!(vec([PCs_p[1, i]]), vec([PCs_p[2, i]]), color=colors[c],
        label="diverse parents",
        markersize=7, markerstrokewidth=1, alpha=1)
end
scatter!(PCs_c[1, :], PCs_c[2, :], color=:black, label="common parents",
    markersize=10, markerstrokewidth=1, alpha=1, markerstrokecolor=:white)
