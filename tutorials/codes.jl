# ====== load genome and phenome ====== ====== ====== ====== ====== ====== ======


build_genome(species="cattle")
build_demo_small()
n_qtl = [2, 2]
Vg    = [ 1 .5
         .5  1]
build_phenome(n_qtl, Vg)
effects = Array([0 0 0 .5 0 .3 0 0 0 0
                .3 0 0 0 .8 0 0 .1 0 0]')
Vg    = [ 1 .5
         .5 1]
build_phenome(effects, Vg)





n_chr = 2
n_loci_chr = 5
n_loci = n_chr * n_loci_chr

chromosome = [i        for i in 1:n_chr for j in 1:n_loci_chr]
bp         = [10 * j   for i in 1:n_chr for j in 1:n_loci_chr]
cM         = [1.5 * j  for i in 1:n_chr for j in 1:n_loci_chr]
maf        = fill(0.5, n_loci)
rate_mutation = 0.0
rate_error    = 0.0
build_genome(chromosome, bp, cM, maf, rate_mutation, rate_error)


reference = XSim.data("cattle_map")
build_genome(reference)



build_phenome(3, [3, 6.4])
build_phenome([3, 2], 5.0)


# ====== Haplotype ====== ====== ====== ====== ====== ====== ======
using XSim
import Random
Random.seed!(95616)
CLEAR()
build_demo_small()

genotypes = XSim.data("genotypes")

cohort = Founders(genotypes)

root      = dirname(dirname(pathof(XSim)))
filepath  = joinpath(root, "data", "demo_genotypes.csv")
cohort    = Founders(filepath)

get_genotypes(cohort)
get_QTLs(cohort)
get_BVs(cohort)
get_pedigree(cohort)

args = Dict(:n_per_mate      => 4,
            :n_gens          => 5,
            :ratio_malefemale=> 2,
            :h2              => [.8, .8])

males, females = breed(cohort, cohort; args...)

summary(males + females)

get_QTLs(males + females) |> get_MAF
get_QTLs(cohort) |> get_MAF



# ====== Basic ====== ====== ====== ====== ====== ====== ======
using XSim
# import Random
# Random.seed!(95616)

build_demo_()
cohort = Founders(5)
genetic_evaluation(cohort)


n_sires = 5
n_dams  = 2
sires   = Founders(n_sires)
dams    = Founders(n_dams)
get_pedigree(sires, "JWAS")
get_phenotypes(sires)




args_mate     = Dict(:n_per_shared => n_dams,
                     :n_per_mate   => 2)
progenies     = mate(sires, dams; args_mate...)

args_select   = Dict(:h2 => [.5, .5])
progenies_sel = select(progenies, 10; args_select...)

args_breed    = Dict(:n_gens       => 5,
                     :n_select     => 10)
sires, dams   = breed(sires, dams; args_breed..., args_mate..., args_select...)

# for i in 1:5
#     progenies = mate(sires, dams; args_mate...)
#     progenies = select(progenies, 10; args_select...)
#     sires, dams = progenies, progenies
# end


# ====== Simple ====== ====== ====== ====== ====== ====== ======
using XSim
XSim.CLEAR()

import Random
Random.seed!(95616)

build_demo()

n_sires       = 20
dams_per_sire = 10
n_dams        = n_sires * dams_per_sire

args_mate   = Dict(:n_per_shared       => dams_per_sire,
                   :n_per_mate         => 2,
                   :ratio_malefemale   => 1,
                   :replace_shared     => false,
                   :replace_per_shared => false)
args_select = Dict(:h2                 => [.5, .5],
                   :is_random          => false)
args_breed  = Dict(:n_gens             => 10,
                   :n_select_males     => n_sires)

sires       = Founders(n_sires)
dams        = Founders(n_dams)

@time sires, dams = breed(sires, dams; args_breed..., args_mate..., args_select...)


@time XSim.set_BV!(dams[1])
@time XSim.set_BV2!(dams[1])

dams[1].genome_sire


# Results
summary(sires)
summary(dams)
summary(sires + dams)
get_MAF(get_QTLs(sires+dams))

# ====== Cross breeds ====== ====== ====== ====== ====== ====== ======
using XSim
build_demo()

# Small breed
n_sires       = 50
dams_per_sire = 10
n_dams        = n_sires * dams_per_sire
args          = Dict(# Mating
                     :n_per_shared     => dams_per_sire,
                     :n_per_mate       => 2,
                     :ratio_malefemale => 1,
                     # Selection
                     :h2               => [.8, .2],
                     :is_random        => false,
                     # Breeding
                     :n_gens           => 10,
                     :n_select_males   => n_sires)
# Breed A
sires_A         = Founders(n_sires)
dams_A          = Founders(n_dams)
sires_A, dams_A = breed(sires_A, dams_A; args...)

# Large breeds
n_sires        = 100
dams_per_sire  = 20
n_dams         = n_sires * dams_per_sire
args[:n_per_shared]   = dams_per_sire
args[:n_select_males] = n_sires

# Breed B
sires_B         = Founders(n_sires)
dams_B          = Founders(n_dams)
sires_B, dams_B = breed(sires_B, dams_B; args...)

# Breed C
sires_C         = Founders(n_sires)
dams_C          = Founders(n_dams)
sires_C, dams_C = breed(sires_C, dams_B; args...)

# Rotation
args_rotate          = Dict(:n_pop            => 2000,
                            :n_per_mate       => 2,
                            :ratio_malefemale => 1)
# Rotation (G1)
males_G1, females_G1 = mate(sires_B, dams_C; args_rotate...)

# Rotation (G2)
males_G2, females_G2 = mate(sires_A, females_G1; args_rotate...)

# Rotation (G3)
males_G3, females_G3 = mate(sires_C, females_G2; args_rotate...)

