using XSim

build_demo_small()

# Initialize a population with 1,500 founders
founders = Founders(1500)
# Let founders random mate with each other
# for 1,000 generations
for _ in 1:1000
    founders = mate(founders)
end
# Drop the population size to 100 individuals and
# continue the random mating for another 15 generations
for _ in 1:15
    founders = mate(founders[1:100])
end
sires_base = dams_base = founders


#Simulate three pure breeds
args_A  = Dict(# Mating
               :nA               => 50,
               :nB_per_A         => 10,
               :n_per_mate       => 2,
               :ratio_malefemale => 1,
               # Selection
               :criteria         => "random",
               # Breeding
               :n_gens           => 10,
               :n_select_A       => 50,
               :n_select_B       => 500)
args_BC = Dict(# Mating
               :nA               => 100,
               :nB_per_A         => 20,
               :n_per_mate       => 2,
               :ratio_malefemale => 1,
               # Selection
               :criteria         => "random",
               # Breeding
               :n_gens           => 10,
               :n_select_A       => 100,
               :n_select_B       => 2000)
# Breed A, B, and C
sires_A, dams_A = breed(sires_base, dams_base; args_A...)
sires_B, dams_B = breed(sires_base, dams_base; args_BC...)
sires_C, dams_C = breed(sires_base, dams_base; args_BC...)


# Rotation parameters
args_X = Dict(:n_pop            => 2000,
              :n_per_mate       => 2,
              :ratio_malefemale => 1)
args_A[:n_gens]  = 1
args_BC[:n_gens] = 1
# Rotation (G1)
sires_A1, dams_A1    = breed(sires_A, dams_A; args_A...)
sires_B1, dams_B1    = breed(sires_B, dams_B; args_BC...)
sires_C1, dams_C1    = breed(sires_C, dams_C; args_BC...)
males_G1, females_G1 = mate(sires_B, dams_C;  args_X...)
# Rotation (G2)
sires_A2, dams_A2    = breed(sires_A1, dams_A1; args_A...)
sires_B2, dams_B2    = breed(sires_B1, dams_B1; args_BC...)
sires_C2, dams_C2    = breed(sires_C1, dams_C1; args_BC...)
males_G2, females_G2 = mate(sires_A1, females_G1;
                            args_X...)
# Rotation (G3)
sires_A3, dams_A3    = breed(sires_A2, dams_A2; args_A...)
sires_B3, dams_B3    = breed(sires_B2, dams_B2; args_BC...)
sires_C3, dams_C3    = breed(sires_C2, dams_C2; args_BC...)
males_G3, females_G3 = mate(sires_C2, females_G2;
                            args_X...)



build_genome("map.csv")
build_phenome("map.csv",
              vg = [ 1 .5; .5  1],
              h2 = [0.3, 0.7])
# Initialize a population with 1,500 founders
founders = Founders(1500)
# Let founders random mate with each other
# for 1,000 generations 
for _ in 1:1000
    founders = mate(founders)
end
# Drop the size to 100 individuals and
# continue the random mating for another 15 generations
for _ in 1:15
    founders = mate(founders[1:100])
end
# Subset 26 founders to become the base population
common_parents  = founders[1]
diverse_parents = founders[2:26]


# Cross each diverse parent with the common parent
F1 = Cohort()
for parent in diverse_parents
    F1 += common_parents * parent
end
# Each family produce 200 progenies by selfing
args = Dict(# Mating
            :n_per_mate => 200,
            :is_selfing => true,
            # Selection
            :criteria         => "random",
            # Breed
            :n_gens     => 6,
            :n_select   => 200)
NAM = Cohort()
for family in F1
    NAM += breed(family, args...)
end



v_src = XSim.handle_diagonal([1.0], 3)



# Load XSim
using XSim

# Simulate genome with 10 chromosomes, and 100 markers are located on each chromosome.
build_genome(n_chr=10, n_marker=100)
# Simulate two independent traits controlled by 3 and 8 QTLs, respectively.
build_phenome([3, 8])

# Initialize founders
n_sires = 3
n_dams  = 20
sires   = Founders(n_sires)
dams    = Founders(n_dams)

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
                :n_select_A       => 3,
                :n_select_B       => 20)

# Breeding program
sires_new, dams_new   = breed(sires, dams; args...)

# Inspect the results
summary(sires + dams)
summary(sires_new + dams_new)



