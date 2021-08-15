using XSim, CSV, DataFrames
# Genome --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Option 1: Defined by literature-reported species
build_genome(species="cattle")

# Option 2: Defined by dataframes
build_genome("map.csv")
# or refer missing genetic positions by literatures
build_genome("map.csv"; species="cattle")

# Option 3: Defined by dataframe
using CSV, DataFrames
dataframe = CSV.read("map.csv", DataFrame)
build_genome(dataframe)
# or refer missing genetic positions by literatures
build_genome(dataframe; species="cattle")

# Option 4: Defined manually
ch  = [1,    1,     2,    2,    2]
bp  = [130,  205,   186,  503,  780]
cM  = [85.7, 149.1, 37.4, 83.6, 134.3]
maf = [0.5,  0.5,   0.5,  0.5,  0.5]
build_genome(ch  = [1,    1,     2,    2,    2],
             bp  = [130,  205,   186,  503,  780],
             cM  = [85.7, 149.1, 37.4, 83.6, 134.3],
             maf = [0.5,  0.5,   0.5,  0.5,  0.5])

# Phenome --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Option 1: Defined by Vg and n QTLs
n_qtl = [1, 2]
Vg    = [ 1 .5
         .5  1]
build_phenome(n_qtl, Vg)

# Option 2: Defined by Vg and QTL effects
effects = Array([0 .5 0 0  0
                .3  0 0 0 .8]')
build_phenome(effects, Vg)

# Option 3: Defined by dataframes
build_phenome("map.csv", Vg)

# Option 4: Defined by dataframes
using CSV, DataFrames
dataframe = CSV.read("map.csv", DataFrame)
build_phenome(dataframe, Vg)

### Manually define genome and Phenome
# Define a phenome by a user-provided file by numbers of QTL for each trait (n_qtl) and the genetic variances (Vg); Phenotypes are simulated based on the specified heritability or residual variance.
n_qtl = [5, 10]
Vg    = [ 1 .5
         .5  1]

build_genome(ch =[1,    1,     2,    2,    2],
             bp =[130,  205,   186,  503,  780],
             cM =[85.7, 149.1, 37.4, 83.6, 134.3],
             maf=[0.5,  0.5,   0.5,  0.5,  0.5])

build_phenome(n_qtl=[2, 3],
              Vg   =[1 .5
                    .5  1])


# Founders --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Option 1: Generate founders by population size n
n        = 10
founders = Founders(n)

# Option 2: Defined by file paths
founders  = Founders("genotypes.csv")
# Option 3: Defined by dataframes
dataframe = CSV.read("genotypes.csv", DataFrame)
founders  = Founders(dataframe)

# Mating --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
args = Dict(:n_shared           => 5,
            :n_per_shared       => 10,
            :replace_shared     => false,
            :replace_per_shared => false,
            :n_per_mate         => 2,
            :ratio_malefemale   => 1)
male, female = mate(sires, dams, args...)

# default
# All mating settings
args = Dict(:nA               => cohort_A.n,
            :nB               => 1,
            :replace_A        => false,
            :replace_B        => false,
            :n_per_mate       => 1,
            :ratio_malefemale => -1)
progenies = mate(cohort_A, cohort_B, args...)
# or the overloaded operators
progenies = cohort_A * cohort_B

# Random mating settings
args = Dict(:nA               => n,
            :nB               => 1,
            :replace_A        => true,
            :replace_B        => true,
            :n_per_mate       => 1,
            :ratio_malefemale => -1)
progenies = mate(cohort_A, cohort_B, args...)

# Selfing mating settings
args = Dict(:nA               => 10,
            :replace_A        => false,
            :n_per_mate       => 50,
            :ratio_malefemale => -1,
            :is_selfing       => true)
progenies = mate(cohort_A, args...)




# Historical founders
# Simulate LD for a base population in linkage and Hardy–Weinberg equilibria
founders = Founders(1500)
for _ in 1:1000
    founders = mate(founders[1:1500])
end
for _ in 1:15
    founders = mate(founders[1:100])
end

# DH
args = Dict(:n_shared           => 5,
            :n_per_shared       => 10,
            :replace_shared     => false,
            :replace_per_shared => false,
            :n_per_mate         => 100,
            :is_DH              => true)
DHs = mate(parents, args...)


# Selection --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
args = Dict(:h2       => [0.5, 0.8],
            :cofactor => cofactors,
            :criteria => "GEBV",
            :equation => "y1 = genotypes
                          y2 = genotypes",
            :method   => "GBLUP")
progenies = select(cohort, 50, args...)

# Breed --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Introduction to Breed()
n_sires        = 10
dams_per_sires = 5
n_dams         = n_sires * dams_per_sires
sires = Founders(n_sires)
dams  = Founders(n_dams)
args  = Dict(:n_per_shared     => dams_per_sires,
             :n_per_mate       => 2,
             :ratio_malefemale => 1.0,
             :is_random        => true,
             :n_gens           => 5,
             :n_select_male    => n_sires)
sires, dams = breed(sires, dams; args...)

# And it's equivalent to the for-loop function:
for _ in 1:5
    males, females = mate(sires, dams, args...)
    sires          = select(males, n_sires, args...)
    dams           = females
end


#Simulate LD from a base population in linkage and Hardy–Weinberg equilibria

# Build Genome and Phenome
build_genome("map.csv")
Vg    = [ 1 .5
         .5  1]
build_phenome("map.csv", Vg)
# Initialize the population with 1,500 founders
founders = Founders(1500)
for _ in 1:1000
    founders = mate(founders[1:1500])
end
for _ in 1:15
    founders = mate(founders[1:100])
end
sires_base = dams_base = founders

# Case --- --- --- --- ---

#Simulate three pure breeds
args_A  = Dict(# Mating
               :n_per_shared     => 10,
               :n_per_mate       => 2,
               :ratio_malefemale => 1,
               # Selection
               :is_random        => true,
               # Breeding
               :n_gens           => 10,
               :n_select_A   => 50)

args_BC = Dict(# Mating
               :n_per_shared     => 20,
               :n_per_mate       => 2,
               :ratio_malefemale => 1,
               # Selection
               :is_random        => true,
               # Breeding
               :n_gens           => 10,
               :n_select_A   => 100)

# Breed A, B, and C
sires_A, dams_A = breed(sires_base, dams_base;  args_A...)
sires_B, dams_B = breed(sires_base, dams_base;  args_BC...)
sires_C, dams_C = breed(sires_base, sires_base; args_BC...)

# NAM
# Base population
common_parents  = Founders(1)
diverse_parents = Founders(25)

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
            :is_random  => true,
            # Breed
            :n_gens     => 6,
            :n_select   => 200)
NAM = Cohort()
for family in F1
    NAM += breed(family, args...)
end


# Overloaded --- --- --- --- ---
cohort_D = cohort_A[1:5] + cohort_B + cohort_C[1:5]

progenies = Cohort()
for animal in cohort
	progenies += mate(animal, args...)
end

# Sample 5 non-duplicated individual from a cohort 
new_cohort = sample(cohort, 5, replace=false)
# By default, sort cohort by their breeding values
sort_cohort = sort(cohort)
# We can sort the cohort based on their ID/pedigree
sort_cohort = sort(cohort; by="ID")


#### Appendix
# haplotypes.csv
1,1,0,0,0,0,1,0
0,0,0,0,1,0,0,0
0,0,0,1,0,0,1,1
1,0,1,0,0,0,1,1
1,1,0,0,1,1,0,0

# genotypes.csv
2,0,0,1
0,0,1,0
0,1,0,2
1,1,0,2
2,0,2,0

# map.csv
id,chr,bp,cM,eff_1,eff_2
snp_1,1,1204,50.8,1.03,31.6
snp_2,1,2938,80.3,-2.34,25.7
snp_3,2,1035,39.2,0.32,30.9
snp_4,2,3653,66.3,0.56,35.1



#
# Alternatively, users can provide the tabular information by \lstinline{DataFrame} instead of the path to the file:
# \begin{lstlisting}
# # Option 3: Defined by dataframes
# using CSV, DataFrames
# dataframe = CSV.read("map.csv", DataFrame)
# build_genome(dataframe)
# # or refer missing genetic positions by literatures
# build_genome(dataframe; species="cattle")
# \end{lstlisting}

# Lastly, manually provide genome information is also supported:
# \begin{lstlisting}
# # Option 4: Defined manually
# ch  = [1,    1,     2,    2,    2]
# bp  = [130,  205,   186,  503,  780]
# cM  = [85.7, 149.1, 37.4, 83.6, 134.3]
# maf = [0.5,  0.5,   0.5,  0.5,  0.5]
# build_genome(ch, bp, cM, maf)
# \end{lstlisting}


# \lstinline{Cohort} can also be treated as a regular iterator in \lstinline{for} loop when ones need, for example, perform mating events with each individuals of the cohort:
# \begin{lstlisting}
# progenies = Cohort()
# for animal in cohort
# 	progenies += mate(animal, args...)
# end
# \end{lstlisting}

# build_genome(nch =[1,    1,     2,    2,    2],
#              bp =[130,  205,   186,  503,  780],
#              cM =[85.7, 149.1, 37.4, 83.6, 134.3],
#              maf=[0.5,  0.5,   0.5,  0.5,  0.5])