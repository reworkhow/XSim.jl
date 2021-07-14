using  XSim
using Statistics

n_replicates = 1000
BVs    = zeros(5, n_replicates)
for i in 1:n_replicates
    # Setup genome and phenome
    CLEAR()
    SILENT(true)
    build_demo()
    # Initialize founders
    cohorts = Founders(2)
    # Reproduce progenies
    a3 = Animal(cohorts[1], cohorts[2])
    a4 = Animal(cohorts[1], cohorts[2])
    a5 = Animal(a3, a4)
    # Store breeding values
    BVs[1, i] = get_BVs(cohorts[1])[1]
    BVs[2, i] = get_BVs(cohorts[2])[1]
    BVs[3, i] = get_BVs(a3)[1]
    BVs[4, i] = get_BVs(a4)[1]
    BVs[5, i] = get_BVs(a5)[1]
end

# Summarize statistics for BVs
var(BVs, dims=2)

# Inspect pedigree
final_pop = cohorts + a3 + a4 + a5
get_pedigree(final_pop)

get_QTLs(final_pop)
get_QTLs(final_pop) |> get_MAF

get_QTLs(cohorts) |> get_MAF
get_QTLs(a3 + a4 + a5) |> get_MAF
