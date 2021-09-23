using XSim
using Test

@testset "ReadMe" begin
    CLEAR()
    SILENT(true)

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
end

@testset "Overloaded" begin
    CLEAR()
    SILENT(true)
    build_demo_small()

    #concatenate multiple cohorts
    cohort_A = Cohort(10)
    cohort_B = Cohort(5)
    cohort_C = Cohort(3)
    cohort_D = cohort_A[1:5] + cohort_B + cohort_C

    new_cohort = sample(cohort_A, 5, replace=false)
    # sort cohort by true breeding values (BV). "BV" is the default value.
    sort_cohort = sort(cohort_A, by="BV")
    # or sort the cohort based on their pedigree in an order from the oldest to the youngest. Other options include estimated breeding values (EBV) or phenotypes (e.g., y1). 
    sort_cohort = sort(cohort_A, by="pedigree")
end

@testset "Breed" begin
    CLEAR()
    SILENT(true)
    build_demo_small()
    # Mating and selection cross 5 generations
    args  = Dict(# mating arguments
                :nA               => 10,
                :nB_per_A         => 5,
                :n_per_mate       => 1,
                :ratio_malefemale => 1.0,
                # selection arguments
                :criteria         => "EBV",
                :methods          => "GBLUP",
                # breeding arguments
                :n_gens           => 5,
                :nA_select        => 10)
    # Simulate 10 sires and 50 dams as founders
    sires = Founders(10)
    dams  = Founders(50)
    # Breed cohorts based on the defined arguments 
    sires, dams = breed(sires, dams; args...)
    # # The result is equivalent to the following mate-select iterations:
    # for _ in 1:5
    #     males, females = mate(sires, dams; args...)
    #     sires          = select(males, sires.n; args...)
    #     dams           = females
    # end
end

@testset "JWAS" begin
    CLEAR()
    SILENT(true)
    build_demo_small()
    # Genetic evaluation
    cohort = Cohort(100)
    args = Dict(:criteria => "EBV",
                :methods  => "GBLUP",
                :weights  => [3.0, -2.0])
    offspring = select(cohort, 50; args...)

    genetic_evaluation(cohort);
    summary()
end

@testset "Mating" begin
    CLEAR()
    SILENT(true)
    build_demo_small()

    # Mating
    sires = Cohort(5)
    dams  = Cohort(5 * 10)
    args  = Dict(:nA                 => 5,
                 :nB_per_A           => 10,
                 :replace_A          => false,
                 :replace_B          => false,
                 :n_per_mate         => 1,
                 :ratio_malefemale   => 1)
    male, female = mate(sires, dams; args...)

    cohort_A = Cohort(10)
    cohort_B = Cohort(10)
    args = Dict(:nA               => cohort_A.n,
                :nB_per_A         => 1,
                :replace_A        => false,
                :replace_B        => false,
                :n_per_mate       => 1)
    progenies = mate(cohort_A, cohort_B; args...)
    # Equivalent results without providing any argument
    progenies = mate(cohort_A, cohort_B)
    # Equivalent results by specifying scheme argument
    progenies = mate(cohort_A, cohort_B, scheme = "random")
    # Equivalent results with overloaded operator '*'
    progenies = cohort_A * cohort_B

    # Diallel cross mating scheme
    args = Dict(:nA         => cohort_A.n,
                :nB_per_A   => cohort_B.n,
                :replace_A  => false,
                :replace_B  => false,
                :n_per_mate => 1)
    offspring = mate(cohort_A, cohort_B; args...)
    # Equivalent results by specifying scheme argument
    offspring = mate(cohort_A, cohort_B, 
                    scheme = "diallel cross")
    get_pedigree(offspring)

    # Selfing mating scheme
    args = Dict(:nA         => 10,
                :replace_A  => false,
                :n_per_mate => 50,
                :scheme     => "selfing")
    offspring = mate(cohort_A; args...)
    get_pedigree(offspring)

    DHs = get_DH(cohort_A)

    # Pedigree mating scheme
    founders = mate(DATA("pedigree"))
end


@testset "Load by files or datafrmaes" begin
    # load by file or data frames
    path = PATH("map")
    build_genome(path)
    build_phenome(path, h2 = 0.3)
    build_genome(DATA("map"))
    build_phenome(DATA("map"), h2 = 0.3)
    founders = Founders(10)
    founders = Founders(PATH("haplotypes"))
    founders = Founders(PATH("genotypes"))

    # A quick start of genome with 10k markers on 2 chromosomes __
    build_genome(n_chr = 2,
                n_marker = 100)
    # Phenome with 2 correlated traits controlled by 100 QTLs
    build_phenome([3, 5],
                vg = [1 .5
                    .5  1],
                h2 = [0.3, 0.8])
end


@testset "NAM" begin
    # clear globals
    CLEAR()
    SILENT(true)

    # set up genome and phenome
    build_demo_small()

    # Initialize a population with 1,500 founders
    founders = Founders(1500)
    # Let founders random mate with each other
    # for 15 generations
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
                :n_per_mate => 10,
                :scheme     => "selfing",
                # Selection
                :criteria   => "phenotypes",
                # Breed
                :n_gens     => 4,
                # single-seed decent
                :n_select   => 1)
    NAM = Cohort()
    for family in F1
        F2 = mate(family, n_per_mate=200, scheme="selfing")
        for seed in F2
            NAM += breed(seed; args...)
        end
    end
end

@testset "Rotational cross breed" begin
    # clear globals
    CLEAR()
    SILENT(true)

    # set up genome and phenome
    build_demo_small()

    # Initialize a population with 1,500 founders
    founders = Founders(1500)
    # Let founders random mate with each other
    # for 15 generations
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
    args_XA =  Dict(:nA               => 50,
                    :nB_per_A         => 20,
                    :n_per_mate       => 2,
                    :ratio_malefemale => 1)
    args_XBC = Dict(:nA               => 100,
                    :nB_per_A         => 10,
                    :n_per_mate       => 2,
                    :ratio_malefemale => 1)

    args_A[:n_gens]  = 1
    args_BC[:n_gens] = 1
    # Rotation (G1)
    sires_A1, dams_A1    = breed(sires_A, dams_A; args_A...)
    sires_B1, dams_B1    = breed(sires_B, dams_B; args_BC...)
    sires_C1, dams_C1    = breed(sires_C, dams_C; args_BC...)
    males_G1, females_G1 = mate(sires_B, dams_C;
                                args_XBC...)
    # Rotation (G2)
    sires_A2, dams_A2    = breed(sires_A1, dams_A1; args_A...)
    sires_B2, dams_B2    = breed(sires_B1, dams_B1; args_BC...)
    sires_C2, dams_C2    = breed(sires_C1, dams_C1; args_BC...)
    males_G2, females_G2 = mate(sires_A1, females_G1;
                                args_XA...)
    # Rotation (G3)
    sires_A3, dams_A3    = breed(sires_A2, dams_A2; args_A...)
    sires_B3, dams_B3    = breed(sires_B2, dams_B2; args_BC...)
    sires_C3, dams_C3    = breed(sires_C2, dams_C2; args_BC...)
    males_G3, females_G3 = mate(sires_C2, females_G2;
                                args_XBC...)
end

