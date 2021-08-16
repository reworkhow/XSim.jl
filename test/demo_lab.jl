using XSim

build_demo_small()

cohort_A = Cohort(50)    # 5 parents
cohort_B = Cohort(150)   # 30 dams


args_mating = Dict(:nA         => 5,
                   :nB_per_A   => 30,
                   :n_per_mate => 1)

offsprings = mate(cohort_A, cohort_B; args_mating...)


mate(cohort_A         ::Cohort,
     cohort_B         ::Cohort;
     nA               ::Int64=cohort_A.n,
     nB_per_A         ::Int64=1,
     n_per_mate       ::Int64=1,
     replace_A        ::Bool =false,
     replace_B        ::Bool =false,
     ratio_malefemale ::Union{Float64, Int64}=0,
     scheme           ::String ="none",
     args...)


## Bug list
# 1. index for cohort_A
