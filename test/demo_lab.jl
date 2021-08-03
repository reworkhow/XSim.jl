using XSim

build_demo()

cohort = Founders(10)




XSim.get_genotypes(cohort)


get_BVs(cohort)
get_QTLs(cohort)


sires = Cohort(5)
dams  = Cohort(20)

args_mate = Dict(:nA           => 3,
                :nB_per_A     => 5,
                :n_per_mate   => 2,
                :ratio_malefemale => 1.0)

males, females =  

mate(sires, dams; args_mate...)
mate(sires, dams, nA = 3, nB_per_A = 5, n_per_mate = 2, ratio_malefemale = 1.0)

pool = males + females
pool = males[1:5] + females[5:10]
pool |> get_pedigree


build_demo()
build_phenome([50, 50])
males = Founders(500)
args_case1  = Dict(:h2     => [0.99, .5],
                   :weights=> [.5, .5])
ind = select(males, 50; args_case1...)





args_case2  = Dict(:h2     => [.8, .5],
                   :weights=> [1.0, 0.0])
sel_case2 = select(males, 3; args_case2...)

args_case3  = Dict(:h2     => [.8, .5],
                   :weights=> [0.0, 1.0])
sel_case3 = select(males, 3; args_case2...)

sel_case2 # top-3 on trait_1 
sel_case3 # top-3 on trait_2
sel = sel_case2 + sel_case3  #

summary(sel)