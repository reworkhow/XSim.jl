using XSim

build_genome(n_chr=1, n_marker=50)

data = DATA("map")

### write down: MAF is optional 
### assign additional markers in the same segment  for reference genome

build_genome(data)
build_phenome(data)

build_genome(species="cattle")

### Provide QTL effects
### Save updated map to users' folders (summary)
build_phenome([5, 10], 
             h2=[0.3, .9], 
             vg=[1 .5
                 .5 1])

GLOBAL("effects_QTLs")
GLOBAL("effects")

DATA("map")


# Cohort
data_g = DATA("genotypes")

cohort = Cohort(data_g, n=3, random=false)

get_genotypes(cohort)


# Setup --- --- --- --- --- --- --- --- --- --- --- --- --- ---
using XSim
QTL_effects = [1.0   0
               0   1.0
               0   1.0
               1.0   0]
QTL_freq = [0.5 for i in 1:4]
Vg_goal = [1 0.5; 0.5 1]

# Scaled effects (XSim function) --- --- --- --- --- --- --- ---
XSim.scale_effects(QTL_effects, QTL_freq, Vg_goal)

# Scaled effects (line by line) --- --- --- --- --- --- --- ---
# Compute Vg for input QTL_effects
Vg_ori    = XSim.get_Vg(QTL_effects, QTL_freq)

# Decompose original variance
Vg_ori_U  = XSim.cholesky(Vg_ori).U'
Vg_ori_Ui = inv(Vg_ori_U)

# Decompose goal variance
Vg_goal_U = XSim.cholesky(Vg_goal).U

# m by t = m by t * t by t
QTL_effects * Vg_ori_Ui'Vg_goal_U


