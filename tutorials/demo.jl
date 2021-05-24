# ====== Demo1 ====== ====== ====== ====== ====== ====== ======
using XSim
build_demo()

n_sires       = 20
dams_per_sire = 10
n_dams        = n_sires * dams_per_sire

args_mate   = Dict(:n_per_shared       => dams_per_sire,
                   :n_per_mate         => 2,
                   :ratio_malefemale   => 1,
                   :replace_shared     => false,
                   :replace_per_shared => false)
args_select = Dict(:h2                 => [.8, .2],
                   :is_random          => false)
args_breed  = Dict(:n_gens             => 10,
                   :n_select_males     => n_sires)

sires       = Founders(n_sires)
dams        = Founders(n_dams)
sires, dams = breed(sires, dams; args_breed..., args_mate..., args_select...)

# Custom
for i in 1:10
    males, females = mate(sires, dams; args_mate...)
    sires          = select(males, n_sires; args_select...)
    dams           = select(females, n_dams; args_select...)
end

# Results
summary(sires)
summary(dams)
summary(sires + dams)
get_MAF(get_QTLs(sires+dams))

# ====== Demo2 ====== ====== ====== ====== ====== ====== ======
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

