mutable struct Cohort
    animals      ::Array{Animal, 1}
    n            ::Int64

    ```Initiate cohorts with sampled genotypes```
    function Cohort(n          ::Int64=0)

        cohort = Array{Animal}(undef, n)
        for i in 1:n
            cohort[i] = Animal(Animal(), Animal())
        end

        return new(cohort, n)
    end

    ```Initiate cohorts with given genotypes```
    function Cohort(genotypes   ::Union{DataFrame, Array{Int64}};
                    n           ::Int64=-1,
                    alter_maf   ::Bool=true)

        if isa(genotypes, DataFrame)
            n_founders = nrow(genotypes)
            n_loci     = ncol(genotypes)
            genotypes = Array(genotypes)

        elseif isa(genotypes, Array)
            n_founders = size(genotypes, 1)
            n_loci     = size(genotypes, 2)

        else
            LOG("Input genotypes not supported", "error")
        end

        if n_loci != GLOBAL("n_loci")
            LOG("Number of loci mismatches the given genome", "error")
        end

        if alter_maf
            LOG("MAF has been updated based on provided haplotypes/genotypes")
            SET("maf", get_MAF(genotypes))
        end

        if n == -1
            n = n_founders
        end

        cohort        = Array{Animal}(undef, n)
        for i in 1:n
            hap       = vcat(genotypes[i, :]...)
            cohort[i] = Animal(Animal(), Animal(), haplotypes=hap)
        end

        return new(cohort, n)
    end

    function Cohort(filename    ::String; args...)
        dt = CSV.read(filename, DataFrame,
                       header=false,
                       missingstrings=["-1", "9"])

        return Cohort(dt; args...)
    end

    # Constructor for non-founder
    function Cohort(animals     ::Array{Animal, 1})
        n = length(animals)
        return new(animals, n)
    end

    function Cohort(animal      ::Animal)
        return new([animal], 1)
    end
end

Founders(genotypes::Union{DataFrame, Array{Int64}}; args...) = Cohort(genotypes; args...)
Founders(filename ::String; args...)                         = Cohort(filename; args...)
Founders(n        ::Int64)                                   = Cohort(n)

function Base.summary(cohort::Cohort; is_return=true)
    bvs        = get_BVs(cohort)
    mu_g       = round.(XSim.mean(bvs, dims=1), digits=3)
    var_g      = round.(XSim.var(bvs,  dims=1), digits=3)

    if is_return
        return Dict("n"     => cohort.n,
                    "mu_g"  => mu_g,
                    "var_g" => var_g)
    else
        LOG("Cohort ($(cohort.n) individuals)")
        LOG()
        LOG("Mean of breeding values: ")
        LOG("$mu_g")
        LOG()
        LOG("Variance of breeding values: ")
        LOG("$var_g")
    end
end


function genetic_evaluation(cohort         ::Cohort;
                            cofactors      ::DataFrame=DataFrame(),
                            model_equation ::String="",
                            covariate      ::String="",
                            random_iid     ::String="",
                            random_str     ::String="")

    jwas_ped = get_pedigree(cohort,   "JWAS")
    jwas_P   = get_phenotypes(cohort, "JWAS"; cofactors=cofactors)
    genotypes= get_genotypes(cohort,  "JWAS") # 0 1 2


    # Step 3: Build Model Equations
    if model_equation == ""
        traits = names(jwas_P)
        array_eq = ["$(trait) = intercept + genotypes" for trait in traits]
        model_equation = join(array_eq, "\n")
    end

    model = JWAS.build_model(model_equation);

    # Step 4: Set Factors or Covariates
    if covariate != ""
        JWAS.set_covariate(model, covariate);
    end

    # Step 5: Set Random or Fixed Effects
    if random_iid != ""
        JWAS.set_random(model, random_iid);
    end

    if random_str != ""
        JWAS.set_random(model, random_str, jwas_ped);
    end

    # Step 6: Run Analysis
    out = JWAS.runMCMC(model, jwas_P);

    return out
end



get_QTLs(cohort::Cohort) = get_genotypes(cohort)[:, GLOBAL("is_QTLs")]
get_IDs(cohort::Cohort)  = (animal->animal.ID).(cohort)

function get_genotypes(cohort::Cohort, option::String="XSim")
    genotypes_2d_tmp = (animal->get_genotypes(animal)).(cohort)
    genotypes        = hcat(genotypes_2d_tmp...)'

    if option == "XSim"
        return genotypes

    elseif option == "JWAS"
        dt_G = hcat(get_IDs(cohort), genotypes) |> XSim.DataFrame
        # dt_G = hcat("a".* string.(get_IDs(cohort)), genotypes) |> XSim.DataFrame
        CSV.write("jwas_g.csv", dt_G)
        return JWAS.get_genotypes("jwas_g.csv")
    end
end

function get_BVs(cohort::Cohort)
    bv_2d = (animal->get_BVs(animal)).(cohort)
    return hcat(bv_2d...)'
end

# available types: phenotypic, genotypic, estimated
function get_phenotypes(cohort   ::Cohort,
                        option   ::String="XSim";
                        cofactors::DataFrame=DataFrame(),
                        h2       ::Union{Array{Float64}, Float64}=.5,
                        Ve       ::Union{Array{Float64}, Float64}=get_Ve(GLOBAL("n_traits"),
                                                                         GLOBAL("Vg"),
                                                                         h2),
                        return_Ve::Bool=false)

    phenotypes_tmp = (animal->get_phenotypes(animal; h2=h2, Ve=Ve)).(cohort)
    phenotypes     = hcat(phenotypes_tmp...)'
    Ve             = handle_diagonal(Ve, GLOBAL("n_traits"))

    if option == "XSim"
        if return_Ve
            return phenotypes, Ve
        else
            return phenotypes
        end

    elseif option == "JWAS"
        jwas_P = hcat(get_IDs(cohort), phenotypes) |> XSim.DataFrame
        rename!(jwas_P, vcat(["ID"], ["y$i" for i in 1:GLOBAL("n_traits")]))
        if nrow(cofactors) != 0
            jwas_P = hcat(jwas_P, cofactors)
        end
        return jwas_P
    end
end

function get_pedigree(cohort::Cohort, option::String="XSim")

    # return a 3-column matrix: ID, SireID, DamID
    ped_tmp  = (animal->[animal.ID, animal.sire.ID, animal.dam.ID]).(cohort)
    ped_array = hcat(ped_tmp...)'

    if option == "XSim"
        return ped_array

    elseif option == "JWAS"
        # Get pedigree in dataframe
        df = ped_array |> Array |> XSim.DataFrame

        # Cast columns to strings
        df[!,1] = strip.(string.(df[!,1]))
        df[!,2] = strip.(string.(df[!,2]))
        df[!,3] = strip.(string.(df[!,3]))

        # Instantiate Pedigree
        ped = JWAS.PedModule.Pedigree(
                1,
                Dict{AbstractString, JWAS.PedModule.PedNode}(),
                Dict{Int64, Float64}(),
                Set(), Set(), Set(), Set(),
                Array{String, 1}())

        # JWAS things
        JWAS.PedModule.fillMap!(ped,df)
        for id in keys(ped.idMap)
         JWAS.PedModule.code!(ped, id)
        end

        for id in keys(ped.idMap)
          JWAS.PedModule.calcInbreeding!(ped,id)
        end

        ped.IDs = JWAS.PedModule.getIDs(ped)

        JWAS.PedModule.get_info(ped)
        JWAS.PedModule.writedlm("IDs_for_individuals_with_pedigree.txt",ped.IDs)

        return ped
    end
end


function get_DH(parents::Cohort, n::Int64)
    animals = Array{Animal}(undef, n)
    select_idx = sample(parents.n, n, replace=true)
    for i in 1:nDHs
        parent = parents[select_idx[i]]
        animals[i] = get_DH(parent)
    end
    return cohort(animals)
end

function get_haplotype_founder(cohort::Cohort)
    haps = []
    for animal in cohort
        for chromosome in [animal.genome_sire; animal.genome_dam]
            for ori in chromosome.ori
                push!(haps, ori)
            end
        end
    end
    haps
end

# ped, mme, out is from JWAS get_pedigree(), build_model(), and solve()
function putEBV(cohort::Cohort, ped, mme, out)
    # transfer ebv from mme to XSim
    trmAnimal = mme.modelTermDict["1:Animal"]
    for animal in cohort
        id = animal.ID
        strID = string(id)
        mmePos = ped.idMap[strID].seqID + trmAnimal.startPos - 1
        animal.traits[1].estimated = out[mmePos, 2]
    end
end

function sample(cohort ::Cohort,
                n      ::Int64;
                replace::Bool=true)

    select = sample(1:cohort.n, n, replace=replace)
    return cohort[select]
end

function print(cohort::Cohort, option::String="None")
    if option == "None"
        if cohort.n != 0
            Base.summary(cohort, is_return=false)
        end

    elseif option == "ID"
        print("Individual: [ ")
        for animal in cohort
             print(animal.ID, " ")
        end
        println("]")

    elseif option == "Pedigree"
        return get_pedigree(cohort)
    end
end

function getindex(cohort::Cohort, I...)
    if length(I...) == 1
        return cohort.animals[I[1]]
    else
        return Cohort(getindex(cohort.animals, I...))
    end
end

Base.:+(x::Cohort, y::Cohort)       = Cohort(vcat(x.animals, y.animals))
Base.:+(x::Cohort, y::Animal)       = Cohort(vcat(x.animals, y))
Base.:+(x::Animal, y::Cohort)       = Cohort(vcat(x, y.animals))
Base.length(cohort::Cohort)         = length(cohort.animals)
Base.show(io::IO, cohort::Cohort)   = GLOBAL("silent") ? nothing : print(cohort)
Base.iterate(cohort::Cohort, i...)  = Base.iterate(cohort.animals, i...)
Base.lastindex(cohort::Cohort)      = length(cohort)

Base.setindex!(cohort::Cohort, animal::Animal, i::Int64) =
    Base.setindex!(cohort.animals, animal, i)




# function get(cohort::Cohort,
#              item  ::String,
#              option::Any)

#     if item == "Traits"
#         return get_phenotypes(cohort, option)

#     elseif item == "ID"
#         return get_IDs(cohort)

#     elseif item == "Pedigree"
#         return get_pedigree(cohort)

#     elseif item == "DH"
#         return get_DH(cohort, option)

#     else
#         println("""
#             The available options are: 'Triats', 'ID', 'Pedigree', and 'DH'
#         """)
#     end
# end/