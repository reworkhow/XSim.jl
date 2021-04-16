#One mating produces one child by default.
"""
    children::Cohort = get_children(fathers::Cohort,mothers::Cohort,
                                    nmatings::Int64 = length(fathers.animalCohort)*length(mothers.animalCohort),
                                    nchildren_per_mating::Int64=1;
                                    sex_ratio = [1,1],
                                    fathers_replace::Bool=true, mothers_replace::Bool=true,
                                    strategy="random mating",#"random mating","all matings","embryo transfer")

* Mate a population **fathers** to another population **mothers** with **nmatings** matings
  (one father to one mother) and **nchildren_per_mating** children from each mating
* **sex_ratio**=[nfathers,nmothers], defaulting to [1,1], represents every **nfathers** are mated to **nmothers**.
  For example, every 1 sire is mated to 20 dams in conventional pig breeding (**sex_ratio**=[1,20]);
  1 dam is mated to 10 sires (**sex_ratio**=[20,1]) with embryo transfer; Now we require that at least one element
  in sex_ratio equals 1.
* We have implemented multiple strategies including "random mating" (default), "all matings" (all combinations of parents),
  "embryo transfer" (e.g., sex_ratio=[10,1])
"""
function get_children(fathers::Cohort,
                      mothers::Cohort,
                      nmatings::Int64 = length(fathers.animalCohort)*length(mothers.animalCohort);
                      nchildren_per_mating::Int64=1,
                      sex_ratio = [1,1],
                      fathers_replace::Bool=true, mothers_replace::Bool=true,
                      strategy="all mating" #"random mating","all matings"
                      )

    if strategy == "embryo transfer"
      mothers_replace = false
      nchildren_per_mating = 1
      if sex_ratio[2] != 1
          error("1 dam is mated to multiple sires with embryo transfer, thus the 2nd element of sex_ratio should be 1.")
      end
    end

    #find whether father or mother = 1, thus 1 father (mother) mate to multiple mothers (fathers)
    matingtype = ["onefather_multiplemothers", "onemother_multiplefathers"][findfirst(x->x==1, sex_ratio)]

    nchildren  = nmatings*nchildren_per_mating
    my         = Cohort(Array{Animal}(undef,0),Array{Int64}(undef,0,0))
    resize!(my.animalCohort,nchildren)

    #a dictionary used to count in how many matings an indiviudal is used
    fatherDict = Dict{Int64,Int64}()
    motherDict = Dict{Int64,Int64}()

    println("Get ",nchildren," offspring into the next generation.")


    if strategy == "random mating"
        animali = 1
        for i in 1:nmatings
          if matingtype == "onefather_multiplemothers"
              father = fathers_replace ? get_random_ind(fathers) : get_random_ind_without_replace(fathers,i)
          else
              mother = mothers_replace ? get_random_ind(mothers) : get_random_ind_without_replace(mothers,i)
          end
          for j in 1:sex_ratio[sex_ratio.!=1][1]
              if matingtype == "onefather_multiplemothers"
                  mother = mothers_replace ? get_random_ind(mothers) : get_random_ind_without_replace(mothers,j)
              else
                  father = fathers_replace ? get_random_ind(fathers) : get_random_ind_without_replace(fathers,j)
              end
              for k in 1:nchildren_per_mating
                  my.animalCohort[animali] = get_child(father,mother)
                  animali += 1

                  motherDict[mother.myID] = get(motherDict, mother.myID, 0) + 1
                  fatherDict[father.myID] = get(fatherDict, father.myID, 0) + 1
              end
          end
        end
    elseif strategy == "all matings"
        #nchildren  = length(fathers.animalCohort)*length(mothers.animalCohort)*nchildren_per_mating
        println("Producing ",nchildren," offspring per mating into the next generation.")
        println("There are a total of ",nmatings,"  matings across all parent combinations and thus ", nchildren, "  offspring in the next generation.")
        animali = 1
        for f in fathers.animalCohort #progressMeter
            for m in mothers.animalCohort
                animal = get_child(fathers.animalCohort[f],mothers.animalCohort[m])
                my.animalCohort[animali] = animal
                animali += 1
                motherDict[mother.myID] = get(motherDict, mother.myID, 0) + 1
                fatherDict[father.myID] = get(fatherDict, father.myID, 0) + 1
            end
        end
    end
    println("Number of fathers used: ", length(fatherDict))
    println("Number of mothers used: ", length(motherDict))
    println("Number of offspring generated: ", length(my.animalCohort))
    return my
end

"""
    DHs = get_double_haploids(parents::Cohort, nDHs::Int64)

* Produce *nDHs* double-haploids from a population **parents**.
"""
function get_double_haploids(parents::Cohort, nDHs::Int64)
    println("Producing $nDHs double-haploids from parents randomly selected from a population of size ",length(parents.animalCohort))
    offspring = Cohort(Array{Animal}(undef,0),Array{Int64}(undef,0,0))
    resize!(offspring.animalCohort,nDHs)

    for i in 1:nDHs
        parent = get_random_ind(parents)
        offspring.animalCohort[i] = get_double_haploid(parent)
    end
    return offspring
end

# Get a subset (1 individual) from a population
function get_random_ind(my::Cohort)
    cohort_size = length(my.animalCohort)
    thisone     = rand(1:cohort_size)
    return my.animalCohort[thisone]
end

# Get one individual from a population without replacement (need to be used in a for loop)
function get_random_ind_without_replace(my::Cohort,i::Int64)
   cohort_size=length(my.animalCohort)
   if i >= cohort_size
     println("No animal left in cohort to sample.")
     exit(1)
   end
   sampled_one                  = rand(i:cohort_size)
   current_ith_ind              = my.animalCohort[i]
   my.animalCohort[i]           = my.animalCohort[sampled_one] #move the sampled ind to the i_th element
   my.animalCohort[sampled_one] = current_ith_ind
   return my.animalCohort[i]
end
