#select animals from a pool of candidates
function select(n::Int, animals::Cohort;
                weights=false,
                criteria="phenotypes", selection_criteria="positive") #"ebv"...#"negetive"
    if weights == false
        weights = ones(size(common.LRes, 2))
    end

    if selection_criteria == "positive"
        direction = 1
    end

    candidates   = deepcopy(animals)
    selected_animals = Cohort(Array{Animal}(undef, 0),
                              Array{Int64}(undef, 0, 0))

    y = getOurPhenVals(candidates) * weights * direction
    # select top-n animals
    selected_animals.animalCohort = candidates.animalCohort[sortperm(y)][(end - n + 1):end]
    return selected_animals
end

#male female ratio (e.g.,1 male to 10 females)
#mate sires and dams to generate offsprings
function mating(n::Int, sires::Cohort, dams::Cohort; strategy="random")
    boys = get_children(sires, dams, round(Int, n / 2))
    gals = get_children(sires, dams, round(Int, n / 2))
    return boys, gals
end


function mating(n::Int, plants::Cohort; strategy="random")
    boys, gals = mating(n, plants, plants, strategy=strategy)
    return concatCohorts(boys, gals)
end

function selection_for_ngenerations(noffspring, nsires, ndams, sires, dams;
                                    ngenerations=5, strategy="random")
    maleCandidates   = deepcopy(sires)
    femaleCandidates = deepcopy(dams)
    boys  = Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))
    gals  = Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))

    for i = 1:ngenerations
        selected_sires = select(nsires, sires, criteria="phenotypes",
                                selection_criteria="positive")
        selected_dams  = select(ndams, dams, criteria="phenotypes",
                                selection_criteria="positive")
        boys, gals     = mating(noffspring, selected_sires, selected_dams,
                                strategy="random")
        maleCandidates.animalCohort   = [sires.animalCohort; boys.animalCohort]
        femaleCandidates.animalCohort = [dams.animalCohort;  gals.animalCohort]
    end
    return boys, gals
end

##########
#genetic_evaluation(sires,dams)
########




