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
    selected_animals = Cohort()

    y = getOurPhenVals(candidates) * weights * direction
    # select top-n animals
    selected_animals.animals = candidates.animals[sortperm(y)][(end - n + 1):end]
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
    boys  = Cohort()
    gals  = Cohort()

    for i = 1:ngenerations
        selected_sires = select(nsires, sires, criteria="phenotypes",
                                selection_criteria="positive")
        selected_dams  = select(ndams, dams, criteria="phenotypes",
                                selection_criteria="positive")
        boys, gals     = mating(noffspring, selected_sires, selected_dams,
                                strategy="random")
        maleCandidates.animals   = [sires.animals; boys.animals]
        femaleCandidates.animals = [dams.animals;  gals.animals]
    end
    return boys, gals
end

##########
#genetic_evaluation(sires,dams)
########




