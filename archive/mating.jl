function sampleAllMatingsSel(numOffPerMating, nSires::Int64, nDams::Int64,
    nGen::Int64, maleParents, femaleParents;
    gen=1,fileName="", weights = false, direction=1)

    maleCandidates   = deepcopy(maleParents)
    femaleCandidates = deepcopy(femaleParents)
    sires = Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))
    dams = Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))
    offspring  = Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))
    boys  = Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))
    gals  = Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))

    if weights == false
        weights = ones(size(GLOBAL.LRes,2))
    end
    println(weights)

    for i=1:nGen
        # @printf "Generation %5d: sampling %5d offspring per mating by crossing %5d male parents to each of %5d female parents\n" gen+i numOffPerMating nSires nDams
        y = getOurPhenVals(maleCandidates)*weights*direction
        sires.animals = maleCandidates.animals[sortperm(y)][(end-nSires+1):end]
        # @printf "Phenotypically best %5d fathers selected from cohort of male parents of size %5d\n" nSires length(maleCandidates.animals)
        y = getOurPhenVals(femaleCandidates)*weights*direction
        dams.animals = femaleCandidates.animals[sortperm(y)][(end-nDams+1):end]
        # @printf "Phenotypically best %5d mothers selected from cohort of female parents of size %5d\n" nDams length(femaleCandidates.animals)
        offspring = sampleOffAllMatings(sires,dams,numOffPerMating)
        # @printf "Dividing offspring into half males and females by sampling %5d males randomly from %5d offspring\n" round(Int,length(offspring.animals)/2) length(offspring.animals)
        numMaleOff = round(Int,length(offspring.animals)/2)
        #offspring = getRandomSampleOfIndWoReplacement!(offspring, numMaleOff)
        boys.animals = offspring.animals[1:numMaleOff]
        gals.animals = offspring.animals[numMaleOff+1:end]
        @printf "boys and gals done\n"
        # if fileName!=""
        #     outputPedigree(boys,fileName)
        #     outputPedigree(gals,fileName)
        # end
        maleCandidates.animals   = [sires.animals; boys.animals]
        femaleCandidates.animals = [dams.animals;  gals.animals]
    end
    gen += nGen
    return boys,gals, gen
end


function sampleBLUPSel(popSize::Int64, nSires::Int64, nDams::Int64, nGen::Int64,maleParents,femaleParents,varRes=GLOBAL.varRes,varGen=1;gen=1,fileName="XSim", direction=1)
    error("sampleBLUPSel() not adapted to multitrait selection yet")
    # initial BLUP evaluation--- parents
    run(`\rm -f $fileName.ped`)
    run(`\rm -f $fileName.phe`)
    run(`\rm -f $fileName.brc`)
    run(`\rm -f $fileName.gen`)
    outputPedigree(maleParents,fileName)
    outputPedigree(femaleParents,fileName)
    pedfile   = fileName*".ped"
    phenofile = fileName*".phe"
    colNames = ["Animal";"y";"bv"]
    dfPhen = CSV.read(phenofile, DataFrame, delim = ' ',header=false,names=colNames)
    ped = get_pedigree(pedfile)
    mme = build_model("y = intercept + Animal",varRes)
    set_random(mme,"Animal",ped,varGen)
    out = solve(mme,dfPhen,solver="GaussSeidel",printout_frequency=40)
    # trasfer BLUP-EBV to animals
    putEBV(maleParents,ped,mme,out)
    putEBV(femaleParents,ped,mme,out)

    maleCandidates   = copy(maleParents)
    femaleCandidates = copy(femaleParents)
    sires = Cohort(Array{Animal}(undef,0), Array{Int64}(undef,0,0))
    dams  = Cohort(Array{Animal}(undef,0), Array{Int64}(undef,0,0))
    boys  = Cohort(Array{Animal}(undef,0), Array{Int64}(undef,0,0))
    gals  = Cohort(Array{Animal}(undef,0), Array{Int64}(undef,0,0))
    for i = 1:nGen
        @printf "Generation %5d: sampling %5d males and %5d females\n" gen+i round(Int,popSize/2) round(Int,popSize/2)
        y = direction*[animal.ebv for animal in maleCandidates.animals]
        sires.animals = maleCandidates.animals[sortperm(y)][(end-nSires+1):end]
        y = direction*[animal.ebv for animal in femaleCandidates.animals]
        dams.animals = femaleCandidates.animals[sortperm(y)][(end-nDams+1):end]
        boys = get_children(sires, dams, round(Int, popSize / 2))
        gals = get_children(sires, dams, round(Int, popSize / 2))
        outputPedigree(boys,fileName)
        outputPedigree(gals,fileName)
        maleCandidates.animals   = [sires.animals; boys.animals]
        femaleCandidates.animals = [dams.animals;  gals.animals]

        # BLUP ebvs
        dfPhen = CSV.read(phenofile,DataFrame,delim = ' ',header=false,names=colNames)
        ped = get_pedigree(pedfile)
        mme = build_model("y = intercept + Animal",varRes)
        set_random(mme,"Animal",ped,varGen)
        out = solve(mme,dfPhen,solver="GaussSeidel",printout_frequency=40)
        # trasfer BLUP-EBV to animals
        putEBV(maleCandidates,ped,mme,out)
        putEBV(femaleCandidates,ped,mme,out)
    end
    gen += nGen
    return boys,gals, gen
end
