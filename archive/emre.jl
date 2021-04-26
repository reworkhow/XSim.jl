repl = isempty(ARGS) ? "arg1" : ARGS[1]

#repl = isdefined(:newARGS) ? newARGS : ARGS
println(typeof(parse(Int32,repl[1])))

using DataFrames
using DelimitedFiles
using Distributions
using XSim
using Random
using CSV

thisBreed       = "JER"
nFounderIndJER  = 1305
nFounderXSimJER = 1305 #will be 1305 for haplotypes or 1.000 or 10.000 for frequencies
nFounderSire    = 100
nFounderDam     = nFounderXSimJER-nFounderSire
nGen            = 10
thisChr      = 24
nFounderLoci = 487880
regionSize   = 0.25 #used to select a region from the chromosome
nQTL         = 1000

#repl = 1
Random.seed!(parse(Int,"123$repl"))

path2GenoData = "/usr/home/qgg/emre/zexiData/"
outputFolder  = path2GenoData*"mySimData/$thisBreed/repl$repl"

if isdir(outputFolder)==true
    println("$outputFolder exists. Removing its content")
    run(`rm -rf $outputFolder/*`)
else
    println("$outputFolder has been created")
    run(`mkdir -p $outputFolder`)
end

cd(outputFolder)

println("Working directoty set to: $outputFolder")

finalData = Array{Int64}(undef, 0, nFounderLoci)
f = open(path2GenoData*"selected$(thisBreed)Ind_chr$(thisChr).ped")
for i in 1:nFounderIndJER
    print(i,',')
    x = readline(f)
    x = split(x)
    x = hcat(x...)
    x = parse.(Int64,x[7:end])
#    x = DataFrames.recode(x[7:end], "A"=>1,"C"=>2,"G"=>3,"T"=>4)
    global finalData  = vcat(finalData ,reshape(x, 2, :))
end
eof(f)
finalData .-= 1;

mapData = readdlm(path2GenoData*"recoded_chr$(thisChr).map",header=false,'\t')
length_cm = maximum(mapData[:,4])/(10^6)
length_morgan = length_cm / 100
mapInMorgan = round.(mapData[:,4]./(10^8),digits=8)
mapData = hcat(mapData,mapInMorgan)

selected01 = []
selCount = 0
for i in 1:1000
    midPos = rand(Uniform(0,length_morgan))
    firstPos,lastPos = midPos-(regionSize/2),midPos+(regionSize/2)
    if (firstPos > 0) & (lastPos <= length_morgan)
    global selected01 = push!(selected01,[firstPos midPos lastPos])
        global selCount += 1
    end
end
if isempty(selected01) == true
    selected01 = [minimum(mapData[:,5]) minimum(mapData[:,5])+(regionSize*0.5) minimum(mapData[:,5])+regionSize]
    println("Default region size was used: $selected01")
    else selected01 = sample(selected01,1)
    println("Selected region from $selCount available: $selected01")
end

selectedLoci = findall((selected01[][1] .<= mapData[:,5]) .& (selected01[][3] .>= mapData[:,5]))

mapData   = mapData[selectedLoci,:]
CSV.write("usedMap_$thisChr", convert(DataFrame,mapData))

haploIDs = vcat([["p$i"; "m$i"] for i in 1:nFounderIndJER]...)

finalData = finalData[:,selectedLoci];
writedlm("haploData4XSim_$(thisBreed)_chr$(thisChr).txt",[haploIDs finalData])

alleleFreq = vec(mean(finalData,dims=1)./2)
CSV.write("usedAlleleFreq_$(thisBreed)$(thisChr)", convert(DataFrame,[collect(1:length(alleleFreq)) alleleFreq]))

numChr,numLoci,chrLength,mapPos,mutRate = 1,size(mapData,1),length_morgan,mapInMorgan,0.0
build_genome(numChr,chrLength+0.0000000000000001,numLoci,alleleFreq,mapPos,mutRate)

nTotInd      = nFounderXSimJER
baseAnimals  = Cohort(nFounderXSimJER,"haploData4XSim_$(thisBreed)_chr$(thisChr).txt");

baseSireID   = sample(1:nFounderXSimJER,nFounderSire,replace=false)
baseDamID    = setdiff(1:nFounderXSimJER,baseSireID)
baseSires = cohortSubset(baseAnimals,baseSireID)
baseDams  = cohortSubset(baseAnimals,baseDamID);
outputPedigree(baseAnimals, "output$(thisBreed)_$(thisChr).txt")

generationCounter = Int.(zeros(nFounderXSimJER));

currentCohort = deepcopy(baseAnimals)
nLastInd      = deepcopy(nFounderXSimJER)
generationCounterA = deepcopy(generationCounter) #for the hist ind data

    nBreedAMale    = nFounderSire
    nBreedAFemale  = nFounderDam
    nAOffspring    = 2
    nBreedA        = nBreedAMale + nBreedAFemale
    nGenA          = nGen

for gen in 1:nGenA

    lastGenInd     = (nTotInd-nLastInd+1):nTotInd
    println(lastGenInd)
    breedASireID   = sample(lastGenInd,nBreedAMale,replace=false)
    breedANonSire  = setdiff(lastGenInd,breedASireID)
    breedADamID    = sample(breedANonSire,nBreedAFemale,replace=false)

    myBreedASire = repeat(breedASireID,outer=[Int(ceil(nBreedA/nBreedAMale))])[1:nBreedAFemale,:]
    myBreedADam  = repeat(breedADamID,outer=[Int(ceil(nBreedA/nBreedAFemale))])[1:nBreedAFemale,:]
    myBreedAPed  = [myBreedASire myBreedADam]; ##not pedigreee actually... just a mating because some matings produce 2 offsprings

    ######## added to remove some animals. Each meting produces 2 offspring, but I remove some to keep popSize constant
    # removed animals have ID assigned, but simply removed from the pedigree. So I have some IDs missing in ped, which is OK.
    keepNoff = shuffle([fill(1,nBreedAFemale-nBreedAMale);fill(2,nBreedAMale)]) #I use nBreedAMale just as a number, no consequence
    for i in 1:size(myBreedAPed,1)
        global    nowSires,nowDams,nowGen = sampleRan(nAOffspring, 1, cohortSubset(currentCohort,findall(x->x==myBreedAPed[i,1],lastGenInd)), cohortSubset(currentCohort,findall(x->x==myBreedAPed[i,2],lastGenInd))) ##always one more generation!!!
        if keepNoff[i] == 1
            nowAnimals = nowSires #could be dams... actually no sex assigned
#            println("I used sires")
            elseif keepNoff[i] == 2
#            println("I used both")
            nowAnimals = concatCohorts(nowSires,nowDams)
                else println("THERE IS A MISTAKE")
        end
        outputPedigree(nowAnimals, "output$(thisBreed)_$(thisChr).txt")
        global    currentCohort = concatCohorts(currentCohort,nowAnimals) ## add to existing cohort
    end

    global    currentCohort = cohortSubset(currentCohort,collect((length(lastGenInd)+1):(length(lastGenInd)+(nBreedAMale+nBreedAFemale)))) #clean the cohort from the old ones

    global nLastInd = nBreedAMale+nBreedAFemale
    global nTotInd += nLastInd
    global generationCounterA = vcat(generationCounterA,fill(gen,nBreedAMale+nBreedAFemale))
end

pedigree = readdlm("output$(thisBreed)_$(thisChr).txt.ped",header=false)
writedlm("output$(thisBreed)_$(thisChr).txt.FINALped",Int.([generationCounterA pedigree]))
