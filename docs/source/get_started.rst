Get Started
===========

An example to run ``XSim`` is shown below.  

.. code-block:: julia

	using XSim

    #set genome information
    chrLength, numChr, numLoci, mutRate = 1.0, 1, 100, 0.0
    locusInt  = chrLength/numLoci
    mapPos    = [0:locusInt:(chrLength-0.0001)];
    geneFreq  = fill(0.5,numLoci);

    XSim.init(numChr,numLoci,chrLength,geneFreq,mapPos,mutRate)
    pop1 = startPop()

    #generate populations
    ngen,popSize    = 10,10

    pop1.popSample(ngen,popSize)
    pop2 = pop1.popNew(10);
    pop3 = popCross(5,pop1,pop2);

    #generate genotypes
    M = pop3.getGenotypes()
