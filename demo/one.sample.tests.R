    #load the leaflife dataset:
    data(leaflife)

    #consider only the low rainfall sites:
    leaf.low.rain=leaflife[leaflife$rain=='low',]

    #construct a plot
    plot(log10(leaf.low.rain$lma), log10(leaf.low.rain$longev), xlab='leaf mass per area [log scale]', ylab='leaf longevity [log scale]')

    #test if the SMA slope amongst species at low rainfall sites is 1, for log (base 10) transformed data:
    slope.test(log10(longev), log10(lma), data=leaf.low.rain)
    
    #test if the MA slope is 2/3
    slope.test(log10(longev), log10(lma), data=leaf.low.rain, test.value = 2/3, method = 'MA')

    #produce CI's for MA slope and elevation:
    line.cis(log10(longev),log10(lma),data=leaf.low.rain, method=2)
