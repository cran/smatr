    #load the leaflife dataset:
    data(leaflife)

    #construct a scatterplot, with different colours for different sites
    plot(leaflife$lma,leaflife$longev,type="n",xlab="leaf mass per area [log scale]",ylab="leaf longevity [log scale]",log="xy")
    colours <- c("blue", "red", "green", "yellow")
    points(leaflife$lma,leaflife$longev,col=colours[as.numeric(leaflife$site)])
    legend(55,5,as.character(unique(leaflife$site)),col=colours,pch=rep(1,4))

    #test for a common SMA slope amongst species from sites with different rainfall/nutrients:
    fit.slopes <- slope.com(log10(longev), log10(lma), site, data = leaflife)

    #Test for common SMA slope amongst species at low rainfall sites with different levels of soil nutrients
    leaf.low.rain=leaflife[leaflife$rain=="low",]
    slope.com(log10(longev), log10(lma), soilp, data=leaf.low.rain)
    
    #Now test for common elevation of the groups fitted with an axis of common slope, at low rainfall sites:
    elev.com(log10(longev), log10(lma), soilp, data = leaf.low.rain)

    #Now test for no shift along the axes of common slope, for sites with different soil nutrient levels but low rainfall:
    shift.com(log10(longev), log10(lma), soilp, data=leaf.low.rain)

    #Test for common major axis slope, and construct 90% confidence intervals for common slope and each separate slope:
    slope.com(log10(longev), log10(lma), site, data=leaflife, method="MA", alpha=0.1)

    #Test for common elevation amongst the MA's of common slope, for low rainfall sites, and construct 99% confidence intervals for all elevation estimates:
    elev.com(log10(longev), log10(lma), soilp, method="MA", data = leaf.low.rain, alpha=0.01)

