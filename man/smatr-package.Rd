\name{smatr}

\alias{smatr-package}
\alias{smatr}

\title{(Standardised) Major Axis Estimation and Testing Routines}

\description{This package provides methods of fitting bivariate lines in allometry using the major axis (MA) or standardised major axis (SMA), and for making inferences about such lines. The available methods of inference include confidence intervals and one-sample tests for slope and elevation, testing for a common slope or elevation amongst several allometric lines, constructing a confidence interval for a common slope or elevation, and testing for no shift along a common axis, amongst several samples.}

\details{
The key functions available in this package are the following.

\bold{For one sample:}
    \describe{
      \item{line.cis}{ fit allometric lines and construct confidence intervals for the true slope and elevation }
      \item{slope.test}{ test if the slope of an allometric line equals a specific value }
      \item{elev.test}{ test if the elevation of an allometric line equals a specific value }
    }
  
\bold{For several samples:}
    \describe{
      \item{slope.com}{ test if several major axis or standardised major axis lines share a common slope, and construct a confidence interval for the true common slope }
      \item{elev.com}{ test if several major axis or standardised major axis lines share a common elevation, and construct a confidence interval for the true common elevation }
      \item{shift.com}{ test if several groups of observations have no shift in location along major axis or standardised major axis lines with a common slope }
    }
 
Options are available for forcing lines through the origin (where appropriate), and all procedures have the option of correcting for measurement error (although only in an approximate fashion, valid in large samples).

Additional features of this package are listed below.

\bold{Measurement error estimation:}
    \describe{
      \item{meas.est}{ Estimates the average variance matrix of measurement error for a set of subjects with repeated measures }
    }

\bold{Example datasets:}
    \describe{
      \item{leaflife}{ leaf longevity and leaf mass per area for plant species from different sites. Used to demonstrate the functionality of one-sample and multi-sample procedures }
      \item{leafmeas}{ leaf mass per area and photosynthetic rate for plant species from different sites. Used to demonstrate the meas.est function }
    }

For more details, see the documentation for any of the individual functions listed above.

}
\references{ Warton D. I. and Weber N. C. (2002) Common slope tests for bivariate structural relationships.  \emph{Biometrical Journal} \bold{44}, 161--174.

Warton D. I., Wright I. J., Falster D. S. and Westoby M. (2006) A review of bivariate line-fitting methods for allometry.  \emph{Biological Reviews} (in press)
}

\author{ Warton, D. \email{David.Warton@unsw.edu.au}, translated to R by Ormerod, J. 2005-12-08 }

\seealso{\code{\link{line.cis}}, \code{\link{slope.test}}, \code{\link{elev.test}}, \code{\link{slope.com}}, \code{\link{elev.com}}, \code{\link{shift.com}}, \code{\link{meas.est}}, \code{\link{leaflife}}, \code{\link{leafmeas}} }

\examples{
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

    #now construct a scatterplot, with different colours for different sites
    plot(leaflife$lma,leaflife$longev,type="n",xlab="leaf mass per area [log scale]",ylab="leaf longevity [log scale]",log="xy")
    colours <- c("blue", "red", "green", "yellow")
    points(leaflife$lma,leaflife$longev,col=colours[as.numeric(leaflife$site)])
    legend(55,5,as.character(unique(leaflife$site)),col=colours,pch=rep(1,4))

    #test for a common SMA slope amongst species from sites with different rainfall/nutrients:
    fit.slopes <- slope.com(log10(longev), log10(lma), site, data = leaflife)

    #Test for common SMA slope amongst species at low rainfall sites with different levels of soil nutrients
    slope.com(log10(longev), log10(lma), soilp, data=leaf.low.rain)
    
    #Now test for common elevation of the groups fitted with an axis of common slope, at low rainfall sites:
    elev.com(log10(longev), log10(lma), soilp, data = leaf.low.rain)

    #Now test for no shift along the axes of common slope, for sites with different soil nutrient levels but low rainfall:
    shift.com(log10(longev), log10(lma), soilp, data=leaf.low.rain)

    #Test for common major axis slope, and construct 90% confidence intervals for common slope and each separate slope:
    slope.com(log10(longev), log10(lma), site, data=leaflife, method="MA", alpha=0.1)

    #Test for common elevation amongst the MA's of common slope, for low rainfall sites, and construct 99% confidence intervals for all elevation estimates:
    elev.com(log10(longev), log10(lma), soilp, method="MA", data = leaf.low.rain, alpha=0.01)
}
    
\docType{ package }
\keyword{ documentation }
