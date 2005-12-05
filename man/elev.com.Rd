\name{elev.com}
\alias{elev.com}
\title{ Test for equal elevation among several lines fitted with (standardised) major axes of common slope }
\description{
  Test if several major axis or standardised major axis lines share a common elevation. 
}
\usage{
elev.com(y, x, groups, data = NULL, method = 'SMA', alpha = 0.05,
  V = array(0, c(2, 2, length(unique(groups)))), group.names = sort(unique(groups)))
}
\arguments{
  \item{y}{ The Y-variable for all observations (as a vector) }
  \item{x}{ The X-variable for all observations (as a vector) }
  \item{groups}{ Coding variable identifying which group each observation belongs to (as a factor or vector) }
  \item{data}{ (optional) data frame containing the data }
  \item{method}{ The line fitting method: 
    \describe{
      \item{'SMA' or 1}{ standardised major axis (this is the default) }
      \item{'MA' or 2}{ major axis }
    }
  }
  \item{alpha}{ The desired confidence level for the 100(1-alpha)\% confidence interval for the common elevation. (Default value is 0.05, which returns a 95\% confidence interval.) }
  \item{V}{ The estimated variance matrices of measurement error, for each group. This is a 3-dimensional array with measurement error in Y in the first row and column, error in X in the second row and column,and groups running along the third dimension. Default is that there is no measurement error. }
  \item{group.names}{ (optional: rarely required). A vector containing the labels for `groups'. (Only actually useful for reducing computation time in simulation work). }
}
\value{
  \item{stat }{The Wald statistic testing for no shift along the common axis}
  \item{p }{The P-value of the test. This is calculated assuming that stat has a chi-square distribution with (g-1) df, if there are g groups}
  \item{a }{The estimated common elevation}
  \item{ci }{A 100(1-alpha)\% confidence interval for the true common elevation}
  \item{as }{Separate elevation estimates for each group}
}
\details{
  Calculates a Wald statistic to test for equal elevation of several MA's or SMA's with a common slope. This is done by testing for equal mean residual scores across groups.

  Note that this test is only valid if it is reasonable to assume that the axes for the different groups all have the same slope.

  The test assumes the following:
  \itemize{
    \item each group of observations was independently sampled
    \item the axes fitted to all groups have a common slope
    \item y and x are linearly related within each group
    \item residual scores independently follow a normal distribution with equal variance at all points along the line, within each group
  }
  Note that we do not need to assume equal variance across groups, unlike in tests comparing several linear regression lines.

  The assumptions can be visually checked by plotting residual scores against fitted axis scores, and by constructing a Q-Q plot of residuals against a normal distribution.  Residual scores can be obtained as y-bx, and for fitted axis scores y+bx (for SMA) or by+lx (for MA or `lamest'), where b represents the common slope estimate, and l the error variance ratio. If there is a distinct increasing or decreasing trend within any of the groups, this suggests that all groups do not share a common slope.

  A plot of residual scores against fitted axis scores can also be used as a visual test for common elevation. If residual scores generally differ across groups (with some groups generally having larger residual scores than others) then this is evidence that the groups do not share a common elevation.

  The common slope (\eqn{\hat{\beta}}{b}) is estimated from a maximum of 100 iterations, convergence is reached when the change in \eqn{\hat{\beta} < 10^{-6}}{b < 10^-6}.
}
\references{ Warton D. I., Wright I. J., Falster D. S. and Westoby M. (2006) A review of bivariate line-fitting methods for allometry.  \emph{Biological Reviews} (in press) }
\author{ Warton, D. \email{David.Warton@unsw.edu.au}, translated to R by Ormerod, J. 2005-12-08 }

\seealso{ \code{\link{line.cis}}, \code{\link{slope.com}}, \code{\link{shift.com}} }

\examples{
#load leaf longevity data
data(leaflife)

#Test for common SMA slope amongst species at low rainfall sites
#with different levels of soil nutrients
leaf.low.rain=leaflife[leaflife$rain=='low',]
slope.com(log10(longev), log10(lma), soilp, data=leaf.low.rain)

#Now test for common elevation of the groups fitted with an axis
#of common slope, at low rainfall sites:
elev.com(log10(longev), log10(lma), soilp, data = leaf.low.rain)

#Or test for common elevation amongst the MA's of common slope,
#for low rainfall sites, and construct 99% a confidence interval
#for the common elevation:
elev.com(log10(longev), log10(lma), soilp, method='MA',
    data = leaf.low.rain, alpha=0.01)  
}

\keyword{ htest }
