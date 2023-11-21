\name{smof-package}
\alias{smof-package}
% \alias{smof}
\docType{package}
\title{
% \packageTitle{smof}
Scoring Methodology for Ordered Factors
}
\description{
% \packageDescription{smof}
Starting from a given object representing a fitted model (within a
  certain set of model classes) whose linear predictor includes some ordered 
  factor(s) among the explanatory variables, a new model is constructed and 
  fitted where each named factor is  replaced by a single numeric score, 
  suitably chosen so that the new variable produces a fit comparable 
  with the standard methodology based on a set of polynomial contrasts. 
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{smof}
\packageIndices{smof}
% An overview of how to use the package, including the most important functions  
}

\author{
% \packageAuthor{smof}
% Maintainer: \packageMaintainer{smof}
Author: Adelchi Azzalini [aut, cre] (<https://orcid.org/0000-0002-7583-1269>)                   
Maintainer: Adelchi Azzalini <adelchi.azzalini@unipd.it>
}
\references{
Azzalini, A. (2023). On the use of ordered factors as explanatory variables.
  \emph{Stat} \bold{12}, e624. \doi{10.1002/sta4.624}
}

\keyword{ package }

\examples{
library(datasets)
data(esoph)
contrasts(esoph$agegp, 2) <- contr.poly(6) # optional
contrasts(esoph$tobgp, 1) <- contr.poly(4) # optional
obj1 <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp + alcgp, family=binomial(), data=esoph)        
out0 <- smof(obj1, esoph, "alcgp", distr.type="gh")
print(summary(out0$object)) 
}
