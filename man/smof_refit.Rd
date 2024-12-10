\name{smof_refit}
\alias{smof_refit}

\title{Re-fitting an existing \code{smof} model for improved optimization}

\description{Given a model fitted by \code{smof}, this function helps 
to improve the achieved fitting level (represented by the value of the target
criterion) by launching new numerical optimizations 
of the underlying target function.
The search process is initiated with parameters randomly chosen
in the vicinity of those of the \code{object} provided.}

\usage{
smof_refit(object, searches = 10, sd = 1, opt.control=list(), verbose = 1)
}
 
\arguments{
  \item{object}{an object returned by a call to \code{smof};
        that call must include the argument \code{original=TRUE}. }
  
  \item{searches}{a positive integer representing the number of attempted searches.}
  
  \item{sd}{the standard deviation of the zero-mean normal variates used 
    for earch search to generate the initial set of parameters from which 
    to initiate a numerical optimization.}
  \item{opt.control}{a list passed to the optimization function (the same used 
    to create \code{object}) as its \code{control} argument. 
    It must not be used to turn the minimization problem into a maximization one.}      
  \item{verbose}{an integer regulating the amount of messages displayed;
    it can be 0 (no messages), 1 (default value) or larger than 1 for 
    a more verbose outcome.}
}

\details{
Section \sQuote{Computational aspects} of the \code{smof} documentation
provides some simple suggestions to overcome numerical difficulties arising  
in the fitting process.
For harder problems, \code{smof_refit} aims at improving the results obtained 
by an initial call to \code{smof} via a sequence of calls to \code{smof} itself, 
each time starting from the last best fit obtained up to that point.

In \sQuote{regular} situations, it is expected that the main usage of this 
function is to  improve the fitting of models produced with 
\code{scoring$type="spline"}, since this option typically generates a model 
with more parameters than an equivalent model with \code{scoring$type="distr"},
hence generally more problematic at the optimization step. 

Since the initial values of each optimization step are randomly generated
(in the vicininy of the best known parameter set), 
it is advisable to make a preliminary call to \code{set.seed}, to ensure
replicability of the results.
}
 
\value{an object of class \code{smof}}
 
\seealso{\code{\link[smof]{smof}},  \code{\link{set.seed}}}

\examples{
library(datasets)
data(esoph)
fit <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp + alcgp, family=binomial(), data=esoph)
smof2 <- smof(fit, esoph, c("agegp", "alcgp"), 
             scoring=list(type="spline", in.knots=c(2,1)), original=TRUE)
set.seed(1)             
smof2a <- smof_refit(smof2, searches=4, opt.control=list(maxit=50))
# smof2b <- smof_refit(smof2a)  # further improvement can be attempted
}

\keyword{models}
\keyword{classes}
\keyword{category}
\concept{ordered factor}
\concept{categorical data}
 
  