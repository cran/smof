\name{smof-methods}
\alias{print.smof}
\alias{plot.smof}
\alias{predict.smof}
\alias{summary.smof}
\alias{print.summary.smof}
\title{Methods for \code{smof} objects}

\description{The list of methods that apply to \code{smof} objects}

\usage{
\method{print}{smof}(x, ...)
\method{plot}{smof}(x, which, ...)
\method{summary}{smof}(object,  ...)
\method{print}{summary.smof}(x, ...)
\method{predict}{smof}(object, newdata, ...)
}

\arguments{
  \item{object}{an object returned by \code{smof}.}
  \item{x}{an object returned by \code{smof}, except for \code{print.summary.smof} 
    where \code{x} is the outcome of \code{summary.smof}.}
  \item{which}{either a vector or a two-elements list; see \sQuote{Details}
    for full specification.}
  \item{newdata}{a data frame which includes ordered factors with the same names
    and levels as those in the data frame used to produce \code{object};
    see \sQuote{Details} for additional information.}
 \item{...}{arguments passsed through to other methods.}   
}

\details{
There are two main parts in the outcome of \code{summary.smof}. One is the
outcome of the selection of the factor(s) transformation(s); the other part,
denoted \kbd{Final fitting call}, represents the newly fitted model  
with the original ordered factors  replaced  by numeric scores. 
The corresponding \code{data.frame}, denoted \code{new.data}, is one 
of the components of the object returned by \code{smof}.

If \code{which} is a character vector, its components are interpreted as
names of the \code{factors} in the calling statement of the object to be 
plotted, producing a set of
graphs where the numeric scores of each named factor are plotted versus
the equally spaced ticks associated to the original level names.
The same effect is obtained when \code{which} is a numeric vector of integers,
which then select the corresponding components of the \code{factors} sequence.

If \code{which} is a list, its first element is assumed to be a vector having
the meaning just described. After the pertaining set of graphs has been 
completed,  the second element of the list is passed to the plotting method 
for the object produced by the model fitting procedure.
Currently this option operates only for objects which inherits from class
\code{lm}; specifically, it works for objects originated by a call to 
\code{lm} or to \code{glm}.

With \code{predict.smof}, the outcome of a \code{smof} fit is applied 
to new data frame which includes ordered factors analogous to those used to 
compute the \code{smof} transformation. 
Only the factors of \code{newdata} with the same name as those processed by 
\code{smof} to produce \code{object} are examined. The levels of these
factors must coincide with or be a subset of those of the original data frame.

Note the difference between the role played by \code{predict.smof} and  
the one of \code{predict} methods for most other classes. Usually the 
returned values pertain to the response variable, or to some related entity, 
while here the outcome refers to explanatory variables.
}
 
\value{For \code{summary.smof}, a list of class \code{summary.smof}. 
  For \code{predict.smof}, a data frame.
  For \code{plot.smof}, \code{NULL} with graphical side effects. 
}

\author{Adelchi Azzalini}

% \note{  ~~further notes~~ }

 
\seealso{\code{\link{smof}}, \code{\link[stats]{lm}}, \code{\link[stats]{glm}} }

\examples{
library(datasets)
data(esoph)
contrasts(esoph$agegp, 2) <- contr.poly(6)  
contrasts(esoph$tobgp, 1) <- contr.poly(4)  
fit1 <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp + alcgp, family=binomial(), data=esoph)      
smof1 <- smof(fit1, esoph, "alcgp")
print(smof1)
print(summary(smof1))
plot(smof1, type="b", pch=19, col="blue")
plot(smof1, which=list(1, 1:4))
predict(smof1, newdata=esoph[seq(1, 88, by=8), ]) 
}

\keyword{models}
\keyword{classes}
\keyword{category}
\concept{ordered factor}
\concept{categorical data}