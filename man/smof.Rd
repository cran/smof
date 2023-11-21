\name{smof}
\alias{smof}

\title{Scoring Methodology for Ordered Factors}

\description{Starting from an object representing a fitted model whose linear
  predictor includes some ordered factor(s) among the explanatory variables, a
  new model is constructed where each named factor is replaced by a single
  numeric score, suitably chosen so that the new variable produces a data
  fit comparable with the standard methodology based on a set of polynomial
  contrasts. 
}
\usage{
smof(object, data, factors, distr.type = "gh", fast.fit = FALSE, trace = FALSE)
}
\arguments{
  \item{object}{an object produced by a fitting function; see \sQuote{Details}
    below for specification of the admissible classes of objects. }
  \item{data}{the data frame used for producing \code{object}.}
  \item{factors}{a character vector with the names of the ordered factors of 
    \code{data} which must be converted to numeric scores.}
  \item{distr.type}{a character string with the name of the parametric family 
    of distributions used to construct the numeric scores.  
    See \sQuote{Details} for the set of admissible choices; 
    default value: \code{"gh"}.}
  \item{fast.fit}{a logical value (default value: \code{FALSE}) indicating 
    whether a fast-fitting procedure must be used. 
    This option is available only under certain circumstances specified 
    in the \sQuote{Details} below.}
  \item{trace}{a logical value (default value: \code{FALSE}) indicating 
    whether details of the iterative procedure for parameter estimation must 
    be printed at each iteration.}
}
\details{
Function \code{smof} implements the methodology proposed by Azzalini (2023),
briefly summarized in the \sQuote{Background} section. 
It is advisable to read at least this section in case the referenced paper 
is not examined. The published paper is open access.  

Start from an \code{object} obtained as the outcome from some fitting 
procedure, whose linear predictor includes one or more ordered factor(s)
among the explanatory variables. 
For each ordered factor whose name is included in the vector \code{factors}, 
a suitable vector of numeric scores is constructed.
The selection process examines the quantiles of the members of a specified 
parametric  class of distributions  and selects the member with optimizes 
(i.e. minimizes) a pertaining target criterion.

The admissible parametric families whose quantiles are used to construct 
the scores of the factors are all obtained by monotonic transformations of
a standard normal variate.  
Specifically, the admissible families and corresponding strings to be
specified in \code{distr.type} are as follows:
\tabular{lll}{
  Johnson's \eqn{S_U}            \tab \tab\code{"SU"} \cr
  Tukey's \emph{g}-and-\emph{h}  \tab \tab \code{"gh"},  \code{"g-and-h"} \cr
  Jones and Pewsey's \emph{sinh-arcsinh}
                            \tab \tab \code{"sinh-arcsinh"}, \code{"SAS"} \cr
}% end tabular
where either string name can be used when two of them are indicated. 
Each of these families involve two parameters for shape regulation.
The same \code{distr.type} is employed for all components of \code{factors}; 
only the chosen parameters differ.

The admissible classes for \code{object} are currently as follows,
listed along the corresponding target criteria:
\tabular{lll}{
  \kbd{class} \tab \kbd{fitting function (package)} \tab \kbd{target criterion}\cr
  \kbd{------------}\tab \kbd{----------------------------}\tab \kbd{-------------}\cr
   lm           \tab lm (stats)         \tab sum of squared residuals \cr
   mlm          \tab lm (stats)         \tab [see below] \cr
   glm          \tab glm (stats)        \tab deviance \cr
   survreg      \tab survreg (survival) \tab \eqn{-}loglikelihood \cr 
   coxph        \tab coxph (survival)   \tab \eqn{-}loglikelihood \cr
   coxph.penal  \tab coxph (survival)   \tab \eqn{-}loglikelihood 
} % end tabular
For an object on class \code{mlm}, the target function is formed by summing
terms where the contribution  from the \eqn{j}-th response variable is
\eqn{(1-R^2_j)}, where \eqn{R^2_j} is the \code{r-squared} statistic for 
that component of the fitted model. 
Note that, in the case of a single response variable,
its \eqn{(1-R^2)} value is equivalent, up to an algebraic transformation, 
to the sum of squared residuals used for \code{lm} objects;
hence the chosen target criterion for \code{mlm} models is a direct extension
of the one for \code{lm}'s.
The above list of classes may be expanded in the future, depending on feedback.

The rest of this section is slightly of more technical nature, 
and it may be not of interest to the casual user, especially if
the option \code{fast.fit=TRUE} is not activated.
Operationally, estimation  of the  \code{distr.type}  parameters 
is performed via optimization of the pertaining target criterion,
as indicated by the table above.
For each candidate set of parameters, each factor included in 
\code{factors} is replaced by values determined by the 
quantiles of  \code{distr.type}  and the current parameters.
The name of the new constructed variable is formed by adding \code{.score} 
to the original name.
For instance, an ordered factor called \code{ordfac} is replaced 
by the numeric variable \code{ordfac.score} both 
in the linear predictor of \code{object} and in the \code{data} frame.
A call to \code{update} using the modified linear predictor and data delivers
a new fitting, with attached a value of the  target criterion.
An interative optimization process the target criterion leads to 
the estimated parameters of \code{distr.type}  with a corresponding 
fitted model.

There are in fact two variants of the procedure. 
What has been just described refers to the more \sQuote{general} variant form.
However, in the prominent cases of an \code{object} of class \code{lm} 
or \code{glm}, the procedure can be speeded-up by setting \code{fast.fit=TRUE}, 
provided the fitted model is of a basic form, that is, a model specification
via  a formula, and a \code{family} in the \code{glm} case, without  
non-basic arguments such as \code{offset}, \code{subset} and alike.
If these non-basic arguments are included in the \code{object} call, 
they are ignored for estimation of the \code{distr.type} parameter.
However, they are included for producing the final \code{object} 
returned by the function. 
With this option, the sequence of calls to \code{lm} and \code{glm}
involved by the iterative search procedure is replaced by faster calls
to  \code{lm.fit} and \code{glm.fit}.
Correspondingly, the internal target function (\code{target.fit}) is slightly
different from the one used on the more general case (\code{target.gen}). 
Since the selection of the parameters involves an iterative process with
dimensionality equal to twice the length of \code{factors} and each iteration
involves a new data fitting process, the saving in execution time can be 
substantial in some cases.

} % end details

\value{
A list with the following components
  \item{object}{an updated version of the original \code{object}; 
        this is itself a list, whose structure depends on its class.}
  \item{data}{the new data frame with the  ordered factors 
        replaced by numeric score variables.}
  \item{distr}{a list with two components: (1) \code{type}, the name of
     of the distribution type; (2) \code{param}, a matrix of fitted parameters  
     with as many rows as \code{length(factors)}.}
  \item{factor.scores}{a list of numeric vectors with the scores assigned 
      to the levels of each factor.}    
}

\section{Background}{
The methodology proposed in the reference below deals with the presence of 
ordered factors used as explanatory variables, hence included in the linear 
predictor of some model under consideration.
For any given ordered factor with \emph{K} levels, say, a set of \emph{K} 
numeric scores is introduced, with a certain value assigned to each factor level.
In the end, the original factor is effectively replaced by a numeric variable.
This scheme represents a refinement of the elementary scoring system based
on the basic sequence \emph{1, ..., K}, which constitutes a simple 
time-honoured option to deal with ordered factors, 
but it is not always appropriate.

The actual construction of numeric scores proceeds by selecting \emph{K} 
quantiles of a distribution belonging to some parametric family.
The adoption of a sufficiently flexible parametric family helps to find
a scoring system best suited for the data under consideration,
hence improving upon the basic sequence \emph{1, ..., K}.
A concomitant product of this scheme is the identification of numeric values
which indicate how the \emph{K} levels are \dQuote{really} spaced.
Combining these two features, the key feature of the proposal is 
interpretability of the construction.

The proposed method represents an alternative to the use of polynomial contrasts, 
which is the default action taken by \R for ordered factors;
see the documentation of \code{contr.poly}.

In the proposed logic, the constructed scores are intended to be used, and
interpreted, without further manipulation. Hence, for instance, building a 
polynomial form using one such variable would diverge somewhat from the 
proposed logic, although still conceivable.
With a single numeric variable to represent a given factor, 
one cannot expect to achieve the same numerical fit to the data 
as obtained the polynomial contrasts built for the original factor, 
when these constrasts involve high degrees polynomials, and 
correspondingly several parameters.
However, a range of numerical explorations has indicated that in many cases
the resulting fit is equal or similar to the one achieved via polynomial 
constrasts, with non-negligible simplification in the model specification, 
and easier interpretation, 

In a nutshell, the aim of the approach is to achieve a satisfactory data fit
while improving an model parsimony, with simple interpretability of the 
score system.

For a more comprehensive exposition and discussion, see the reference below.
}


\references{
Azzalini, A. (2023). On the use of ordered factors as explanatory variables.
  \emph{Stat} \bold{12}, e624. \doi{10.1002/sta4.624}
}

\author{Adelchi Azzalini}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[stats]{contr.poly}}, \code{\link[stats]{update}}, 
  \code{\link[stats]{lm.fit}}, \code{\link[stats]{glm.fit}} }

\examples{
# Example 1, reconstructs Table 2 (fist part) of the reference
message("--- Example 1: esoph data ---")
library(datasets)
data(esoph)
contrasts(esoph$agegp, 2) <- contr.poly(6) # optional
contrasts(esoph$tobgp, 1) <- contr.poly(4) # optional
obj1 <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp + alcgp, family=binomial(), data=esoph) 
message("original fit:") 
print(summary(obj1))      
out1 <- smof(obj1, esoph, "alcgp", distr.type="SU")
message("smof fit:") 
print(summary(out1$object)) 
print(out1$distr)
print(out1$factors.scores)
#
# Example 2 , reconstructs Table 4 (fist part) of the reference
if(require(ggplot2, quietly=TRUE)) {
message("--- Example 2: diamonds data ---")
data(diamonds)  
dmd <- data.frame(diamonds[seq(1, 53940, by=100),]) # use a subset of the data
dmd <- dmd[-c(518, 519, 523),] # remove three outliers
contrasts(dmd$cut, 1) <- contr.poly(5) 
obj2 <- lm(sqrt(price) ~ carat + clarity + color + cut, data=dmd)
message("original fit:") 
print(summary(obj2))      
out2 <- smof(obj2, dmd,  c("color", "clarity"), distr.type="gh")
message("smof fit:") 
print(summary(out2$object))
print(out2$distr)
print(out2$factors.scores)
} # end diamonds example
#
# Example 3
if(require(survival, quietly=TRUE)) {
message("--- Example 3: lung data ---")
lung0 <- lung
lung0$ph.karno <- ordered(lung0$ph.karno)
contrasts(lung0$ph.karno, 3) <- contr.poly(6)
lfit3 <- survreg(Surv(time, status) ~ ph.karno, data=lung0)
out3 <- smof(lfit3, lung0, "ph.karno")
print(summary(out3$object))
print(out3$factors.scores) # not equally spaced scores
} # end lung example
}
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
\keyword{classes}
\keyword{category}
\concept{ordered factor}
\concept{categorical data}
 
