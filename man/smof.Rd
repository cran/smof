\name{smof}
\alias{smof}

\title{Scoring Methodology for Ordered Factors}

\description{Starting from an object representing a fitted model whose (non-)linear
  predictor includes some ordered factor(s) among the explanatory variables, 
  a new model is constructed where each named factor is replaced by a single
  numeric score, suitably chosen so that the new variable produces a data
  fit comparable with the standard methodology based on a set of polynomial
  contrasts. 
}
\usage{
smof(object, data, factors, scoring, fast.fit = FALSE, original = FALSE,
       f.tail=".score", opt.control=list(), verbose = 1)
}
\arguments{
  \item{object}{an object produced by a fitting function; see \sQuote{Details}
    below for specification of the admissible classes of objects. }
  \item{data}{the data frame used for producing \code{object}.}
  \item{factors}{a character vector with the names of the ordered factors of 
    \code{data} which must be converted to numeric scores.}
  \item{scoring}{a list which selects the type of scoring techniques and 
    related ingredients.
    Its key component is the character string \code{scoring$type} with possible
    values \code{"distr"} and \code{"spline"}; the other components of the
    list depend on this value, and are described in the  \sQuote{Details}.
    If \code{scoring} is missing, the default value 
    \code{list(type="distr", family="gh")} is used.}
  \item{fast.fit}{a logical value (default value: \code{FALSE}) indicating 
    whether a fast-fitting procedure must be used. 
    This option is available only under certain circumstances specified 
    in the \sQuote{Details} below.}
  \item{original}{a logical value (default value: \code{FALSE}) indicating
    whether the original \code{object} must be included in the returned object.} 
  \item{f.tail}{a character string representing the suffix to be added to
    each factor name when its numeric score variable is constructed (default
    value: \code{".score"}). It may be suitably adjusted in case the supplied 
    string leads to an invalid variable name.}
  \item{opt.control}{a list passed to \code{optim} as its \code{control}
    argument. It must not be used to turn the minimization problem into a
    maximization one.}   
  \item{verbose}{an integer regulating the amount of messages displayed;
    it can be 0 (no messages), 1 (default value) or larger than 1 for 
    a more verbose outcome.}   
}
\details{
In its original formulation,
function \code{smof} implements  the methodology proposed by 
Azzalini (2023), briefly summarized in the \sQuote{Background} section.
It is recommended to read at least that section in case the referenced paper 
is not examined. The published paper has open access.  
Later on, in version 1.2.0 of \pkg{smof}, a variant methodology has been 
included, based on the use of splines instead of quantile functions,
presented in Azzalini (2024).

Start from an \code{object} obtained as the outcome from some fitting 
procedure, whose linear predictor includes one or more ordered factor(s)
among the explanatory variables. 
For each ordered factor whose name is included in vector \code{factors}, 
a suitable vector of numeric scores is constructed.
The selection process examines the quantiles of the members of a specified 
parametric  class of distributions  and selects the member with optimizes 
(i.e. minimizes) a suitable target criterion.
To avoid trivialities, each factor in vector \code{factors} must have at 
least three levels.

There are two quite different options to build the numeric scores assigned
to an ordered factor. The selection of one of these options is made via
the component \code{type} of the list \code{scoring}, which 
can be either \code{"distr"} or \code{"spline"}. 
The other components of \code{scoring} depends on the chosen \code{type} 
and are described below.

If \code{scoring$type="distr"}, the numeric scores are obtained as quantiles
of a probability distribution belonging to a certain parametric family;
this route corresponds to the original construction of \code{smof},
following Azzalini (2023).
The admissible parametric families are all obtained by monotonic 
transformations of a standard normal variate.  
Specifically, the admissible families and corresponding strings to be
specified in \code{scoring$family} are as follows:
\tabular{lll}{
  Johnson's \eqn{S_U}            \tab \eqn{\quad}{  } \tab \code{"SU"} \cr
  Tukey's \emph{g}-and-\emph{h}  \tab \tab \code{"gh"}, \code{"g-and-h"} \cr
  Jones and Pewsey's \emph{sinh-arcsinh}
                            \tab \tab \code{"sinh-arcsinh"}, \code{"SAS"}  
}% end tabular
where either string name can be used when two of them are indicated. 
All these families involve two parameters for shape regulation; location and 
scale parameters are not considered, because irrelevant for our purposes.
Of the two shape parameters, the first one regulates asymmetry and can take
any value, while the second one regulates tail thickness and must be positive.
In each case, the adopted parameterization is the \sQuote{standard} one,
but explicit specifications are provided in the reference below.  
The same \code{family} is employed for all the components of \code{factors}.

If \code{scoring$type="spline"}, the numeric scores are obtained via a
monotonic spline function. Specifically, the scores are generated using
\code{splinefun} with \code{method="monoH.FC"}. Since currently this is 
the only admissible form of spline, it does not need to be specified.
What must be specified is \code{scoring$in.knots}, the number of internal knots
between the fixed extremal knots \code{1} and \code{K}, if \code{K} denotes
the number of levels of any given factor.  Hence \code{scoring$in.knots}
should be an integer vector with as many components as \code{factors};
if a shorter vector is supplied, its values will be recycled.
Since each internal knot involves the selection of two numeric values,
the total number of fitted parameters will equal the sum of twice the
\code{in.knots} values, summed over the components of \code{factors}.

Estimation of the transformation parameters is performed by minimization,
performed by the optimizer \code{optim}, of a suitable target criterion,
which depends on \code{class(object)}.
The admissible classes for \code{object} are currently as follows,
listed along the corresponding target criteria:
\tabular{lll}{
  \kbd{class} \tab \kbd{fitting function (package)} \tab \kbd{target criterion}\cr
  \eqn{\rule[0.8ex]{5em}{0.02ex}}{------------}\tab 
       \eqn{\rule[0.8ex]{12em}{0.02ex}}{----------------------------}\tab 
       \eqn{\rule[0.8ex]{10em}{0.02ex}}{-------------}\cr
   lm           \tab lm (stats)         \tab sum of squared residuals \cr
   mlm          \tab lm (stats)         \tab [see below] \cr
   glm          \tab glm (stats)        \tab deviance \cr
   survreg      \tab survreg (survival) \tab \eqn{-}loglikelihood \cr 
   coxph        \tab coxph (survival)   \tab \eqn{-}loglikelihood \cr
   coxph.penal  \tab coxph (survival)   \tab \eqn{-}loglikelihood \cr
   gnm          \tab gnm (gnm)          \tab deviance
} % end tabular
For an object of class \code{mlm}, the target function is formed by summing
terms where the contribution  from the \eqn{j}-th response variable is
\eqn{(1-R^2_j)}, where \eqn{R^2_j} is the \code{r-squared} statistic for 
that component of the fitted model. 
Note that, in the case of a single response variable,
its \eqn{(1-R^2)} value is equivalent, up to an algebraic transformation, 
to the sum of squared residuals used for \code{lm} objects;
hence the chosen target criterion for \code{mlm} models is a direct extension
of the one for \code{lm}'s.
The above list of classes may be expanded in the future, depending on feedback.

In certain cases, models involving  a large number of parameters for 
regulation of the  scores transformation(s) can lead to difficulties 
at the optimization step, which is performed by \code{optim}.
This situation may be flagged by the message
\preformatted{Non-zero error code returned by optim:...}
In these circumstances, one direction for improvement is to make use of
the argument \code{opt.control} for better tuning the working of \code{optim}. 
In more problematic cases, the auxiliary function \code{smof_refit} is
provided to improve upon an initial fitting.

The rest of this section is slightly of more technical nature, 
and it may be not of interest to the casual user, especially if
the option \code{fast.fit=TRUE} is not selected.
Operationally, estimation of the \code{scoring$family} parameters is performed 
via optimization of the pertaining target criterion,
as indicated by the table above.
For each candidate set of parameters, each factor included in 
\code{factors} is replaced by values determined by the 
quantiles of  \code{scoring$family}  and the current parameters.
The name of the new constructed variable is formed by adding \code{.score} 
to the original name.
For instance, an ordered factor called \code{ordfac} is replaced 
by the numeric variable \code{ordfac.score} both 
in the linear predictor of \code{object} and in the \code{data} frame.
A call to \code{update} using the modified linear predictor and data delivers
a new fitting, with attached a value of the  target criterion.
An iterative optimization process of the target criterion leads to 
the estimated parameters of \code{scoring$family}  with a corresponding 
fitted model.

If the component \code{scoring$param} is not \code{NULL}, it will be used for
parameter initialization in the parameter search process. In this case,
if $\code{scoring$type="distr"}, the even ordered components of \code{scoring$param} 
must be positive, since they represent tail-weight parameters.
If $\code{scoring$type="distr"}, they must be a vector of length \code{2*in.knots}.

There are in fact two variants of the procedure. 
What has been just described refers to the more \sQuote{general} variant form.
However, in the prominent cases of an \code{object} of class \code{lm} 
or \code{glm}, the procedure can be speeded-up by setting \code{fast.fit=TRUE}, 
provided the fitted model is of a basic form, that is, a model specification
via  a formula, and a \code{family} in the \code{glm} case, without  
non-basic arguments such as \code{offset}, \code{subset} and alike.
If these non-basic arguments are included in the \code{object} call, 
they are ignored for estimation of the \code{scoring$family} parameter.
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
A list with the following components:
  \item{call}{the calling statement}
  \item{new.object}{an updated version of the original \code{object}, with 
    the components of \code{factors} in the model replaced by new variables; 
    this \code{object} is itself a list, whose structure depends on its class.}
  \item{new.data}{a new data frame where the ordered factors are replaced 
     by numeric variables representing scores.}
  \item{scoring}{a list similar to the input argument with the addition of
     the estimates of the parameters.}
  \item{factor.scores}{a list of numeric vectors with the scores assigned 
    to the levels of each factor.}   
  \item{original.factors}{a list with the names and the levels of the original 
    \code{factors}.} 
  \item{target.criterion}{the final value of the target criterion used for fitting.}
  \item{opt}{the list returned by \code{optim}.}
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

There are two variants of the methodology, selected with the value of
the component \code{type} of the list \code{scoring}.
Here we summarize the working of the original formulation, selected by
setting \code{type="distr"}. 
The actual construction of numeric scores proceeds by selecting \emph{K} 
quantiles of a distribution belonging to some parametric family.
The adoption of a sufficiently flexible parametric family helps to find
a scoring system best suited for the data under consideration,
hence improving upon the basic sequence \emph{1, ..., K}.
A concomitant product of this scheme is the identification of numeric values
which indicate how the \emph{K} levels are \dQuote{really} spaced.
Combining these two features, the key feature of the proposal is 
interpretability of the construction.

The proposed method represents an alternative to the use of polynomial 
contrasts, which is the default action taken by \R for ordered factors;
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

The alternative variant of the methodology is selected by setting  
\code{scoring$type="spline"}. The underlying principle is similar to the
one just described, but it makes use of splines instead of quantile functions.
Its operational working is described in the \sQuote{Details}.

For a more comprehensive exposition and discussion, see the references below.
}

\note{
For subsequent computations on the object returned by \code{smof}, difficulties 
may arise if the call to the fitting function does not set \code{model=TRUE}.
This is not a problem with \code{lm} and \code{glm}, if their default setting
\code{model=TRUE} has not been modified.
The  default setting of \code{coxph} is instead \code{model=FALSE}.
This implies, for instance, that issuing the \code{survival} command   
\code{survfit(smof4$new.object)}, right after running 
the code of Example 4 below, would cause an error. 
The main route to avoid this issue is to set \code{model=TRUE} in the call to
the fitting function, that is, \code{coxph} or whaterever function is used.
Alternatively, if one does not want to refit an already existing object,
there exist various ways to overcome this snag; the simplest one is to write  
\tabular{rl}{
  \tab \code{new.data <- smof4$new.data} \cr
  \tab \code{s <- survfit(smof4$new.object)}
} 
This indication is temporary and it may be superseded by a different design 
in future versions of the package.
}

\references{
Azzalini, A. (2023). On the use of ordered factors as explanatory variables.
  \emph{Stat} \bold{12}, e624. \doi{10.1002/sta4.624}
  
Azzalini, A. (2024). On the use of splines for representing ordered factors.
  \emph{arXiv:2406.15933}, \doi{10.48550/arXiv.2406.15933}
}

\author{Adelchi Azzalini}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{smof_refit}}, \code{\link[stats]{optim}},
  \code{\link[stats]{contr.poly}}, \code{\link[stats]{update}}, 
  \code{\link[stats]{splinefun}}, \code{\link[stats]{plogis}},
  \code{\link[stats]{lm}}, \code{\link[stats]{lm.fit}}, 
  \code{\link[stats]{glm}}, \code{\link[stats]{glm.fit}} 
}

\examples{
# Example 1, reconstructs Tables 1 and 2 (second part) of the reference
message("--- Example 1: esoph data ---")
library(datasets)
data(esoph)
contrasts(esoph$agegp, 2) <- contr.poly(6)  
contrasts(esoph$tobgp, 1) <- contr.poly(4)  
fit1 <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp + alcgp, family=binomial(), data=esoph) 
message("original fit:") 
print(summary(fit1))      
smof1 <- smof(fit1, esoph, "alcgp")
# to select the Johnson's S_U family of distributions, write instead:
# smof1 <- smof(fit1, esoph, "alcgp", scoring=list(type="distr", family="SU"))
print(smof1)
print(summary(smof1))
plot(smof1, type="b", pch=20, col=4)
#
# Example 2 , reconstructs Tables 3 and 4 (first part) of the reference
if(require(ggplot2, quietly=TRUE)) {
message("--- Example 2: diamonds data ---")
data(diamonds, package="ggplot2")  
dmd <- data.frame(diamonds[seq(1, 53940, by=100),]) # use a subset of the data
dmd <- dmd[-c(518, 519, 523),] # remove three outliers
contrasts(dmd$clarity, 3) <- contr.poly(8)  
contrasts(dmd$color, 4) <- contr.poly(7)
contrasts(dmd$cut, 1) <- contr.poly(5) 
fit2 <- lm(sqrt(price) ~ carat + clarity + color + cut, data=dmd)
smof2 <- smof(fit2, dmd,  c("color", "clarity"))
message("smof fit:") 
print(smof2)
print(summary(smof2))
plot(smof2, which="clarity")
} # end diamonds example
#
# Example 3
if(require(survival, quietly=TRUE)) {
message("--- Example 3: lung data ---")
lung0 <- lung
lung0$ph.karno <- ordered(lung0$ph.karno)
contrasts(lung0$ph.karno, 3) <- contr.poly(6)
fit3 <- survreg(Surv(time, status) ~ ph.karno, data=lung0)
smof3 <- smof(fit3, lung0, "ph.karno")
print(summary(smof3))
plot(smof3)  # Karnofsky scores do not seem to be linearly spaced
# 
message("--- Example 4: PBC data ---")
data(pbc, package="survival")
pbc$stage <- ordered(pbc$stage)
fit4 <- coxph(Surv(time) ~ strata(status) + stage, data=pbc)
smof4 <- smof(fit4, data=pbc, factors="stage")
print(summary(smof4))
plot(smof4)
} # end of survival examples
}
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
\keyword{classes}
\keyword{category}
\concept{ordered factor}
\concept{categorical data}
 

  
