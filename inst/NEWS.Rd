\name{NEWS}
\title{News for Package \sQuote{smof}}
\encoding{UTF-8}
%---   
\section{Changes in smof version 1.2.1 (2024-10-30)}{
\itemize{
  \item Objects of class \code{gnm} are now handled.
  \item Methods \code{coef} and \code{anova} are introduced for class \code{smof}.
  \item The \code{smof} argument \code{trace} is replaced by the more flexible
        argument \code{verbose}.
  \item Fixed a bug in the parsing of the object formula which caused problems
        when some factor has a name which is portion of another factor name.
  \item Fixed a bug in \code{smof_refit} affecting the retrival of the 
        parameters estimated in the initial call to \code{smof}.       
  }
}
%---   
\section{Changes in smof version 1.2.0 (2024-09-17)}{
\itemize{
  \item A variant method for assigning numeric scores is introduced, based on 
        monotonic splines.       
  \item An implication of the new variant method is that one argument of the 
        function \code{smof} had to redesigned.
  \item The new function \code{smof_refit} is introduced.    
  \item Additional arguments of \code{smof} have been introduced.  
  \item Fixed a bug in \code{predict.smof} which affected the case with multiple factors.
  }
}  
%---   
\section{Changes in smof version 1.1.0 (2024-03-04)}{
\itemize{
  \item The object returned by \code{smof} now includes additional components
         and is assigned to class \code{smof}.
  \item Various methods are introduced for class \code{smof}: 
        \code{print, summary, plot, predict}.
  \item Fixed a bug which prevented the use of option \code{fast.fit=TRUE}
        when \code{family=binomial(..)} and the response was of binary type.
  }
}
%---   
\section{Initial smof version 1.0.0 (2023-11-21)}{
  \itemize{
  \item First public version of the package and upload to CRAN, with
        only one exported function, namely \code{smof}.
  \item Scoring of factors are produced from distributions among the families
        \code{"gh", "SU", "sinh-arcsinh"} or their equivalente names.
  \item Admissible classes of fitted models are \code{"lm", "mlm", "glm",
         "survreg", "coxph.penal", "coxph"}.  
  }
}