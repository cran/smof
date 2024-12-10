\name{NEWS}
\title{News for Package \sQuote{smof}}
\encoding{UTF-8}
%---   
\section{Changes in smof version 1.2.2 (2024-12-10)}{
\itemize{
  \item The new \code{smof} argument \code{opt.method} is introduced 
        from which a wider choice of numerical optimizers is made available,
        that is, three optimizers from the \pkg{base} package are accessibile.
        In addition, if the package \pkg{nloptr} is installed on the local 
        system, various additional algorithms are made available.
  \item When the transformations are computed using quantiles of distribution,
        the final values of the scores are now mapped to the interval 
        \eqn{(1,K)}, if \eqn{K} is the number of levels of a given facrtor.
        The reasons of this choice are explained in the \sQuote{Details}
        section of the \code{smof} documentation page.       
  \item Generally improved and more detailed documentation.      
  }
}        
%---   
\section{Changes in smof version 1.2.1 (2024-10-30)}{
\itemize{
  \item Objects of class \code{gnm} are now handled.
  \item Methods \code{coef} and \code{anova} are introduced for class \code{smof}.
  \item The \code{smof} argument \code{trace} is replaced by the more flexible
        argument \code{verbose}.
  \item New \code{smof} argument \code{f.tail} allows to select the
        suffix being appended to the names the \code{factors} components.       
  \item Fixed a bug in the parsing of the object formula which caused problems
        when some factor has a name which is portion of another factor name.
  \item Fixed a bug in \code{smof_refit} affecting the retrival of the 
        work parameters estimated by a previous call to \code{smof}.       
  }
}
%---   
\section{Changes in smof version 1.2.0 (2024-09-17)}{
\itemize{
  \item A variant method for assigning numeric scores is introduced, based on  
        a monotonic spline instead of the quantile function of a distribution.
  \item An implication of the new variant method is that \code{smof} argument 
        \code{scoring} had to redesigned.
  \item The new function \code{smof_refit} is introduced.    
  \item Additional arguments of \code{smof} are introduced:
        \code{original} and \code{opt.control}.  
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
  \item Scoring of factors are produced from suitably chosen quantiles of 
        distributions among the families \code{"gh", "SU", "sinh-arcsinh"} 
        or their equivalente names.
  \item Admissible classes of fitted models are \code{"lm", "mlm", "glm",
         "survreg", "coxph.penal", "coxph"}.  
  }
}