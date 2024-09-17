\name{NEWS}
\title{News for Package \sQuote{smof}}
\encoding{UTF-8}
%---   
\section{Changes in smof version 1.2.0 (2024-09-17)}{
\itemize{
  \item A variant method for assigning numeric scores is introduced, based on 
        monotonic splines.       
  \item An implication of the new variant method is that one argument of the 
        function \code{smof} had to redesigned.
  \item New function \code{smof_refit} is introduced.    
  \item Additional arguments of \code{smof} have been introduced.  
  \item Fixed a bug in \code{predict.smof} which affected the case with multiple factors.
  }
}  
%---   
\section{Changes in smof version 1.1.0 (2024-03-04)}{
\itemize{
  \item The object returned by \code{smof} has now more components
         and is assigned class \code{smof}.
  \item Various methods for class \code{smof} are introduced.
  \item Fixed a bug which prevented the use of option \code{fast.fit=TRUE}
        when \code{family=binomial(..)} and the response was of binary type.
  }
}
%---   
\section{Initial smof version 1.0 (2023-11-21)}{
  \itemize{
  \item First public version of the package and upload to CRAN, with
        only one exported function, namely \code{smof}.
  }
}