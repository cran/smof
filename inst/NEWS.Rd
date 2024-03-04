\name{NEWS}
\title{News for Package \sQuote{smof}}
\encoding{UTF-8}
\section{Changes in smof version 1.1.0 (2024-03-04)}{
\itemize{
  \item The object returned by \code{smof} has now more components
         and is assigned class \code{smof}.
  \item Various methods for class \code{smof} are introduced.
  \item Fixed a bug which prevented the use of option \code{fast.fit=TRUE}
        when \code{family=binomial(..)} and the response was of binary type.
  }
}
\section{Initial smof version 1.0 (2023-11-21)}{
  \itemize{
  \item First public version of the package and upload to CRAN, with
        only one exported function, namely \code{smof}.
  }
}