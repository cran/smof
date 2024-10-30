#-------------------------------
#  file smof/R/smof-methods.R   
#  This file is a component of the R package 'smof' 
#  copyright (C) 2023-2024 Adelchi Azzalini
#-------------------------------
#
summary.smof <- function(object, ...) {
  obj <- object
  summ <- obj[match(c("call", "scoring", "factors.scores", "original.factors", 
             "target.criterion"), names(obj), nomatch=0)]
  summ <- c(summ,  list(new.object=summary(obj$new.object, ...)))   
  class(summ) <- "summary.smof"
  summ
  }

print.summary.smof <- function(x, ... ) {
  obj <- x
  dots <- list(...)
  digits <- if(!is.null(dots$digits)) dots$digits else getOption("digits")
  cat("Call:  "); print(obj$call); cat("\n")
  cat("Final fitting call:  ")
  print(obj$new.object$call)
  cat("\n")
  nf <- length(obj$factors)
  cat("Number of ordered factors processed by 'smof':", nf, "\n")
  cat("Ordered factor(s) and numeric scores assigned to the levels:\n")
  orig <- obj$original.factors
  for(j in 1:nf)   {
    u <- c(obj$factors[[j]])
    names(u) <- orig[[j]]
    K <- length(u)
    cat("\n[Factor ", j, "]... " , names(orig[j]), ", with ", K, " levels:\n", sep="") 
    print(u, digits=digits)
    }
  cat("\n")   
  obscor <- obj$scoring 
  cat("Scores obtained by transformation(s) of type:", obscor$type, "\n")
  if(obscor$type == "distr") cat("Family of distributions:", obscor$family, "\n")
  if(obscor$type == "spline") cat("Spline method:", "monoH.FC", "\n")
  cat("Target criterion for selection of the transformation(s):", 
     format(obj$target.criterion, nsmall=2), "\n")
  param <- obscor$param
  if(obscor$type == "distr") { 
    cat("Parameters of the transformation(s) by factor:\n")
    print(param, digits=digits)
    }
  if(obscor$type == "spline")  {
    cat("Splines and parameters of the transformation(s) by factor follow.\n")
    for(j in 1:nf) {
      #if(obj$scoring$type == "spline") {
      knots <- attr(obj$factors[[j]], "knots")
      cat("\nKnots for factor", names(param[j]), ":\n")
      # cat("knots.x:", format(u[1,], digits=digits), "\n")
      # cat("knots.y:", format(u[2,], digits=digits), "\n")
      print(knots)
      #}
    cat("working parameters:", format(param[j], digits=digits), "\n")  
    }}
  if(inherits(obj, "summary.smof"))
    {cat("\nFitted model using factor(s) scores:\n"); print(obj$new.object) } 
  invisible(x)  
}

print.smof <- function(x, ...) {
  print.summary.smof(x, ...)
  invisible(NULL)  
}

plot.smof <- function(x, which, ...) {
  obj <- x
  dots <- list(...)
  f.names <- names(obj$original.factors)
  nf <- length(f.names)
  if(missing(which)) which <- seq_len(nf) 
  if(is.list(which)) {w1 <- which[[1]];  w2 <- which[[2]]} 
    else {w1 <- which; w2 <- NULL}
  if(is.character(w1)) w1 <- match(intersect(w1, f.names), f.names)
  if(is.numeric(w1))  w1 <- intersect(w1, 1:nf)
  if(length(w1) > 0) { 
	ask <- ((length(w1) > 1) && interactive())
	# localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
	void <- function(x) if(is.null(x)) TRUE else is.na(x)  
	for(j in w1) {
	if (ask) {
	  oask <- devAskNewPage(TRUE) 
	  on.exit(devAskNewPage(oask))
	}
	f.name <- f.names[j]
	x.lev <- obj$original.factors[[j]]
	npt <- length(x.lev)
	y.lev <- obj$factors.scores[[j]]
	graphics::plot.default(1:npt, y.lev, xaxt="n", ann=FALSE, ...)
	graphics::axis(side=1, at=1:npt, labels=x.lev)
	xlab <- if(void(dots$xlab[j])) f.name else dots$xlab[j]
	ylab <- if(void(dots$ylab[j])) paste("scores of", f.name) else dots$ylab[j]
	main <- if(void(dots$main[j])) paste("scoring of", f.name) else dots$main[j]
	sub <- if(void(dots$sub[j])) NULL else dots$sub[j]
	graphics::title(main = main, sub = sub, xlab = xlab, ylab = ylab)
	}
  }
  if(!is.null(w2)) { 
    if(inherits(obj$new.object, "lm")) plot(obj$new.object, which=w2, ...)	  
    }
  invisible(NULL)
}


predict.smof <- function(object, newdata=NULL, ...) {
  if(!inherits(object, "smof")) stop("object of wrong class")
  if(is.null(newdata)) return(object$new.data)
  if(!inherits(newdata, "data.frame")) stop("'newdata' must be a data.frame")
  new.data <- data <- newdata
  f.str <- object$original.factors 
  nf <- length(f.str)
  f.names <- names(f.str)  
  scoring <- object$scoring
  param <-  scoring$param 
  for(j in 1:nf) {
    f.name <-  names(f.str)[j]
    if(is.element(f.name, names(data))) {
      f.lev <- f.str[[j]]
      col.j <- which(f.name == names(data))
      new.lev <- levels(data[,col.j]) 
      if(length(setdiff(new.lev, f.lev)) > 0) 
        stop(gettextf("too many levels of factor '%s'", f.name), domain=NA)
      scores <- switch(scoring$type,
        distr = distr.scores(length(f.lev), scoring$family, param[j,]),
        spline = spline.scores(length(f.lev), "monoH.FC", scoring$in.knots[j], param[[j]]),
        NULL 
        )
      ind <- match(data[,col.j], f.lev)  
      new.data[, col.j] <- scores[ind]
      names(new.data)[col.j] <- paste(f.name, ".score", sep="")
      } else 
      warning(gettextf("there is no factor '%s' in 'data'", f.name), domain=NA)
  }
  new.data    
}

anova.smof <- function(object, ...) anova(object$new.object)

coef.smof <- function(object, complete=FALSE, ...) {
  if(!inherits(object, "smof")) stop("object of wrong class")
  par <- object$scoring$param
  if(object$scoring$type == "distr") { 
    out <- c(t(par))
    names(out) <- c(t(outer(dimnames(par)[[1]], dimnames(par)[[2]] , FUN=paste, sep=".")))
    } else  
    out <- unlist(par)
  if(complete) {
     out2 <- coef(object$new.object, ...)
     out <- c(out, out2)  
     }
  out
}