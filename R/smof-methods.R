
summary.smof <- function(object, ...) {
  obj <- object
  summ <- obj[match(c("call", "distr", "factors.scores", "original.factors"),
         names(obj), nomatch=0)]
  # s.new.obj <- summary(obj$new.object)
  summ <- c(summ, list(new.object=summary(obj$new.object, ...)))
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
    u <- obj$factors[[j]]
    names(u) <- orig[[j]]
    K <- length(u)
    cat("[Factor ", j, "]... " , names(orig[j]), ", with ", K, " levels:\n", sep="") 
    print(u, digits=digits)
    }
  cat("\n")    
  cat("Scoring via transformation(s) based on distribution:", obj$distr$type, "\n")
  cat("Parameters of the transformation(s) by factor:\n")
  print(obj$distr$param, digits=digits)
  # for(j in 1:nf)   {
  #     cat("  ", names(obj$original.factors)[j], ": ", sep="") 
  #     cat(format(obj$distr$param[j,], digits=digits), "\n")
  #     } 
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
	ask <- ((nf > 1) && dev.interactive())
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
	plot(1:npt, y.lev, xaxt="n", ann=FALSE, ...)
	axis(side=1, at=1:npt, labels=x.lev)
	xlab <- if(void(dots$xlab[j])) f.name else dots$xlab[j]
	ylab <- if(void(dots$ylab[j])) paste("scores of", f.name) else dots$ylab[j]
	main <- if(void(dots$main[j])) paste("scoring of", f.name) else dots$main[j]
	sub <- if(void(dots$sub[j])) NULL else dots$sub[j]
	title(main = main, sub = sub, xlab = xlab, ylab = ylab)
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
  distr <- object$distr
  f.names <- names(f.str)  
  # new.data <- data 
  for(j in 1:length(f.str)) {
    f.name <-  names(f.str)[j]
    if(is.element(f.name, names(data))) {
      f.lev <- f.str[[j]]
      col.j <- which(f.name == names(data))
      new.lev <- levels(data[,col.j]) 
      if(length(setdiff(new.lev, f.lev)) > 0) 
        stop(gettextf("too many levels of factor '%s'", f.name), domain=NA)
      scores <- transform.scores(length(f.lev), distr$type, distr$param)  
      ind <- match(data[,col.j], f.lev)  
      new.data[, col.j] <- scores[ind]
      names(new.data)[col.j] <- paste(f.name, ".score", sep="")
      } else 
      warning(gettextf("there is no factor '%s' in 'data'", f.name), domain=NA)
  }
  new.data    
}
