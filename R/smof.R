#-------------------------------
#  file smof/R/smof.R   
#  This file is a component of the R package 'smof' 
#  copyright (C) 2023-2024 Adelchi Azzalini
#-------------------------------
smof <- function(object, data, factors, scoring, fast.fit=FALSE, original=FALSE,
   opt.control=list(), trace=FALSE)
{    
  target.gen <- function(par, obj, col, ind, K, data, new.formula, scoring, trace=FALSE) {
    # target function for "general" object 
    new.data <- data
    if(scoring$type == "distr") {
      for(j in 1:length(K)) {
        param.j <- c(par[2*j-1], exp(par[2*j]))
        scores <- distr.scores(K[j], scoring$family, param.j)
        if(any(is.infinite(scores) | is.na(scores))) {
          fn.value <- Inf
          if(trace) cat("target.gen par, value:", par, ", ", fn.value, "\n")
          return(fn.value)
          }  
        new.data[, col[j]] <- scores[ind[,j]]
        } 
    } else {
    if(scoring$type == "spline") { 
      M <- scoring$in.knots
      b <- 0
      for(j in 1:length(K)) {
        np <- 2*M[j]
        param.j <- par[(b+1):(b+np)]
        b <- b + np
        scores <- spline.scores(K[j], scoring$method, scoring$in.knots[j], param.j)
        if(any(is.infinite(scores) | is.na(scores))) {
          fn.value <- Inf
          if(trace) cat("target.gen par, value:", par, ", ", fn.value, "\n")
          return(fn.value)
        }  
      new.data[, col[j]] <- scores[ind[,j]]
      }}
    else stop("unknown scoring type")  }
    new.obj <- try(update(object, new.formula, data=new.data))
    if(inherits(new.obj, "try-error")) {
      value <- NA
      if(trace) cat("target.gen par, value:", par, ", ", fn.value, "\n")
      return(value)
      }
    fn.value <- switch(class(new.obj)[1],
      "lm" = deviance(new.obj), # equivalent to (1- summary(new.obj)$r.squared)
      "mlm" = {
         summ <- summary(new.obj)
         out <- 0
         for(k in 1:length(summ)) out <- out + (1 - summ[[k]]$r.squared)
         out},
      "glm" = deviance(new.obj),
      "survreg" = -logLik(new.obj),
      "coxph" = -logLik(new.obj),
      "coxph.penal" = -logLik(new.obj),
      NA)
    if(trace) cat("target.gen par, value:", par, ", ", fn.value, "\n")
    fn.value  
    }  #  # end of target.gen
    
  target.fit <- function(par, obj, col, ind, K, x, family, trace=FALSE)  {
    # target function for objects of class lm or glm, using (g)lm.fit
    for(j in 1:length(K)) {
      param.j <- c(par[2*j-1], exp(par[2*j]))
      scores <- distr.scores(K[j], family, param.j)
      if(any(is.infinite(scores) | is.na(scores))) {
        fn.value <- Inf
        if(trace) cat("target.fit par, value:", par, ", ", fn.value, "\n")
        return(fn.value)
        }  
      x[, col[j]] <- scores[ind[,j]]
      }
    mf <- model.frame(obj)  
    y <- model.response(mf)  
    obj.class <- class(object)[1]
    if(obj.class  == "glm") {
      if(as.character(obj$family[1]) == "binomial" & is.matrix(y)) y <- y[,1]/rowSums(y) 
      fit <- try(glm.fit(x, y, family=obj$family, weights=obj$prior.weights))
      fn.value <- if(inherits(fit, "try-error")) NA else deviance(fit)
      }
    if(obj.class == "lm")  {
      fit <- try(lm.fit(x, y))
      fn.value <- if(inherits(fit, "try-error")) NA else sum(fit$residuals^2)/fit$df.residual
      }
    if(trace) cat("target.fit par, value:", par, ", ", fn.value, "\n")
    fn.value
    } # end of target.fit   
  #---start of smof function body
  cl <- match.call()
  obj.class <- class(object)[1]   
  if(!(obj.class %in% c("lm", "mlm", "glm", "survreg", "coxph.penal", "coxph"))) 
    stop(gettextf("object class '%s' is not (yet) supported", obj.class), domain=NA)
  # obj.call <- object$call
  obj.formula <- formula(object)
  if(mode(factors) != "character") stop("'factors' must be a character vector")
  formula.char <- as.character(obj.formula)
  new.LP.char <- paste(formula.char[c(1,3)], collapse=" ")
  ind <- K <- sof.names <- col.f <- NULL
  nf <- length(factors)
  factors.str <- NULL
  for(j in 1:nf) {
    f.name <- factors[j]
    if(grep(f.name, formula.char[3]) != 1)  
      stop(gettextf("'%s' is not a model component", f.name), domain=NA)
    f <- data[, f.name]
    if(!is.ordered(f))  
      stop(gettextf("'%s' is not an ordered factor", f.name), domain=NA)
    if(length(levels(f)) < 3)
      stop(gettextf("factor '%s' has less then three levels", f.name), domain=NA)
    eval(str2expression(paste("factors.str$", f.name, "<- levels(f)", sep="")))
    K <- c(K, length(levels(f)))
    ind <- cbind(ind, match(f, levels(f)))
    col.factor <-  which(names(data) == f.name)
    sof.name <- paste(f.name, "score", sep=".")
    new.LP.char <- gsub(f.name, sof.name, new.LP.char)
    sof.names <- c(sof.names, sof.name)
    col.f <- c(col.f, col.factor)
    }
  new.formula <- formula(paste(c(formula.char[2], new.LP.char), collapse= " "))    
  new.data <- data
  for(j in 1:nf) {
    new.data[, col.f[j]] <- 0
    names(new.data)[col.f[j]] <- sof.names[j]
    }
  if(missing(scoring)) scoring <- list(type="distr", family="gh")
  stype <- scoring$type
  if(!(stype %in% c("distr", "spline"))) stop("unknown scoring type") 
  if(stype == "spline") scoring$method <- "monoH.FC"
  if(is.null(scoring$param)) {  
    if(stype == "distr") {
       par <- rep(c(0, 0), nf) 
       family <- scoring$family
       }   
    if(stype == "spline") {
       M <- scoring$in.knots
       if(any(is.na(M)) | any(M<1)) stop("invalid scoring$in.knots")
       M <- scoring$in.knots <- as.integer(numeric(nf) + M)    
       par <- rep(0, sum(2*M))
      }
    }  
  else 
    par <- if(stype == "distr") { 
                  tmp <- matrix(scoring$param, nf, 2)
                  tmp[,2] <- log(pmax(tmp[,2], 1e-12))
                  c(t(tmp))
           }  else unlist(scoring$param)
  if(fast.fit) {
    if(!(class(object)[1] %in% c("lm", "glm"))) stop("wrong object class")
    if(stype != "distr") stop("fast.fit requires scoring$type='distr'")
    x <- model.matrix(new.formula, new.data)
    col.x <- NULL
    for(j in 1:nf) col.x <- c(col.x, which(dimnames(x)[[2]] == sof.names[j]))
    opt <- optim(par, target.fit, obj=object, col=col.x, ind=ind, K=K, x=x,
                 scoring=scoring, trace=trace, control=opt.control)     
    } else 
    opt <- optim(par, target.gen, obj=object, col=col.f, ind=ind, K=K, data=new.data, 
          new.formula=new.formula, scoring=scoring, trace=trace, control=opt.control) 
  if(opt$convergence != 0)              
   message("Non-zero error code returned by optim: ", opt$convergence, 
    ". Consider using opt.control or function smof_refit.")
    # "consider restarting optimization, possibly using function smof_refit"))             
  if(trace) {cat("min optim value:", format(opt$value), "\n")  
             cat("optim$par:", format(opt$par), "\n")  }            
  scores <- param <- xf <- NULL
  b <- 0 
  M <- scoring$in.knots
  f.scores <- list()
  param <- knots <- NULL
  for(j in 1:nf) {
    if(stype == "distr") {
      # browser()
      param <- rbind(param, c(opt$par[2*j-1], exp(opt$par[2*j])))
      scores <- distr.scores(K[j], scoring$family, param[j,]) 
      } 
    if(stype == "spline") { 
      np <- 2*M[j]
      param[[j]] <- opt$par[(b+1):(b+np)]
      b <- b + np   
      scores <- spline.scores(K[j], scoring$method, scoring$in.knots[j], param[[j]], knots=TRUE) 
      }
      f.scores[[j]] <-  scores
      names(f.scores)[j] <- paste(factors[j], ".score", sep="")
      new.data[, col.f[j]] <- scores[ind[, j]]
      } 
  if(stype == "distr") {
    dimnames(param)[[1]] <- sof.names    
    dimnames(param)[[2]] <- switch(scoring$family, 
      "SU" = c("gamma", "delta"),
      "gh" = , "g-and-h" = c("g", "h"),
      "sinh-arcsinh" = , "SAS" = c("epsilon", "delta"),
      NULL) }
  if(stype == "spline") {
    names(param) <- sof.names 
    scoring$knots <- knots 
    }
  scoring$param <- param
  new.obj <- update(object, new.formula, data=new.data) 
  out <- list(call=cl, new.object=new.obj, new.data=new.data, scoring=scoring, 
              factors.scores=f.scores, original.factors=factors.str, 
              target.criterion=opt$value, opt=opt) 
  if(original) {object$data <- data; out$original.object <- object}          
  class(out) <- "smof"
  return(out)          
}

 
distr.scores <- function(K, family, param) 
  qTnorm((1:K)/(1+K), family, param)
  
qTnorm <- function(p, family, param) {
  # quantiles of distributions obtained by transformations of normal variates
  gh.funct <- g_and_h.funct <- function(z, g, h) # Tukey g-and-h distribution
                 if(g == 0) z*exp(h*z^2/2) else (exp(g*z)-1)*exp(h*z^2/2)/g  
  SU.funct <- function(z, gamma, delta) sinh((z-gamma)/delta)  # Johnson S_U
  sinh_arcsinh.funct <- function(z, eps, delta) sinh((asinh(z) + eps)/delta)
  norm.scores <- qnorm(p)
  par <- param
  if(par[2] < 0) stop("(param[2] < 0")
  if(par[2] == 0 & family %in% c("SU", "sinh-arcsinh")) 
        stop("(param[2]=0 not compatible with distribution family")
  switch(family, 
        "gh" = gh.funct(norm.scores, par[1], par[2]),
        "g-and-h" = g_and_h.funct(norm.scores, par[1], par[2]),
        "SU" = SU.funct(norm.scores, par[1], par[2]),
        "sinh-arcsinh" = sinh_arcsinh.funct(norm.scores, par[1], par[2]),
        "SAS" = sinh_arcsinh.funct(norm.scores, par[1], par[2]),
        NA)
  }  

spline.scores <- function(K, method, in.knots, param, knots=FALSE) {
  if(method != "monoH.FC") stop("unsuitable spline method")
  M <- in.knots # in.knots: number of internal knots, inside the interval (1, K)
  if(M < 1) stop("M < 1")
  par.x <- param[1:M]
  x.pos <- cumsum(plogis(par.x))
  x.pos <- x.pos/(0.5 + max(x.pos))
  x.knots <- c(1, 1 + (K-1) * x.pos, K)
  # print(rbind(param[1:M], x.pos, knots))
  par.y <- param[(M+1):(2*M)]
  y.pos <- cumsum(plogis(par.y))
  y.pos <- y.pos/(0.5 + max(y.pos))
  y.knots <- c(1, 1 + (K-1) * y.pos, K)
  spl <- splinefun(x.knots, y.knots, method=method)
  out <- spl(1:K)
  if(knots) attr(out, "knots") <- rbind(x=x.knots, y=y.knots)
  out 
  }

smof_refit <- function(object, searches=10, sd=1, opt.control=list(), verbose=1) 
{
  if(!inherits(object, "smof")) stop("object of wrong class")
  new.scoring <- object$scoring
  # if(new.scoring$type != "spline") stop("scoring type must be 'spline'")
  value.best <- object$opt$value
  par.best <- object$opt$par
  object.best <- object
  new.obj <- object$new.object
  improve <- 0
  if(verbose > 0) {
    cat("initial target value:", format(value.best), "\n")
    cat("initial par vector:", format(par.best), "\n")
    }  
  for(j in 1:searches) {  
    if(verbose > 0) cat("search number:", j, "\n")
    new.scoring$param <- par.best + rnorm(length(par.best), 0, sd)
    new.fit <- smof(object$original.object, data=object$original.object$data, 
                  factors=names(object$original.factors), scoring=new.scoring, 
                  opt.control=opt.control, original=TRUE, trace=(verbose>1))
    if(verbose > 0) { 
      cat("target value:", format(new.fit$opt$value), "\n")
      cat("par vector:", format(new.fit$opt$par), "\n")
      }                
    if(new.fit$opt$value < value.best) {
      value.best <- new.fit$opt$value
      par.best <- new.fit$opt$par
      object.best <- new.fit
      improve <- improve + 1
      if(verbose > 0) cat("* new best value *\n")
      }
    }  
  invisible(object.best)  
  }     
