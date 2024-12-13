#-------------------------------
#  file smof/R/smof.R   
#  This file is a component of the R package 'smof' 
#  copyright (C) 2023-2024 Adelchi Azzalini
#-------------------------------
smof <- function(object, data, factors, scoring, fast.fit=FALSE, original=FALSE, 
      f.tail=".score", opt.method="Nelder-Mead", opt.control=list(), verbose=0)
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
          if(trace) cat("target.gen, work.par=", par, ", value=", fn.value, "\n")
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
          if(trace) cat("target.gen: work.par=", par, ",  value=", fn.value, "\n")
          return(fn.value)
        }  
      new.data[, col[j]] <- scores[ind[,j]]
      }}
    else stop("unknown scoring type")  } 
    new.obj <- try(update(object, new.formula, data=new.data))
    if(inherits(new.obj, "try-error")) {
      fn.value <- NA
      if(trace) cat("target.gen: work.par=", par, ", value=", fn.value, "\n")
      return(fn.value)
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
      "gnm" = deviance(new.obj),
      NA)
    if(trace) cat("target.gen: work.par=", par, ", value=", fn.value, "\n")
    fn.value  
    }  #  # end of target.gen
    
  target.fit <- function(par, obj, col, ind, K, x, family, trace=FALSE)  {
    # target function for objects of class lm or glm, using (g)lm.fit
    for(j in 1:length(K)) {
      param.j <- c(par[2*j-1], exp(par[2*j]))
      scores <- distr.scores(K[j], family, param.j)
      if(any(is.infinite(scores) | is.na(scores))) {
        fn.value <- Inf
        if(trace) cat("target.fit: work.par=", par, ", value=", fn.value, "\n")
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
    if(trace) cat("target.fit: work.par=", par, ", value=", fn.value, "\n")
    fn.value
    } # end of target.fit   
  #---start of smof function body
  cl <- match.call()
  verbose <- as.integer(verbose)
  if(verbose < 0) stop("negative 'verbose'")
  trace <- (verbose > 1)
  if(!inherits(data, "data.frame")) stop("'data' must be a data.frame")
  obj.class <- class(object)[1]   
  if(!(obj.class %in% c("lm", "mlm", "glm", "survreg", "coxph.penal", "coxph", "gnm"))) 
    stop(gettextf("object class '%s' is not (yet) supported", obj.class), domain=NA)
  basic.methods <- c("Nelder-Mead", "BFGS", "nlminb")
  nloptr.methods <- c("newuoa", "bobyqa", "cobyla", "sbplx")
#   if(requireNamespace("nloptr", quietly = TRUE)) {    
#     if(opt.method %in% nloptr.methods)
#        nloptr.method <- get(opt.method, envir = loadNamespace("nloptr"))
#   } 
#   else
#     nloptr.method <- function(...) stop("this optimization needs the nloptr package.")
  if(opt.method %in% nloptr.methods) {
     if(requireNamespace("nloptr", quietly = TRUE))   
       nloptr.method <- get(opt.method, envir = loadNamespace("nloptr"))
     else 
       stop(gettextf("optimization using %s needs the nloptr package.", 
            opt.method), domain=NA)   
  }
  if(!(opt.method %in% c(basic.methods, nloptr.methods))) 
     { warning(gettextf("opt.method '%s' is not supported, replaced by 'Nelder-Mead'",
              opt.method), domain = NA)
      opt.method <- 'Nelder-Mead'}
  obj.formula <- formula(object)
  if(mode(factors) != "character") stop("'factors' must be a character vector")
  if(any(duplicated(factors))) stop("duplicated factors")
  new.rhs <- rhs <- deparse(obj.formula[[3]]) 
  ind <- K <- sof.names <- col.f <- NULL
  nf <- length(factors)
  factors.str <- NULL
  for(j in 1:nf) {
    f.name <- factors[j]  
    if(grep(f.name, rhs) != 1)  
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
    sof.name <- make.names(paste(f.name, f.tail, sep=""))
    pattern <- paste(f.name, "{1}[^._A-Za-z0-9]|", f.name, "{1}$", sep= "")
    f.x <- c(gregexpr(pattern, new.rhs)[[1]])
    nx <- length(f.x)
    for(i in 1:nx) {
      leftside <- if(f.x[i]==1) NULL else substr(new.rhs, 1,  f.x[i]-1)
      rightside <- substr(new.rhs, f.x[i] + nchar(f.name), nchar(new.rhs))
      new.rhs <- paste(leftside, sof.name, rightside, sep="")
      if(i < nx) f.x[(i+1):nx] <- f.x[(i+1):nx] + nchar(f.tail)
      }
    sof.names <- c(sof.names, sof.name)
    col.f <- c(col.f, col.factor)
    }
  new.formula <- formula(paste(deparse(obj.formula[[2]]), deparse(obj.formula[[1]]),
       new.rhs), collapse= " ")
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
    {tmp <- scoring$param 
    par <- if(stype == "distr") {
             if(!all(dim(tmp) == c(nf,2))) stop("wrong 'scoring$param'") 
             tmp[,2] <- log(pmax(scoring$param[,2], 1e-12))
             c(t(tmp))
           }  else {
             if(!is.list(tmp) | length(tmp) != nf) stop("wrong 'scoring$param'") 
             unlist(tmp)
           }
     }      
  if(fast.fit) {
    if(!(class(object)[1] %in% c("lm", "glm"))) stop("wrong object class")
    if(stype != "distr") stop("fast.fit requires scoring$type='distr'")
    x <- model.matrix(new.formula, new.data)
    col.x <- NULL
    for(j in 1:nf) col.x <- c(col.x, which(dimnames(x)[[2]] == sof.names[j]))
    opt <- optim(par, fn=target.fit, method="Nelder-Mead", obj=object, 
                 col=col.x, ind=ind, K=K, x=x,
                 scoring=scoring, trace=trace, control=opt.control)     
    } 
  else {
    if(opt.method %in% c("Nelder-Mead", "BFGS") ) {
      opt <- optim(par, fn=target.gen, gr=NULL, method=opt.method, obj=object, 
               col=col.f, ind=ind, K=K, data=new.data, new.formula=new.formula,
               scoring=scoring, trace=trace, control=opt.control) 
      target.value <- opt$value  
      }  
    if(opt.method == "nlminb") {
          opt <- nlminb(par, objective=target.gen,  obj=object, 
               col=col.f, ind=ind, K=K, data=new.data, new.formula=new.formula,
               scoring=scoring, trace=trace, control=opt.control)
         target.value <- opt$objective
         }  
     if(opt.method %in% nloptr.methods) {    
        opt <- nloptr.method(x0=par, fn=target.gen,  obj=object, 
               col=col.f, ind=ind, K=K, data=new.data, new.formula=new.formula,
               scoring=scoring, trace=trace, control=opt.control)
        target.value <- opt$value    
        }
    }   
  bad.conv <- (((opt.method %in% basic.methods) & opt$convergence != 0) |
              ((opt.method %in% nloptr.methods) & ((opt$convergence < 0) | 
                 (opt$convergence > 4))))
  if(verbose > 0 | bad.conv )   {
    if(bad.conv) warning("Possibly unsatisfacory outcome from optimization function.")
    cat("Exit from optimization function '", opt.method, "'\n", sep="")
    if(!is.null(opt$message)) cat("with message '", opt$message, "'\n", sep="")
    if(!is.null(opt$convergence)) cat("convergence code: ", opt$convergence, "\n")
    cat("value of the target criterion:", format(target.value, nsmall=2), "\n")  
    cat("work parameters:", format(opt$par), "\n")

    }              
  scores <- param <- xf <- NULL
  b <- 0 
  M <- scoring$in.knots
  if(is.null(scoring$mapping)) scoring$mapping <- "(1,K)"
  if(!(scoring$mapping %in% c("none", "(1,K)"))) stop('illegal scoring$mapping value')
  f.scores <- list()
  param <- knots <- NULL
  for(j in 1:nf) {
    if(stype == "distr") {
      param <- rbind(param, c(opt$par[2*j-1], exp(opt$par[2*j])))
      scores <- distr.scores(K[j], scoring$family, param[j,]) 
      if(scoring$mapping == "(1,K)") 
        scores <- 1 + (K[j] - 1)*(scores -min(scores))/diff(range(scores))
      } 
    if(stype == "spline") { 
      np <- 2*M[j]
      param[[j]] <- opt$par[(b+1):(b+np)]
      b <- b + np   
      scores <- spline.scores(K[j], scoring$method, scoring$in.knots[j], param[[j]], knots=TRUE) 
      }
    f.scores[[j]] <-  scores
    names(f.scores)[j] <- paste(factors[j], f.tail, sep="")
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
  opt. <- c(list(opt.method=opt.method), opt)
  out <- list(call=cl, new.object=new.obj, new.data=new.data, scoring=scoring, 
              factors.scores=f.scores, original.factors=factors.str, 
              target.value=target.value, opt=opt.) 
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
  if(anyNA(par)) return(rep(NA, length(p)))
  if(par[2] < 0) stop("(param[2] < 0")
  if(par[2] == 0 & family %in% c("SU", "sinh-arcsinh")) 
        stop("(param[2]=0 not compatible with family of distributions")
  switch(family, 
        "gh" = gh.funct(norm.scores, par[1], par[2]),
        "g-and-h" = g_and_h.funct(norm.scores, par[1], par[2]),
        "SU" = SU.funct(norm.scores, par[1], par[2]),
        "sinh-arcsinh" = sinh_arcsinh.funct(norm.scores, par[1], par[2]),
        "SAS" = sinh_arcsinh.funct(norm.scores, par[1], par[2]),
        NA)
  }  

spline.scores <- function(K, method, in.knots, param, knots=FALSE) {
  if(missing(method)) method <- "monoH.FC"
  if(method != "monoH.FC") stop("unsuitable spline method")
  incr <- 0.5
  M <- in.knots   # number of internal knots, inside the interval (1, K)
  if(M < 1) stop("M < 1")
  if(M > K-2) stop("M exceeds K-2")
  par.x <- param[1:M]
  x.pos <- cumsum(plogis(par.x))
  x.pos <- x.pos/(incr + max(x.pos))
  x.knots <- c(1, 1 + (K-1) * x.pos, K)
  par.y <- param[(M+1):(2*M)]
  y.pos <- cumsum(plogis(par.y))
  y.pos <- y.pos/(incr+ max(y.pos))
  y.knots <- c(1, 1 + (K-1) * y.pos, K)
  spl <- splinefun(x.knots, y.knots, method=method)
  out <- spl(1:K)
  if(knots) attr(out, "knots") <- rbind(x=x.knots, y=y.knots)
  out 
  }

smof_refit <- function(object, searches=10, sd=1, opt.control=list(), verbose=1) 
{
  if(!inherits(object, "smof")) stop("object of wrong class")
  best.object <- object
  best.value <- object$target.value
  best.par <- object$scoring$param
  new.obj <- object$new.object
  new.scoring <- object$scoring
  nf <- length(object$factors.scores)
  searches <- as.integer(searches)
  if(searches < 1) stop("meaningless number of 'searches'")
  improve <- 0
  if(verbose > 0) 
    cat("initial target value:", format(best.value, nsmall=2), "\n")
  new.par <- best.par
  # basic.methods <- c("Nelder-Mead", "BFGS", "nlminb")
  for(search in 1:searches) {  
    if(verbose > 0) cat(">> search number:", search, "\n")
    if(object$scoring$type == "distr") {
       new.par[,1] <- best.par[,1] + rnorm(nf, 0, sd)
       new.par[,2] <- best.par[,2] * exp(rnorm(nf, 0, sd) - 1)
       }  
    if(object$scoring$type == "spline") 
      for(j in 1:nf) new.par[[j]] <- new.par[[j]] + rnorm(length(new.par[[j]]), 0, sd)
    new.scoring$param  <- new.par 
    # opt.method <- basic.methods[1 + (search-1) %% 3]
    if(verbose > 0) {
      cat("param vector:", format(new.scoring$param), "\n")   
      # cat("opt.method:", opt.method,"\n")
      } 
    new.fit <- try(smof(object$original.object, data=object$original.object$data, 
                  factors=names(object$original.factors), scoring=new.scoring, 
                  opt.method=object$opt$opt.method, 
                  opt.control=opt.control, original=TRUE, verbose=(verbose-1)))
    if(inherits(new.fit, "try-error")) 
      { if(verbose > 0) cat("Problems with current search. Skip it.\n");  next}
    new.value <- new.fit$target.value       
    if(verbose > 0) 
      cat("current target value:", format(new.value, nsmall=2), "\n")
    if(new.fit$target.value < best.value) {
      best.object <- new.fit
      best.value <-  new.value
      best.par <-  best.object$scoring$param
      improve <- improve + 1
      if(verbose > 0) cat("** new best value **\n")
      }
    }  
  if(verbose > 0) {
    cat(">> End of search.\n")
    if(improve == 0) cat("No improvement of the target criterion\n")
    if(improve > 0) cat("Target criterion has improved", improve, "time(s).",
      "Best target value is", format(best.value, nsmall=2), "\n")
    }  
  invisible(best.object)  
  }     
