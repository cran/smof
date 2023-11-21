#  file smof/R/smof.R   
#  This file is a component of the R package 'smof' 
#  copyright (C) 2023 Adelchi Azzalini

smof <-
function(object, data, factors, distr.type="gh", fast.fit=FALSE, trace=FALSE)
{
  transform.scores <- function(K, distr.type, param) {
    gh.funct <- g_and_h.funct <- function(z, g, h) # Tukey g-and-h distribution
                   if(g == 0) z*exp(h*z^2/2) else (exp(g*z)-1)*exp(h*z^2/2)/g  
    SU.funct <- function(z, gamma, delta) sinh((z-gamma)/delta)  # Johnson S_U
    sinh_arcsinh.funct <- function(z, eps, delta) sinh((asinh(z) + eps)/delta)
    p <- param
    if(p[2] < 0) stop("(param[2] < 0")
    if(p[2] == 0 & distr.type %in% c("SU", "sinh-arcsinh")) 
          stop("(param[2]=0 not compatible with distr.type")
    norm.scores <- qnorm((1:K)/(1+K))
    switch(distr.type, 
          "gh" = gh.funct(norm.scores, p[1], p[2]),
          "g-and-h" = g_and_h.funct(norm.scores, p[1], p[2]),
          "SU" = SU.funct(norm.scores, p[1], p[2]),
          "sinh-arcsinh" = sinh_arcsinh.funct(norm.scores, p[1], p[2]),
          "SAS" = sinh_arcsinh.funct(norm.scores, p[1], p[2]),
          NA)
    }
    
  target.gen <- function(par, obj, col, ind, K, data, new.formula, distr.type, trace=FALSE) {
    # target function for "general" object 
    new.data <- data
    for(j in 1:length(K)) {
      param.j <- c(par[2*j-1], exp(par[2*j]))
      scores <- transform.scores(K[j], distr.type, param.j)
      if(any(is.infinite(scores) | is.na(scores))) {
        fn.value <- Inf
        if(trace) cat("target.gen par, value:", par, ", ", fn.value, "\n")
        return(fn.value)
        }  
      new.data[, col[j]] <- scores[ind[,j]]
      }
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
    
  target.fit <- function(par, obj, col, ind, K, x, distr.type, trace=FALSE)  {
    # target function for objects of class lm or glm, using (g)lm.fit
    for(j in 1:length(K)) {
      param.j <- c(par[2*j-1], exp(par[2*j]))
      scores <- transform.scores(K[j], distr.type, param.j)
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
      if(as.character(obj$family[1]) == "binomial") y <- y[,1]/rowSums(y)
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
  obj.class <- class(object)[1]   
  if(!(obj.class %in% c("lm", "mlm", "glm", "survreg", "coxph.penal", "coxph"))) 
    stop(gettextf("object class '%s' is not (yet) supported", obj.class), domain=NA)
  obj.call <- object$call
  obj.formula <- formula(object)
  if(mode(factors) != "character") stop("'factors' must be a character vector")
  formula.char <- as.character(obj.formula)
  new.LP.char <- paste(formula.char[c(1,3)], collapse=" ")
  ind <- K <- sof.names <- col.f <- NULL
  nf <- length(factors)
  for(j in 1:nf) {
    factor <- factors[j]
    f <- data[, factor]
    if(!is.ordered(f))  
      stop(gettextf("'%s' is not an ordered factor", factor), domain=NA)
    K <- c(K, length(levels(f)))
    ind <- cbind(ind, match(f, levels(f)))
    col.factor <-  which(names(data) == factor)
    sof.name <- paste(factor, "score", sep=".")
    new.LP.char <- gsub(factor, sof.name, new.LP.char)
    sof.names <- c(sof.names, sof.name)
    col.f <- c(col.f, col.factor)
    }
  new.formula <- formula(paste(c(formula.char[2], new.LP.char), collapse= " "))    
  new.data <- data
  for(j in 1:nf) {
    new.data[, col.f[j]] <- 0
    names(new.data)[col.f[j]] <- sof.names[j]
    }
  
  par <- rep(c(0, 0), nf)  
  if(fast.fit) {
    if(!(class(object)[1] %in% c("lm", "glm"))) stop("wrong object class")
    x <- model.matrix(new.formula, new.data)
    col.x <- NULL
    for(j in 1:nf) col.x <- c(col.x, which(dimnames(x)[[2]] == sof.names[j]))
    opt <- optim(par, target.fit, obj=object, col=col.x, ind=ind, K=K, x=x,
                 distr.type=distr.type, trace=trace)     
    } 
  else 
    opt <- optim(par, target.gen, obj=object, col=col.f, ind=ind, K=K, data=new.data, 
                new.formula=new.formula, distr.type=distr.type, trace=trace)      
  scores <- param <- xf <- NULL
  f.scores <- list()
  for(j in 1:nf) {
      param <- cbind(param, c(opt$par[2*j-1], exp(opt$par[2*j])))
      scores <- transform.scores(K[j], distr.type, param[,j])
      f.scores[[j]] <-  scores
      names(f.scores)[j] <- paste(factors[j], ".score", sep="")
      new.data[, col.f[j]] <- scores[ind[, j]]
      } 
  dimnames(param)[[2]] <- sof.names    
  distr <- list(type=distr.type, param=t(param))
  new.obj <- update(object, new.formula, data=new.data) 
  list(object=new.obj, data=new.data, distr=distr, factors.scores=f.scores)   
}
