library(laGP)
library(MASS)

# setOldClass("laGP")

# matching the covariance structure
setClass("covlaGP2covkm", 
         slots = list(sd2 = "numeric", covtype = "character",
                      theta = "numeric", nugget.flag = "logical"))

covlaGP2covkm <- function(model){
  return(new(Class = "covlaGP2covkm", sd2 = model$phi, covtype = "Gaussian",
             theta = model$d$start, nugget.flag = FALSE))
}

if(!isGeneric("covMat1Mat2")) {
  setGeneric(name    = "covMat1Mat2",
             def     = function(object, X1, X2, nugget.flag=FALSE) standardGeneric("covMat1Mat2")
  )
}

covMat1Mat2 <- function (object, X1, X2, nugget.flag) {
  UseMethod("covMat1Mat2", object, X1, X2, nugget.flag)
}

covMat1Mat2.covlaGP2covkm <- function(object, X1, X2, nugget.flag=FALSE){
  return(exp(-distance(X1, X2) / object@theta))
}

setMethod("covMat1Mat2", "covlaGP2covkm", 
          definition = function(object, X1, X2, nugget.flag = FALSE){
            covMat1Mat2.covlaGP2covkm(object, X1, X2, nugget.flag)
          }
)

# setOldClass("prelaGP")
setClass("laGP2km", slots = list(X = "matrix", y = "matrix", d = "integer", n = "integer", start = "integer", end = "integer",
                                 noise.var = "numeric", noise.flag = "logical", method = "character", da = "list", ga = "list",
                                 covariance = "covlaGP2covkm", trend.formula = "formula", trend.coef = "numeric",
                                 M = "matrix", T = "matrix", F = "matrix", z = "numeric"))

# newmod <- new(Class = "homGP2km", model = mhom1, X = mhom1$X0)

# model: pre-laGP object containing mle (with both d and g), X, Z, phi, start, end fields (basically all what is needed to call laGP)
# phi <- crossprod(Z, chol2inv(T)) %*% Z

laGP2km <- function(model){
  ## These are not relevant (or very costly for large X, Z)
  M <- matrix(NA); T <- matrix(NA); F <- matrix(NA); m <- matrix(NA); z <- c(0)
  
  res <- new(Class = "laGP2km", X = model$X, y = matrix(model$Z, ncol = 1),
             d = ncol(model$X), n = nrow(model$X), noise.var = model$g$start * model$phi, noise.flag = TRUE,
             covariance = covlaGP2covkm(model), trend.formula = ~1, trend.coef = 0, method = model$method,
             start = model$start, end = model$end, da = model$da, ga = model$ga, 
             F = F, M = M, T = T, z = z
  )
  return(res)
}

setMethod("predict", "laGP2km", function(object, newdata, type, se.compute = TRUE, 
                                         cov.compute = FALSE, light.return = FALSE,
                                         bias.correct = FALSE, checkNames = TRUE, ...){
  if(class(newdata) == "data.frame") newdata <- as.matrix(newdata)
  
  preds <- laGP(Xref = newdata, start = object@start, end = object@end, X = object@X, Z = object@y,
                d = object@da, g = object@ga, method = object@method, lite = !cov.compute)

  # Ensure positivity of s2
  preds$s2 <- pmax(0, preds$s2)
  
  if(cov.compute){
    sd <- sqrt(diag(preds$s2))
    cov <- preds$s2
  }else{
    sd <- sqrt(preds$s2)
    cov <- NULL
  }
  # if(!is.null(xprime)) preds$cov <- 1/2*(preds$cov + t(preds$cov)) # to ensure symmetry
  return(list(mean = preds$mean, sd = sd, cov = cov))
})


setMethod("update", "laGP2km", function(object, newX, newy, newX.alreadyExist = FALSE,
                                        cov.reestim = TRUE, trend.reestim = TRUE, nugget.reestim = FALSE, 
                                        newnoise.var = NULL, kmcontrol = NULL, newF = NULL,...){
  # res <- update(object@model, Xnew = newX, Znew = newy)
  # return(homGP2km(res))
  object@X <- rbind(object@X, newX)
  object@y <- rbind(object@y, as.numeric((newy)))
  object@n <- object@n + length(newy)
  object@da <- darg(object@da, object@X)
  object@ga <- darg(object@ga, object@y)
  
  return(object)
})


setMethod("simulate", "laGP2km", function(object, nsim=1, seed=NULL, newdata=NULL,
                                          cond=FALSE, nugget.sim=0, checkNames=TRUE, ...){
  if(cond==FALSE){
    mu <- rep(0, nrow(newdata))
    K <-  exp(-distance(newdata)/ model@d$start) + diag(model@g$start, nrow(newdata))
  }else{
    preds <- laGP(Xref = newdata, start = object@start, end = object@end, X = object@X, Z = object@y, 
                  d = object@da, g = object@ga, method = object@method, lite = FALSE) 
    mu <- preds$mean
    K <- preds$se
  }
  y <- mvrnorm(n = nsim, mu = mu, Sigma = K)
})

################################################################################
### Same for laGPsep

# matching the covariance structure
setClass("covlaGPsep2covkm", 
         slots = list(sd2 = "numeric", covtype = "character",
                      theta = "numeric", nugget.flag = "logical"))

covlaGPsep2covkm <- function(model){
  return(new(Class = "covlaGPsep2covkm", sd2 = model$phi, covtype = "Gaussian",
             theta = model$d$start, nugget.flag = FALSE))
}

covMat1Mat2.covlaGsepP2covkm <- function(object, X1, X2, nugget.flag=FALSE){
  return(exp(-distance(X1, X2) %*% diag(1/object@theta)))
}

setMethod("covMat1Mat2", "covlaGPsep2covkm", 
          definition = function(object, X1, X2, nugget.flag = FALSE){
            covMat1Mat2.covlaGPsep2covkm(object, X1, X2, nugget.flag)
          }
)

setClass("laGPsep2km", slots = list(X = "matrix", y = "matrix", d = "integer", n = "integer", start = "integer", end = "integer",
                                 noise.var = "numeric", noise.flag = "logical", method = "character", da = "list", ga = "list",
                                 covariance = "covlaGPsep2covkm", trend.formula = "formula", trend.coef = "numeric",
                                 M = "matrix", T = "matrix", F = "matrix", z = "numeric"))

# model: pre-laGP object containing mle (with both d and g), X, Z, start, end fields (basically all what is needed to call laGP)

laGPsep2km <- function(model){
  ## These are not relevant (or very costly for large X, Z)
  M <- matrix(NA); T <- matrix(NA); F <- matrix(NA); m <- matrix(NA); z <- c(0)
  
  res <- new(Class = "laGPsep2km", X = model$X, y = matrix(model$Z, ncol = 1),
             d = ncol(model$X), n = nrow(model$X), noise.var = model$g$start * model$phi, noise.flag = TRUE,
             covariance = covlaGPsep2covkm(model), trend.formula = ~1, trend.coef = 0, method = model$method,
             start = model$start, end = model$end, da = model$da, ga = model$ga, 
             F = F, M = M, T = T, z = z
  )
  return(res)
}

setMethod("predict", "laGPsep2km", function(object, newdata, type, se.compute = TRUE, 
                                         cov.compute = FALSE, light.return = FALSE,
                                         bias.correct = FALSE, checkNames = TRUE, ...){
  if(class(newdata) == "data.frame") newdata <- as.matrix(newdata)
  
  preds <- laGPsep(Xref = newdata, start = object@start, end = object@end, X = object@X, Z = object@y,
                d = object@da, g = object@ga, method = object@method, lite = !cov.compute)
  
  # Ensure positivity of s2
  preds$s2 <- pmax(0, preds$s2)
  
  if(cov.compute){
    sd <- sqrt(diag(preds$s2))
    cov <- preds$s2
  }else{
    sd <- sqrt(preds$s2)
    cov <- NULL
  }
  # if(!is.null(xprime)) preds$cov <- 1/2*(preds$cov + t(preds$cov)) # to ensure symmetry
  return(list(mean = preds$mean, sd = sd, cov = cov))
})


setMethod("update", "laGPsep2km", function(object, newX, newy, newX.alreadyExist = FALSE,
                                        cov.reestim = TRUE, trend.reestim = TRUE, nugget.reestim = FALSE, 
                                        newnoise.var = NULL, kmcontrol = NULL, newF = NULL,...){
  # res <- update(object@model, Xnew = newX, Znew = newy)
  # return(homGP2km(res))
  object@X <- rbind(object@X, newX)
  object@y <- rbind(object@y, as.numeric((newy)))
  object@n <- object@n + length(newy)
  object@da <- darg(object@da, object@X)
  object@ga <- darg(object@ga, object@y)
  
  return(object)
})


setMethod("simulate", "laGPsep2km", function(object, nsim=1, seed=NULL, newdata=NULL,
                                          cond=FALSE, nugget.sim=0, checkNames=TRUE, ...){
  if(cond==FALSE){
    mu <- rep(0, nrow(newdata))
    K <-  exp(-distance(newdata) %*% diag(1 / model@d$start)) + diag(model@g$start, nrow(newdata))
  }else{
    preds <- laGPsep(Xref = newdata, start = object@start, end = object@end, X = object@X, Z = object@y, 
                  d = object@da, g = object@ga, method = object@method, lite = FALSE) 
    mu <- preds$mean
    K <- preds$se
  }
  y <- mvrnorm(n = nsim, mu = mu, Sigma = K)
})

## Examples

## with GPareto
library(GPareto)
library(DiceDesign)

set.seed(25468)

d <- 2 
fname <- P1
n.grid <- 21
test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
test.vals <- t(apply(test.grid, 1, fname))
nappr <- 15
design.grid <- maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design
response.grid <- t(apply(design.grid, 1, fname))
Front_Pareto <- t(nondominated_points(t(response.grid)))

mf1 <- km(~1, design = design.grid, response = response.grid[, 1], lower=c(.1,.1))
mf2 <- km(~., design = design.grid, response = response.grid[, 2], lower=c(.1,.1))
model <- list(mf1, mf2)

nsteps <- 10
lower <- rep(0, d)
upper <- rep(1, d)

# Optimization 1: EHI with pso
optimcontrol <- list(method = "pso", maxit = 20)
critcontrol <- list(refPoint = c(1, 10))
omEGO1 <- GParetoptim(model = model, fn = fname, crit = "EHI", nsteps = nsteps,
                      lower = lower, upper = upper, critcontrol = critcontrol,
                      optimcontrol = optimcontrol)

# Same with laGP

prelaGP <- function(X, Z, start, end, phi = 1, method = "nn", d = NULL, g = NULL, ...){
  if(is.null(d)) d <- darg(NULL, X)
  if(is.null(g)) g <- garg(list(mle = TRUE), Z)
  if(length(g) == 1) g <- list(mle = FALSE, start = g, max = g, min = sqrt(.Machine$double.eps), ab = c(0, 0))
  
  phi <- 1 # Global variance is not relevant here
  res <- list(X = X, Z = Z, da = d, ga = g, phi = phi, 
              start = as.integer(start), end = as.integer(end), method = method)
  class(res) <- 'prelaGP'
  return(res)
}

mla1 <- prelaGP(design.grid, response.grid[,1], start = nappr-5, end = nappr-1, g = 1e-7)
mla2 <- prelaGP(design.grid, response.grid[,2], start = nappr-5, end = nappr-1, g = 1e-7)


laaskm1 <- laGP2km(mla1)
laaskm2 <- laGP2km(mla2)


omEGObis <- GParetoptim(model = list(laaskm1, laaskm2), fn = fname, crit = "EHI", nsteps = nsteps,
                        lower = lower, upper = upper, critcontrol = critcontrol,
                        optimcontrol = optimcontrol, reinterpolation = FALSE)

prelaGPsep <- function(X, Z, start, end, phi = 1, method = "nn", d = rep(0.1, ncol(X)), g = NULL, ...){
  if(is.null(d)) d <- darg(NULL, X)
  if(!is.list(d)) d <- darg(d, X)
  print(d)
  if(is.null(g)) g <- garg(list(mle = TRUE), Z)
  if(length(g) == 1) g <- list(mle = FALSE, start = g, max = g, min = sqrt(.Machine$double.eps), ab = c(0, 0))
  
  phi <- 1 # Global variance is not relevant here
  res <- list(X = X, Z = Z, da = d, ga = g, phi = phi, 
              start = as.integer(start), end = as.integer(end), method = method)
  class(res) <- 'prelaGPsep'
  return(res)
}

mlasep1 <- prelaGPsep(design.grid, response.grid[,1], start = nappr-5, end = nappr-1, g = 1e-7)
mlasep2 <- prelaGPsep(design.grid, response.grid[,2], start = nappr-5, end = nappr-1, g = 1e-7)


lasepaskm1 <- laGPsep2km(mlasep1)
lasepaskm2 <- laGPsep2km(mlasep2)

omEGOter <- GParetoptim(model = list(lasepaskm1, lasepaskm2), fn = fname, crit = "EHI", nsteps = nsteps,
                        lower = lower, upper = upper, critcontrol = critcontrol,
                        optimcontrol = optimcontrol, reinterpolation = FALSE)


plotGPareto(omEGObis, control = list(PF.line.col = "green", col = "gray"))
plotGPareto(omEGO1, add = TRUE)
plotGPareto(omEGOter, add = TRUE, control = list(PF.line.col = "yellow", col = "orange"))
points(test.vals, pch = '.')

# ## Not yet working for GPGame
# library(GPGame)
# 
# ## Not run: 
# # To use parallel computation (turn off on Windows)
# library(parallel)
# parallel <- FALSE # TRUE # 
# if(parallel) ncores <- detectCores() else ncores <- 1
# 
# ##############################################
# # x.to.obj indicates that P1 chooses x1 and P2 chooses x2
# x.to.obj   <- c(1,2)
# 
# ##############################################
# # Define a discretization of the problem: each player can choose between 21 strategies
# # The ensemble of combined strategies is a 21x21 cartesian grid
# 
# # n.s is the number of strategies (vector)
# n.s <- rep(21, 2)
# # gridtype is the type of discretization
# gridtype <- 'cartesian'
# 
# integcontrol <- list(n.s=n.s, gridtype=gridtype)
# 
# ##############################################
# # Run solver with 6 initial points, 14 iterations
# n.init <- 6 # number of initial points (space-filling)
# n.ite <- 10 # number of iterations (sequential infill points)
# 
# 
# res <- solve_game(fname, model = model,
#                   equilibrium = "NE", crit = "sur", n.init=n.init, n.ite=n.ite,
#                   d = 2, nobj=2, x.to.obj = x.to.obj, integcontrol=integcontrol,
#                   ncores = ncores, trace=1, seed=1)
# 
# resbis <- solve_game(fname, model = list(laaskm1,laaskm2),
#                      equilibrium = "NE", crit = "sur", n.init = n.init, n.ite = n.ite,
#                      d = 2, nobj = 2, x.to.obj = x.to.obj, integcontrol = integcontrol,
#                      ncores = ncores, trace = 1, seed=1)
# 
# 
# par(mfrow = c(1, 3))
# plotGameGrid(P1, n.grid = 21, graphs = "objective")
# plotGame(res, add = TRUE)
# plotGameGrid(P1, n.grid = 21, graphs = "objective")
# plotGame(resbis, add = TRUE)
# plotGameGrid(P1, n.grid = 21, graphs = "objective")
# plotGame(rester, add = TRUE, )
# par(mfrow = c(1, 1))




