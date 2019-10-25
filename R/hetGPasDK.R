## Try using hetGP with GPareto and GPGame
## Note: for covMat1Mat2 to work, the definition code somehow needs to be ran twice
library(hetGP)

setOldClass("homGP")

# matching the covariance structure
setClass("covHetGP2covkm", 
         slots = list(sd2 = "numeric", covtype = "character", theta = "numeric", nugget.flag = "logical"))


covHetGP2covkm <- function(model){
  return(new(Class = "covHetGP2covkm", sd2 = model$nu_hat, covtype = model$covtype,
             theta = model$theta, nugget.flag = FALSE))
}


if(!isGeneric("covMat1Mat2")) {
  setGeneric(name    = "covMat1Mat2",
             def     = function(object, X1, X2, nugget.flag=FALSE) standardGeneric("covMat1Mat2")
  )
}

covMat1Mat2 <- function (object, X1, X2, nugget.flag) {
  UseMethod("covMat1Mat2", object, X1, X2, nugget.flag)
}

covMat1Mat2.covHetGP2covkm <- function(object, X1, X2, nugget.flag=FALSE){
  if(!nugget.flag) cov_gen(X1 = X1, X2 = X2, theta = object@theta, type = object@covtype)
}

setMethod("covMat1Mat2", "covHetGP2covkm", 
          definition = function(object, X1, X2, nugget.flag = FALSE){
            covMat1Mat2.covHetGP2covkm(object, X1, X2, nugget.flag)
          }
)


setClass("homGP2km", slots = list(model = "homGP", X = "matrix", y = "matrix", d = "integer", n = "integer",
                                  noise.var = "numeric", noise.flag = "logical", 
                                  covariance = "covHetGP2covkm", trend.formula = "formula", trend.coef = "numeric",
                                  M = "matrix", T = "matrix", F = "matrix", z = "numeric"))
# newmod <- new(Class = "homGP2km", model = mhom1, X = mhom1$X0)

homGP2km <- function(model){
  # if(class(model) == "homGP") class(model) <- "list" else stop("Model is not a homGP object")
  C <- cov_gen(X1 = model$X0, theta = model$theta, type = model$covtype)
  T <- chol(C + diag(model$g, nrow(model$X0)))/model$nu_hat
  F <- matrix(1, nrow = nrow(model$X0), 1)
  M <- backsolve(t(T), F, upper.tri = FALSE)
  z <- as.numeric(backsolve(t(T), model$Z0 - F * model$beta0, upper.tri=FALSE))
  
  res <- new(Class = "homGP2km", model = model, X = model$X0, y = matrix(model$Z0, ncol = 1),
             d = ncol(model$X0), n = nrow(model$X0), noise.var = model$g, noise.flag = TRUE,
             covariance = covHetGP2covkm(model), trend.formula = ~1, trend.coef = model$beta0,
             # covariance = new(Class = "covHetGP", sd2 = model$nu_hat, covtype = model$covtype, theta = model$theta, nugget.flag = FALSE)
             F = F, M = M, T = T, z = z
  )
  
  return(res)
}

setMethod("predict", "homGP2km", function(object, newdata, type, se.compute = TRUE, 
                                          cov.compute = FALSE, light.return = FALSE,
                                          bias.correct = FALSE, checkNames = TRUE, ...){
  xprime <- if(cov.compute) xprime <- newdata else xprime <- NULL
  
  if(class(newdata) == "data.frame") newdata <- as.matrix(newdata)
  preds <- predict(object@model, x = newdata, xprime = xprime)
  if(!is.null(xprime)) preds$cov <- 1/2*(preds$cov + t(preds$cov)) # to ensure symmetry
  return(list(mean = preds$mean, sd = sqrt(preds$sd2), cov = preds$cov))
})

setMethod("update", "homGP2km", function(object, newX, newy, newX.alreadyExist = FALSE,
                                         cov.reestim = TRUE, trend.reestim = TRUE, nugget.reestim = FALSE, 
                                         newnoise.var = NULL, kmcontrol = NULL, newF = NULL,...){
  res <- update(object@model, Xnew = newX, Znew = newy)
  return(homGP2km(res))
})


KK <- matrix()
library(MASS)
setMethod("simulate", "homGP2km", function(object, nsim=1, seed=NULL, newdata=NULL,
                                           cond=FALSE, nugget.sim=0, checkNames=TRUE, ...){
  if(cond==FALSE){
    mu <- rep(object$beta0, nrow(newdata))
    K <- cov_gen(newdata, theta = object@covariance@theta, type = object@model$covtype)
  }else{
    preds <- predict(object = object, newdata = newdata, cov.compute = TRUE)
    mu <- preds$mean
    K <- preds$cov
  }
  KK <<- K
  print(K)
  print(eigen(K)$values)
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

# Same with homGP

mhom1 <- mleHomGP(X = design.grid, Z = response.grid[,1])
mhom2 <- mleHomGP(X = design.grid, Z = response.grid[,2])


homaskm1 <- homGP2km(mhom1)
homaskm2 <- homGP2km(mhom2)


omEGObis <- GParetoptim(model = list(homaskm1,homaskm2), fn = fname, crit = "EHI", nsteps = nsteps,
                        lower = lower, upper = upper, critcontrol = critcontrol,
                        optimcontrol = optimcontrol, reinterpolation = FALSE)

plotGPareto(omEGObis, control = list(PF.line.col = "green", col = "gray"))
plotGPareto(omEGO1, add = TRUE)

## With GPGame
library(GPGame)

## Not run: 
# To use parallel computation (turn off on Windows)
library(parallel)
parallel <- FALSE # TRUE # 
if(parallel) ncores <- detectCores() else ncores <- 1

##############################################
# 2 variables, 2 players, Nash equilibrium
# Player 1 (P1) wants to minimize fun1 and player 2 (P2) fun2
# P1 chooses x2 and P2 x2

##############################################
# First, define objective function fun: (x1,x2) -> (fun1,fun2)
fun <- function (x)
{
  if (is.null(dim(x)))    x <- matrix(x, nrow = 1)
  b1 <- 15 * x[, 1] - 5
  b2 <- 15 * x[, 2]
  return(cbind((b2 - 5.1*(b1/(2*pi))^2 + 5/pi*b1 - 6)^2 + 10*((1 - 1/(8*pi)) * cos(b1) + 1),
               -sqrt((10.5 - b1)*(b1 + 5.5)*(b2 + 0.5)) - 1/30*(b2 - 5.1*(b1/(2*pi))^2 - 6)^2-
                 1/3 * ((1 - 1/(8 * pi)) * cos(b1) + 1)))
}

##############################################
# x.to.obj indicates that P1 chooses x1 and P2 chooses x2
x.to.obj   <- c(1,2)

##############################################
# Define a discretization of the problem: each player can choose between 21 strategies
# The ensemble of combined strategies is a 21x21 cartesian grid

# n.s is the number of strategies (vector)
n.s <- rep(21, 2)
# gridtype is the type of discretization
gridtype <- 'cartesian'

integcontrol <- list(n.s=n.s, gridtype=gridtype)

##############################################
# Run solver with 6 initial points, 14 iterations
n.init <- 6 # number of initial points (space-filling)
n.ite <- 10 # number of iterations (sequential infill points)


res <- solve_game(fun, model = model,
                  equilibrium = "NE", crit = "sur", n.init=n.init, n.ite=n.ite,
                  d = 2, nobj=2, x.to.obj = x.to.obj, integcontrol=integcontrol,
                  ncores = ncores, trace=1, seed=1)

resbis <- solve_game(fun, model = list(homaskm1,homaskm2),
                  equilibrium = "NE", crit = "sur", n.init=n.init, n.ite=n.ite,
                  d = 2, nobj=2, x.to.obj = x.to.obj, integcontrol=integcontrol,
                  ncores = ncores, trace=1, seed=1)

par(mfrow = c(1, 2))
plotGameGrid(P1, n.grid = 21, graphs = "objective")
plotGame(res, add = TRUE)
plotGameGrid(P1, n.grid = 21, graphs = "objective")
plotGame(resbis, add = TRUE, )
par(mfrow = c(1, 1))

