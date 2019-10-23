## Try using hetGP with GPareto
library(hetGP)

setOldClass("homGP")

setClass("covHetGP", slots = list(sd2 = "numeric"))
setClass("homGP2km", slots = list(model = "homGP", X = "matrix", y = "matrix", d = "integer", n = "integer",
                                   noise.var = "numeric", noise.flag = "logical", covariance = "covHetGP"))
# newmod <- new(Class = "homGP2km", model = mhom1, X = mhom1$X0)

homGP2km <- function(model){
  # if(class(model) == "homGP") class(model) <- "list" else stop("Model is not a homGP object")
  res <- new(Class = "homGP2km", model = model, X = model$X0, y = matrix(model$Z0, ncol = 1),
             d = ncol(model$X0), n = nrow(model$X0), noise.var = model$g, noise.flag = TRUE,
             covariance = new(Class = "covHetGP", sd2 = model$nu_hat))

  return(res)
}

setMethod("predict", "homGP2km", function(object, newdata, type, se.compute = TRUE, 
                              cov.compute = FALSE, light.return = FALSE,
                              bias.correct = FALSE, checkNames = TRUE, ...){
  xprime <- if(cov.compute) xprime <- newdata else xprime <- NULL
  
  if(class(newdata) == "data.frame") newdata <- as.matrix(newdata)
  preds <- predict(object@model, x = newdata, xprime = xprime)
  return(list(mean = preds$mean, sd = sqrt(preds$sd2), cov = preds$cov))
})

setMethod("update", "homGP2km", function(object, newX, newy, newX.alreadyExist = FALSE,
                                          cov.reestim = TRUE, trend.reestim = TRUE, nugget.reestim = FALSE, 
                                          newnoise.var = NULL, kmcontrol = NULL, newF = NULL,...){
  res <- update(object@model, Xnew = newX, Znew = newy)
  return(homGP2km(res))
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

nsteps <- 20
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






