# b)

beale_rho <- function(theta) {
  theta1 <- theta[1]
  theta2 <- theta[2]
  return((1.5 - theta1 + theta1*theta2)^2 +
         (2.25 - theta1 + theta1*theta2^2)^2 +
         (2.625 - theta1 + theta1*theta2^3)^2)
}

beale_gradient <- function(theta) {
  theta1 <- theta[1]
  theta2 <- theta[2]
  grad1 <- 2*theta1*theta2^6 + 2*theta1*theta2^4 - 4*theta1*theta2^3 -
    2*theta1*theta2^2 - 4*theta1*theta2 + 6*theta1 + 3*theta2 + 4.5*theta2^2 +5.25*theta2^3 - 12.75
  grad2 <- 6*theta1^2*theta2^5 - 6*theta1^2*theta2^2 + 4*theta1^2*theta2^3 -
    2*theta1^2*theta2 + 9*theta1*theta2 - 2*theta1^2 + 3*theta1 + 15.75
  return(c(grad1, grad2))
}

# c)

###############################################################################

# Grid Line Search
gridLineSearch <- function(theta, rhoFn, d, lambdaStepsize = 0.01, lambdaMax = 1) {
  ## grid of lambda values to search
  lambdas <- seq(from = 0, by = lambdaStepsize, to = lambdaMax)
  ## line search
  rhoVals <- sapply(lambdas, function(lambda) {
    rhoFn(theta - lambda * d)
  })
  ## Return the lambda that gave the minimum
  lambdas[which.min(rhoVals)]
}

# Test Convergence
testConvergence <- function(thetaNew, thetaOld, tolerance = 1e-10, relative = FALSE) {
  sum(abs(thetaNew - thetaOld)) < if (relative)
    tolerance * sum(abs(thetaOld)) else tolerance
}

# Gradient Descent
gradientDescent <- function(theta = 0, rhoFn, gradientFn, lineSearchFn, testConvergenceFn,
                            maxIterations = 100, tolerance = 1e-06, relative = FALSE, lambdaStepsize = 0.01,
                            lambdaMax = 0.5) {
  
  converged <- FALSE
  i <- 0
  
  while (!converged & i <= maxIterations) {
    g <- gradientFn(theta)  ## gradient
    glength <- sqrt(sum(g^2))  ## gradient direction
    if (glength > 0)
      d <- g/glength
    
    lambda <- lineSearchFn(theta, rhoFn, d, lambdaStepsize = lambdaStepsize,
                           lambdaMax = lambdaMax)
    
    thetaNew <- theta - lambda * d
    converged <- testConvergenceFn(thetaNew, theta, tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew
    i <- i + 1
  }
  
  ## Return last value and whether converged or not
  list(theta = theta, converged = converged, iteration = i, fnValue = rhoFn(theta))
}

###############################################################################

# i.
gd_res <- gradientDescent(theta = c(3, 3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)
print(gd_res)

# ii.
gd_res <- gradientDescent(theta = c(3, -3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)
print(gd_res)

# iii.
gd_res <- gradientDescent(theta = c(-3, -3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)
print(gd_res)

# iv.
gd_res <- gradientDescent(theta = c(-3, 3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)
print(gd_res)

# d)

# Modified Gradient Descent
gradientDescent <- function(theta = 0, rhoFn, gradientFn, lineSearchFn, testConvergenceFn,
                            maxIterations = 100, tolerance = 1e-06, relative = FALSE, lambdaStepsize = 0.01,
                            lambdaMax = 0.5) {
  converged <- FALSE
  i <- 0
  xpath <- c(theta[1])
  ypath <- c(theta[2])
  
  while (!converged & i <= maxIterations) {
    g <- gradientFn(theta)  ## gradient
    glength <- sqrt(sum(g^2))  ## gradient direction
    if (glength > 0)
      d <- g/glength
    
    lambda <- lineSearchFn(theta, rhoFn, d, lambdaStepsize = lambdaStepsize,
                           lambdaMax = lambdaMax)
    
    thetaNew <- theta - lambda * d
    converged <- testConvergenceFn(thetaNew, theta, tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew
    i <- i + 1
    xpath[i + 1] <- theta[1]
    ypath[i + 1] <- theta[2]
  }
  
  ## Return point of convergence
  list(xpath = xpath, ypath = ypath)
}

###############################################################################

n_pts <- 100 # number of grid points per dimension
t1_surf <- seq(from = -5, to = 5, length.out = n_pts)
t2_surf <- seq(from = -5, to = 5, length.out = n_pts)

cont_mat <- matrix(0, nrow = n_pts, ncol = n_pts)

for (i in 1:n_pts) {
  for (j in 1:n_pts) {
    cont_mat[i, j] <- beale_rho(c(t1_surf[i], t2_surf[j]))
  }
}

levels <- c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5)



# i.
path1 <- gradientDescent(theta = c(3, 3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

# ii.
path2 <- gradientDescent(theta = c(3, -3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

# iii.
path3 <- gradientDescent(theta = c(-3, -3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

# iv.
path4 <- gradientDescent(theta = c(-3, 3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

contour(
  x = t1_surf,
  y = t2_surf,
  z = cont_mat,
  levels = levels,
  col = "darkgrey",
  main = "Contour Plot of Beale Function",
  xlab = "theta1",
  ylab = "theta2",
  cex = 2,
  cex.lab = 2,
  cex.axis = 2,
  cex.main = 2
)

make_segments <- function(path, colour) {
  i <- 1
  points(x = path$xpath[1], y = path$ypath[1], col = colour, pch = 19)
  while (i != length(path$xpath)) {
    x0 <- path$xpath[i]
    y0 <- path$ypath[i]
    x1 <- path$xpath[i + 1]
    y1 <- path$ypath[i + 1]
    segments(x0 = x0, y0 = y0, x1 = x1, y1 = y1, col = colour, lwd = 1, lty = 1)
    i <- i + 1
  }
  points(x = path$xpath[i], y = path$ypath[i], col = colour, pch = 23)
}
############################################## ADD LABELS FOR Z VALUE OF ENDPOINTS
make_segments(path1, "red")
make_segments(path2, "blue")
make_segments(path3, "yellow")
make_segments(path4, "green")
global_min <- c(3, 0.5)
points(x = global_min[1], y = global_min[2], col = "black", pch = 23)
legend("bottomright",
       legend = c("Gradient descent for (3,3)",
                  "Gradient descent for (3,-3)",
                  "Gradient descent for (-3,-3)",
                  "Gradient descent for (-3,3)",
                  "Global maximum"),
       col = c("red", "blue", "yellow", "green", "black"),
       pch = c(NA, NA, NA, NA, 23), lty = c(1, 1, 1, 1, NA), lwd = 1, cex = 1.75)

# e)
# As shown in part c) and d), it is important to choose a starting value close to
# the global minimum.

# f)
beale_psi <- beale_gradient
beale_psi_prime <- function(theta) {
  theta1 <- theta[1]
  theta2 <- theta[2]
  val <- matrix(0, nrow = length(theta), ncol = length(theta))
  val[1, 1] <- 2*theta2^2 + 2*theta2^4 - 4*theta2^3 - 2*theta2^2 - 4*theta2 + 6
  val[1, 2] <- 12*theta1*theta2^5 + 8*theta1*theta2^3 - 12*theta1*theta2^2 -
               4*theta1*theta2 - 4*theta1 + 3 + 9*theta2 + 15.75*theta2^2
  val[2, 1] <- 12*theta1*theta2^5 - 12*theta1*theta2^2 + 8*theta1*theta2^3 -
               4*theta1*theta2 + 9*theta2 - 4*theta1 + 3
  val[2, 2] <- 30*theta1^2*theta2^4 - 12*theta1^2*theta2 +
               12*theta1^2*theta2^2 - 2*theta1^2 + 9*theta1
  return(val)
}

# g)

###############################################################################

# Newton Raphson
NewtonRaphson <- function(theta, psiFn, psiPrimeFn, dim, testConvergenceFn = testConvergence,
                          maxIterations = 100, tolerance = 1e-06, relative = FALSE) {
  if (missing(theta)) {
    ## need to figure out the dimensionality
    if (missing(dim)) {
      dim <- length(psiFn())
    }
    theta <- rep(0, dim)
  }
  converged <- FALSE
  i <- 0
  while (!converged & i <= maxIterations) {
    thetaNew <- theta - solve(psiPrimeFn(theta), psiFn(theta))
    converged <- testConvergenceFn(thetaNew, theta, tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew
    i <- i + 1
  }
  ## Return last value and whether converged or not
  list(theta = theta, converged = converged, iteration = i, fnValue = psiFn(theta))
}

###############################################################################

# i.

nr_res1 <- NewtonRaphson(theta = c(3, 3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime, testConvergenceFn = testConvergence)
print(nr_res1)

# ii.
nr_res2 <- NewtonRaphson(theta = c(3, -3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime, testConvergenceFn = testConvergence)
print(nr_res2)

# iii.
nr_res3 <- NewtonRaphson(theta = c(-3, -3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime, testConvergenceFn = testConvergence)
print(nr_res3)

# iv.
nr_res4 <- NewtonRaphson(theta = c(-3, 3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime, testConvergenceFn = testConvergence)
print(nr_res4)

# h)

###############################################################################

# Modified Newton Raphson
NewtonRaphson <- function(theta, psiFn, psiPrimeFn, dim, testConvergenceFn = testConvergence,
                          maxIterations = 100, tolerance = 1e-06, relative = FALSE) {
  if (missing(theta)) {
    ## need to figure out the dimensionality
    if (missing(dim)) {
      dim <- length(psiFn())
    }
    theta <- rep(0, dim)
  }
  converged <- FALSE
  i <- 0
  xpath <- c(theta[1])
  ypath <- c(theta[2])
  
  while (!converged & i <= maxIterations) {
    thetaNew <- theta - solve(psiPrimeFn(theta), psiFn(theta))
    converged <- testConvergenceFn(thetaNew, theta, tolerance = tolerance,
                                   relative = relative)
    
    theta <- thetaNew
    i <- i + 1
    xpath[i + 1] <- theta[1]
    ypath[i + 1] <- theta[2]
  }
  ## Return paths
  list(xpath = xpath, ypath = ypath)
}

###############################################################################

n_pts <- 100 # number of grid points per dimension
t1_surf <- seq(from = -5, to = 5, length.out = n_pts)
t2_surf <- seq(from = -5, to = 5, length.out = n_pts)

cont_mat <- matrix(0, nrow = n_pts, ncol = n_pts)

for (i in 1:n_pts) {
  for (j in 1:n_pts) {
    cont_mat[i, j] <- beale_rho(c(t1_surf[i], t2_surf[j]))
  }
}

levels <- c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5)

# i.
path1 <- NewtonRaphson(theta = c(3, 3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime)

# ii.
path2 <- NewtonRaphson(theta = c(3, -3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime)

# iii.
path3 <- NewtonRaphson(theta = c(-3, -3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime)

# iv.
path4 <- NewtonRaphson(theta = c(-3, 3), psiFn = beale_psi, psiPrimeFn = beale_psi_prime)

contour(
  x = t1_surf,
  y = t2_surf,
  z = cont_mat,
  levels = levels,
  col = "darkgrey",
  main = "Contour Plot of Beale Function",
  xlab = "theta1",
  ylab = "theta2",
  cex = 2,
  cex.lab = 2,
  cex.axis = 2,
  cex.main = 2
)

make_segments <- function(path, colour) {
  i <- 1
  points(x = path$xpath[1], y = path$ypath[1], col = colour, pch = 19)
  while (i != length(path$xpath)) {
    x0 <- path$xpath[i]
    y0 <- path$ypath[i]
    x1 <- path$xpath[i + 1]
    y1 <- path$ypath[i + 1]
    segments(x0 = x0, y0 = y0, x1 = x1, y1 = y1, col = colour, lwd = 1, lty = 1)
    i <- i + 1
  }
  points(x = path$xpath[i], y = path$ypath[i], col = colour, pch = 23)
}
############################################## ADD LABELS FOR Z VALUE OF ENDPOINTS
make_segments(path1, "red")
make_segments(path2, "blue")
make_segments(path3, "yellow")
make_segments(path4, "green")
global_min <- c(3, 0.5)
points(x = global_min[1], y = global_min[2], col = "black", pch = 23)
legend("bottomright",
       legend = c("Gradient descent for (3,3)",
                  "Gradient descent for (3,-3)",
                  "Gradient descent for (-3,-3)",
                  "Gradient descent for (-3,3)",
                  "Global maximum"),
       col = c("red", "blue", "yellow", "green", "black"),
       pch = c(NA, NA, NA, NA, 23), lty = c(1, 1, 1, 1, NA), lwd = 1, cex = 1.75)

# i)


# newtonr does well if init val is close to min
# newtonr takes into account 2nd derivative so it checks for curvature
# can go up hill

































