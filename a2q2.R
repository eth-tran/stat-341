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
    2*theta1*theta2 - 4*theta1*theta2 + 6*theta1 + 3*theta2 + 4.5*theta2^2 - 12.75
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
gradientDescent(theta = c(3, 3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

# ii.
gradientDescent(theta = c(3, -3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

# iii.
gradientDescent(theta = c(-3, -3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

# iv.
gradientDescent(theta = c(-3, 3), rhoFn = beale_rho, gradientFn = beale_gradient,
                lineSearchFn = gridLineSearch, lambdaStepsize = 1e-03,
                testConvergenceFn = testConvergence, maxIterations = 500)

# d)

# Modified Gradient Descent
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

contour(
  x = t1_surf,
  y = t2_surf,
  z = cont_mat,
  levels = levels,
  col = "darkgrey"
)



for (i in 1:n_pts) {
  for (j in 1:n_pts) {
    theta1 <- t1_surf[i]
    theta2 <- t2_surf[j]
    summary <- gradientDescent(theta = c(theta1, theta2), rhoFn = beale_rho,
               gradientFn = beale_gradient, lineSearchFn = gridLineSearch,
               lambdaStepsize = 1e-03, testConvergenceFn = testConvergence,
               maxIterations = 500)
    cont_mat[i, j] <- summary$fnValue
  }
}

contour(
  x = t1_surf,
  y = t2_surf,
  z = cont_mat,
  nlevels = 25
)




# newtonr does well if init val is close to min
# newtonr takes into account 2nd derivative so it checks for curvature
# can go up hill
