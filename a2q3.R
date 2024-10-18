# a)

# i.

dummy_objective <- function(theta) {
  return(NULL)
}

# ii.

dummy_convergence_test <- function(thetaNew, thetaOld, tolerance = 1e-10, relative = FALSE) {
  return(FALSE)
}

# iii.

decay_learning_rate <- function(theta, rhoFn, d, lambdaStepsize = 0.01, lambdaMax = 1, iteration) {
  return(0.1/(1+0.1*iteration))
}

# iv.

create_least_squares_stochastic_gradient <- function(x, y, n) {
  indices <- seq(from = 1, to = length(x), by = 1)
  indices <- sample(x = indices, n)
  sample_x <- x[indices]
  sample_y <- y[indices]
  sample_xbar <- mean(sample_x)
  function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    -2 * c(sum(sample_y - alpha - beta * (sample_x - sample_xbar)),
           sum((sample_y - alpha - beta * (sample_x -
           sample_xbar)) * (sample_x - sample_xbar)))
  }
}

# b)

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
                           lambdaMax = lambdaMax, iteration = i) ########## stepsize???
    
    thetaNew <- theta - lambda * d
    converged <- testConvergenceFn(thetaNew, theta, tolerance = tolerance,
                                   relative = relative)
    theta <- thetaNew
    i <- i + 1
    xpath[i + 1] <- theta[1]
    ypath[i + 1] <- theta[2]
  }
  
  ## Return path
  list(xpath = xpath, ypath = ypath)
}

###############################################################################

waldo <- read.csv(file = "waldo.csv", header = TRUE)
fn <- create_least_squares_stochastic_gradient(x = waldo$X, y = waldo$Y, n = 20)

path <- gradientDescent(
  theta = c(1, 1),
  rhoFn = dummy_objective,
  gradientFn = fn,
  lineSearchFn = decay_learning_rate,
  testConvergenceFn = dummy_convergence_test,
  maxIterations = 5000,
  tolerance = 1e-06,
  relative = FALSE,
  lambdaStepsize = 0.01,
  lambdaMax = 0.5
)

alpha <- path$xpath
beta <- path$ypath

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

path <- gradientDescent(
  theta = c(1, 1),
  rhoFn = dummy_objective,
  gradientFn = fn,
  lineSearchFn = decay_learning_rate,
  testConvergenceFn = dummy_convergence_test,
  maxIterations = 5000,
  tolerance = 1e-06,
  relative = FALSE,
  lambdaStepsize = 0.01,
  lambdaMax = 0.5
)

par(mfrow = c(2,1))
plot(x = seq(from = 1, to = length(alpha)), y = alpha, type = "l", col = "blue", lwd = 2)
plot(x = seq(from = 1, to = length(alpha)), y = beta, type = "l", col = "red", lwd = 2)
abline(h = )

























