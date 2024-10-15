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
  sample_x <- sample(x, n)
  sample_y <- y[which(x == sample_x)]
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

waldo <- read.csv(file = waldo.csv, header = TRUE)

