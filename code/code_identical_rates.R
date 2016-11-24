# Libraries 
library(deSolve) # integrate ODEs
library(lpSolve) # linear programming
library(ggplot2) # for plotting
library(tidyr) # for organizing data
library(dplyr) # for organizing data
library(far) # for orthogonal vectors

# Global variables
Cutoff <- 25 # Cutoff value, to determine whether a species has 0 density 

# Build a random matrix H such that H_ij + H_ji = 1
random_H <- function(n){
  # build random hypertournament H
  H <- matrix(runif(n * n), n, n)
  return(H / (H + t(H)))
}

# Find the optimal strategy for the two-person game encoded in H
# using linear programming.
# This is also the coexistence equilibrium of the dynamical system.
find_optimal_strategy <- function(H){
  n <- dim(H)[1]
  f.obj <- rep(1, n)
  f.con <- H
  f.rhs <- rep(1, n)
  f.dir <- rep("<=", n)
  z <- lp ("max", f.obj, f.con, f.dir, f.rhs)
  return(z$solution / sum(z$solution))
}

# Functions for integrating the systems of ODEs

# Sample 2 seedlings at a time
eq1 <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  H <- pars$H
  dxidt <- 2 * x * H %*% x - x
  return(list(dxidt))
}

# Sample 3 seedlings at a time
eq3 <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  H <- pars$H
  dxidt <- x * (2 * (H %*% x)^2 + H %*% (2 * x *  H %*% x) - 1)
  return(list(dxidt))
}

# Sample 4 seedlings at a time
sample4 <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  H <- pars$H
  x2 <- 2 * x * H %*% x
  x3 <- x * H %*% x2 + x2 * H %*% x
  x4 <- x * H %*% x3 + x3 * H %*% x
  dxidt <- x4 - x
  return(list(dxidt))
}

# Sample 5 seedlings at a time
sample5 <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  H <- pars$H
  x2 <- 2 * x * H %*% x
  x3 <- x * H %*% x2 + x2 * H %*% x
  x4 <- x * H %*% x3 + x3 * H %*% x
  x5 <- x * H %*% x4 + x4 * H %*% x
  dxidt <- x5 - x
  return(list(dxidt))
}

# Equivalent formulation as replicator equation
# The dynamics match those of eq1
eq2 <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  P2 <- pars$P2
  dxidt <- x *  P2 %*% x
  return(list(dxidt))
}

# Equivalent formulation as replicator equation
# Now payoff is a tensor. Match the dynamics of eq3
eq4 <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  P <- pars$P3
  n <- length(x)
  dxidt <- rep(0, n)
  for (i in 1:n){
    dxidt[i] <- x[i] * x %*% P[i,,] %*% x
  }
  return(list(dxidt))
}

# Master function to run the dynamics
run_dynamics <- function(state = c(1/4, 1/2, 1/4),
                         maxtime = 100, 
                         stepout = 0.1, 
                         model = "eq1", 
                         H = matrix(c(1/2,1,0, 0,1/2,1, 1,0,1/2), 3, 3, byrow = TRUE)
){
  # Build matrix P2 (Payoff matrix, 2-player game)
  P2 <- H - t(H)
  # Build matrix P3 (Payoff tensor, 2-player game)
  n <- nrow(H)
  P3 <- array(0, c(n,n,n))
  for (i in 1:n){
    for (j in 1:n){
      for (k in 1:n){
        P3[i,j,k] <- 2 * H[i,j] * H[i,k] -  H[j,i] * H[j,k] - H[k,i] * H[k,j]
      }
    }
  }
  pars <- list(H = H, 
               P2 = P2, 
               P3 = P3)
  time <- seq(0, maxtime, by = stepout)
  out <- as.data.frame(ode(func = match.fun(model), y = state, parms = pars, times = time, method = "ode45"))
  output <- list(model = model, 
                 time = time, 
                 time_series = out[,-1], 
                 H = H, 
                 equil_approx = apply(round(out[as.integer(3 * nrow(out) / 4):nrow(out) ,-1], 7), 2, mean),
                 equil_lp = find_optimal_strategy(H))
  return(output)
}

# Plot the time-series
plot_dynamics <- function(output){
  time_series <- output$time_series
  time <- output$time
  time_series <- cbind(time, time_series)
  colnames(time_series) <- c("time", paste("sp", colnames(time_series[-1])))
  ts <- time_series %>% gather("species", "density", 2:ncol(time_series))
  pl <- ggplot(ts, aes(x = time, y = density, colour = species)) + geom_line() + theme_bw() + theme(legend.position = "bottom")
}

# Build a matrix H such that x is an equilibrium
# The algorithm is given in the supplement
build_H_given_xstar <- function(x){
  # step 0: normalize x so that this is relative species abundance
  target <- x / sum(x)
  n <- length(x)
  # Odd number of species
  if (n %% 2 == 1){
    # step 1: build a set of n orthonormal vectors
    U1 <- matrix(runif(n * n), n, n)
    # use the target as the first vector (this will be unchanged)
    U1[,1] <- target
    # make orthogonal
    U1 <- orthonormalization(U1)
    # Build complex matrix U of eigenvectors
    U <- matrix(complex(0, 0), n, n)
    # Build complex matrix V of eigenvalues
    V <- matrix(complex(0, 0), n, n)
    # First eigenvector is the target
    U[,1] <- complex(real = U1[,1], imaginary = 0)
    # First eigenvalue is 0
    V[1,1] <- complex(real = 0, imaginary = 0)
    j <- 2
    while(j < n){
      # form eigenvectors
      U[, j] <- complex(real = U1[,j] / sqrt(2), imaginary = U1[, j + 1] / sqrt(2))
      U[, j + 1] <- complex(real = U1[,j] / sqrt(2), imaginary = -U1[, j + 1] / sqrt(2))
      # draw random imaginary eigenvalues
      z <- runif(1)
      V[j, j] <- complex(real = 0, imaginary = z)
      V[j + 1, j + 1] <- complex(real = 0, imaginary = -z)
      j <- j + 2 
    }
    # skew-symmetric matrix
    P <- Re(U %*% V %*% solve(U))
  }
  
  # Even number of species
  if (n %% 2 == 0){
    # step 1: build a set of n orthonormal vectors
    U1 <- matrix(runif(n * n), n, n)
    # use the target as the first vector (this will be unchanged)
    U1[,1] <- target
    # make orthogonal
    U1 <- orthonormalization(U1)
    # Build complex matrix U of eigenvectors
    U <- matrix(complex(0, 0), n, n)
    # Build complex matrix V of eigenvalues
    V <- matrix(complex(0, 0), n, n)
    j <- 1
    while(j < n){
      # form eigenvectors
      U[, j] <- complex(real = U1[,j] / sqrt(2), imaginary = U1[, j + 1] / sqrt(2))
      U[, j + 1] <- complex(real = U1[,j] / sqrt(2), imaginary = -U1[, j + 1] / sqrt(2))
      # draw random imaginary eigenvalues
      if (j < 2){
        z <- 0
      } else{
        z <- runif(1)
      }
      V[j, j] <- complex(real = 0, imaginary = z)
      V[j + 1, j + 1] <- complex(real = 0, imaginary = -z)
      j <- j + 2 
    }
    # skew-symmetric matrix
    P <- Re(U %*% V %*% solve(U))
  }
  H <- (P + 1)/2
  return(H)
}