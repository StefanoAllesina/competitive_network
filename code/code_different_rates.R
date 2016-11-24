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

# Sample 2 seedlings at a time --- include different physiological rates
sample2_full <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  d <- pars$d
  f <- pars$f
  H <- pars$H
  dxidt <- x * (2 * as.numeric(d %*% x ) * (H %*% (f*x) ) * f / as.numeric(f %*% x )^2 - d)
  return(list(dxidt))
}

# Sample 3 seedlings at a time --- include different physiological rates
sample3_full <- function(time, state, pars){
  state <- round(state, Cutoff)
  state[state < 0] <- 0
  x <- state  / sum(state) # This is needed for numerical precision
  d <- pars$d
  f <- pars$f
  H <- pars$H
  y <- f * x
  Fh <- sum(f * x)
  Dh <- sum(d * x)
  n <- nrow(H)
  Y <- matrix(0, n,n)
  diag(Y) <- y
  dxidt <- x * (2 * Dh * ( H %*% Y %*% H %*% y + (H %*% y)^2 ) * f / Fh^3 - d)
  return(list(dxidt))
}

# Master function to run the dynamics
run_dynamics_physrates <- function(state = c(1/4, 1/2, 1/4),
                         maxtime = 100, 
                         stepout = 0.1, 
                         model = "eq1", 
                         H = matrix(c(1/2,1,0, 0,1/2,1, 1,0,1/2), 3, 3, byrow = TRUE),
                         d = c(1,1,1),
                         f = c(1,1,1.01)
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
               P3 = P3,
               d = d,
               f = f)
  time <- seq(0, maxtime, by = stepout)
  out <- as.data.frame(ode(func = match.fun(model), y = state, parms = pars, times = time, method = "ode45"))
  output <- list(model = model, 
                 time = time, 
                 time_series = out[,-1], 
                 H = H, 
                 d = d,
                 f = f,
                 equil_approx = apply(round(out[as.integer(3 * nrow(out) / 4):nrow(out) ,-1], 7), 2, mean)
                 )
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
