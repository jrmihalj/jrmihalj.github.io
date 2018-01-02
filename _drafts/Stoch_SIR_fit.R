# Load some required packages
##############
library(deSolve)
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
##############

# Create a time series over which to integrate.
# Here we have an epidemic that is observed over t_max number of days (or weeks or etc).
t_min = 1
t_max = 50
times <- seq(t_min, t_max, length = (t_max-t_min)*10)

# To simulate the data, we need to assign initial conditions.
# In practice, these will likely be unknown, but can be estimated from the data.
N0 = 100
I0 = 1      # initial fraction infected
S0 = N0 - I0  # initial fraction susceptible
R0 = 0        # initial fraction removed

# Assign transmission and pathogen-induced death rates:
beta_mean = 0.005
gamma = 0.1

# Stochastic variation in beta (daily transmission rate):
stoch_var = 0.35
epsilon = exp(rnorm(t_max-t_min, 0, stoch_var))
beta_adj = rep(epsilon * beta_mean, each=10)
#unique(beta_adj)
#hist(beta_adj)

# We will use the package deSolve to integrate, which requires certain data structures.
# Store parameters and initial values
# Parameters must be stored in a named list.
params <- list(beta_adj = beta_adj,
               #beta = beta_mean,
               gamma = gamma)

# Initial conditions are stored in a vector
inits <- c(S=S0, I=I0, R=R0)

# We must create a function for the system of ODEs.
# See the 'ode' function documentation for further insights.
# SIR <- function(t, y, params) {
#   with(as.list(c(params, y)), {
#     
#     for(i in 1:3){
#       if(is.na(y[i]) | is.infinite(y[i])){
#         y[i] <- 0
#       }
#       if(y[i] <= 0){
#         y[i] <- 0
#       }
#     }
#     
#     beta <- beta_adj[t]
#     
#     dS <- - beta * y[1] * y[2] 
#     dI <- beta * y[1] * y[2] - gamma * y[2]
#     dR <- gamma * y[2]
#     
#     res <- c(dS,dI,dR)
#     list(res)
#   })
# }

# Run the integration:
out <- ode(y = inits, 
           times = times, 
           func = SIR, 
           parms = params, 
           method="ode45")

# Store the output in a data frame:
out <- data.frame(out)

# Quick plot of the epidemic
plot(NA,NA, xlim = c(t_min, t_max), ylim=c(0, N0),
     xlab = "Time", ylab="Density of Host Population")
lines(out$S ~ out$time, col="black")
lines(out$I ~ out$time, col="red")
legend(x = 30, y = 0.8, legend = c("Susceptible", "Infected"),
       col = c("black", "red"), lty = c(1, 1), bty="n")

# Fraction of hosts infected
plot(NA,NA, xlim = c(t_min, t_max), ylim=c(0, 1), 
     xlab = "Time", ylab="Fraction Infected")
lines(out$I / (out$I + out$S) ~ out$time, col="red")

####################################
####################################

sample_days = 20 # number of days sampled throughout the epidemic
sample_n = 25 # number of host individuals sampled per day

# Choose which days the samples were taken. 
# Ideally this would be daily, but we all know that is difficult.
sample_time = sort(sample(1:t_max, sample_days, replace=F))

# Extract the "true" fraction of the population that is infected on each of the sampled days:
sample_I = out[round(out$time, 1) %in% sample_time, 3]
sample_S = out[round(out$time, 1) %in% sample_time, 2]
sample_propinf = sample_I / (sample_I + sample_S)

# Generate binomially distributed data.
# So, on each day we sample a given number of people (sample_n), and measure how many are infected.
# We expect binomially distributed error in this estimate, hence the random number generation.
sample_y = rbinom(sample_days, sample_n, sample_propinf)

# Plot samples:
plot(sample_y/sample_n ~ sample_time, ylim=c(0,1), xlim=c(t_min, t_max))
