

# function declarations ========================================================

standardize <- function(distr) distr / sum(distr)

gridApprox <- function(likelihood, prior){
  if(length(likelihood) != length(prior)) stop("likelihood and prior must have same length")
  standardize(likelihood * prior)
}

# questions ====================================================================

#1 -----------------------------------------------------------------------------
# Suppose the globe tossing data had turned out to be 8 water in 15 tosses.
# Construct the posterior distribution, using grid approximation. 
# Use the same flat prior as before.

xGrid <- seq(0,1,length.out = 1001)
priorInit <- rep(1,length(xGrid))
probData <- dbinom(8,15,prob = xGrid)
posterior1 <- gridApprox(probData, priorInit)
plot(xGrid,posterior1)

# 2 ----------------------------------------------------------------------------
# Start over in 1, but now use a prior that is zero below p = 0:5 and a constant
# above p = 0:5. This corresponds to prior information that a majority
# of the Earth's surface is water. What difference does the better prior make?
# If it helps, compare posterior distributions (using both priors) to the true
# value p = 0:7.

prior2 <- ifelse(xGrid >= 0.5, 1, 0)
posterior2 <- gridApprox(probData, prior2)

# comparing plots
plot(posterior2, lty = 1)
lines(posterior1, col = 'red', lty = 2)
abline(v = 701, col = 'green', lty = 'dotted')


# 3 ----------------------------------------------------------------------------
# This problem is more open-ended than the others. Feel free to collaborate
# on the solution. Suppose you want to estimate the Earth's proportion of
# water very precisely. Specifically, you want the 99% percentile interval of the
# posterior distribution of p to be only 0.05 wide. This means the distance between
# the upper and lower bound of the interval should be 0.05. How many
# times will you have to toss the globe to do this? I won't require a precise
# answer. I'm honestly more interested in your approach.

f <- function(n){
  dbinom(n*0.7, n, xGrid) %>% 
  extract(., 0 < .) %>% 
  length() %>%
  subtract(50) %>% 
  raise_to_power(2)
}

minN <- optim(par = 10000, fn = f,method = 'Brent', 
              lower = 400000, upper = 500000)$par


g <- function(n){
  dbinom(n*0.7, n, seq(0.6,0.8,by = 0.01)) %>% 
    extract(., 0 < .) %>% 
    length()
}

h <- function(n) dbinom(n*0.7, n, xGrid)




################################################################################
set.seed(100)
N <- 2000
p_true <- 0.7
W <- rbinom( 1 , size=N , prob=p_true )
p_grid <- seq( from=0 , to=1 , by=0.001 )
prior <- rep( 1 , length(p_grid) )
prob_data <- dbinom( W , size=N , prob=p_grid )
posterior <- prob_data * prior
posterior <- posterior / sum(posterior)
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )
PI99 <- PI( samples , 0.99 )
as.numeric( PI99[2] - PI99[1] )





f <- function( N ) {
  p_true <- 0.7
  W <- rbinom( 1 , size=N , prob=p_true )
  p_grid <- seq( from=0 , to=1 , length.out=1000 )
  prior <- rep( 1 , 1000 )
  prob_data <- dbinom( W , size=N , prob=p_grid )
  posterior <- prob_data * prior
  posterior <- posterior / sum(posterior)
  samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )
  PI99 <- PI( samples , 0.99 )
  as.numeric( PI99[2] - PI99[1] )
}


























