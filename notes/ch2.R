







# 2.4.1 Grid Approximation

N <- 1000

# define grid
p_grid <- seq( from=0 , to=1 , length.out=N )
# define prior
prior <- rep( 1 , N )
prior <- ifelse( p_grid < 0.5 , 0 , 1 )
prior <- exp( -5*abs( p_grid - 0.5 ) )
# compute likelihood at each value in grid
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( str_c(N," points" ))

# 2.4.2 Quadratic Approximation

library(rethinking)

globe.qa <- map(
  alist(
    w ~ dbinom(9,p) , # binomial likelihood
    p ~ dunif(0,1) # uniform prior
  ) ,
  data=list(w=6) )
# display summary of quadratic approximation
precis( globe.qa )

# analytical calculation
w <- 6
n <- 9
curve( dbeta( x , w+1 , n-w+1 ) , from=0 , to=1 )
# quadratic approximation
curve( dnorm( x , 0.67 , 0.16 ) , lty=2 , add=TRUE )

# 2.4.3 Markov chain Monte Carlo





