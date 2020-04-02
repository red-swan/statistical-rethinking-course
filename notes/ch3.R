


globe_posterior <- function(n, w, d){
  p_grid <- seq( from=0 , to=1 , length.out=n )
  prior <- rep( 1 , n )
  likelihood <- dbinom( w , size=d , prob=p_grid )
  posterior <- likelihood * prior
  list(posterior = posterior / sum(posterior),
       pGrid = p_grid)
}


w3Post <- globe_posterior(1000,3,3)
samples <- sample(x = w3Post$pGrid, prob = w3Post$posterior, size = 1e4, replace = TRUE)

# Loss functions compare point estimates
# Different loss functions create different loss functions
map_dbl(w3Post$pGrid, ~ sum(w3Post$posterior * abs(.x - w3Post$pGrid)))
plot(w3Post$pGrid, .Last.value)

# sampling from a posterior
# grid approximation
hist(samples)




library(rethinking)
dens(samples)

HPDI(samples,0.5)

chainmode(samples, adj = 0.01)


# Posterior Predictive Distribution ============================================
# run the globe tossing experiment for 60% water and see what would
# happen if we ran it 10,000 times
w <- rbinom(1e4, 9, .6)
# What does our model state would happen if we ran the experiment 10,000 times
# remember this is based on n = 3 and catching water each time
w <- rbinom(1e4, 9, samples)




# Exercises ====================================================================
# given
p_grid <- seq( from=0 , to=1 , length.out=1001 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(100)
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

# easy
mean(samples < 0.2)
mean(0.8 < samples)
mean(0.2 < samples & samples < 0.8)
quantile(samples,0.2)
quantile(samples,0.8)
HPDI(samples,0.66)
quantile(samples,c(1/6,5/6))

# medium =========================================
standardize <- function(x) x / sum(x)
likelihood <- dbinom(8,15, p_grid)
posterior <- 
  likelihood %>% 
  multiply_by(prior) %>% 
  standardize()  
samples <- sample(p_grid,1e4,TRUE,posterior)
HPDI(samples, .9)
rbinom(1e4,15,samples) %>% table %>% standardize 
# 6 in 9 waters
rbinom(1e4,9,samples) %>% table %>% standardize
prior <- ifelse(0.5 < p_grid, 1, 0)
# actual values
dbinom(0:9,9,0.7) %>% set_names(0:9)

# hard ===========================================
data(homeworkch3)
prior <- map_dbl(p_grid, ~ 1)
likelihood <- dbinom(sum(birth1) + sum(birth2), 200, p_grid)
posterior <- standardize(likelihood * prior)
p_grid[which.max(posterior)]
samples <- sample(p_grid,1e4,TRUE,posterior)
HPDI(samples,prob = c(0.5,0.89,0.97))
rbinom(1e4,200, samples) %>% dens
  abline(v = 111, col = 'red')

rbinom(1e4,100, samples) %>% dens
  abline(v = 51, col = 'red')
  
rbinom(1e4,sum(birth1 == 0), samples) %>% dens
  abline(v = sum(birth2[birth1 == 0]), col = 'red')
  