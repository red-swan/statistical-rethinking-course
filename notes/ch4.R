library(swanR)
library(stringr)
library(purrr)
library(magrittr)
library(data.table)
library(lubridate)
library(ggplot2)
library(openxlsx)
library(xts)
library(zoo)
library(microbenchmark)
library(rethinking)
data("Howell1")
d <- Howell1 %>% as.data.table()
d2 <- d[age >= 18]
d2[,dens(height,norm.comp = TRUE)]


curve( dnorm( x , 178 , 20 ) , from=100 , to=250 )


N <- 1e6
sample_mu <- rnorm( N , 178 , 20 )
sample_sigma <- runif( N , 0 , 50 )
prior_h <- rnorm( N , sample_mu , sample_sigma )
dens( prior_h )

# building the grid approximation of the posterior
mu.list <- seq( from=153.5, to=155.8 , length.out=1e3 )
sigma.list <- seq( from=7 , to=9 , length.out=1e3 )
post <- expand.grid( mu=mu.list , sigma=sigma.list ) %>% as.data.table
post[ , LL := dnorm(d2$height, .BY$mu, .BY$sigma, log = TRUE) %>% sum,
      by = .(mu,sigma)]
post[ , prod := LL + dnorm(mu, 178, 20, TRUE) + dunif(sigma, 0, 50, TRUE)]
post[, prob := exp( prod - max(prod) )]

#inspecting the distribution
post[, contour_xyz(mu,sigma, prob)]
post[, image_xyz(mu,sigma,prob)]

# sampling the posterior
sample.rows <- post[,sample(1:.N, 1e4, TRUE, prob)]
samples <- post[ sample.rows , .(mu,sigma)]
samples[,plot(mu,sigma,cex = 0.5,pch = 16, col = col.alpha(rangi2,0.1))]
samples[,dens(mu)]
samples[,dens(sigma)]
samples[, HPDI(mu)]
samples[, HPDI(sigma)]

# running with only a small sample
d3 <- sample( d2$height , size=20 )
mu.list <- seq( from=150, to=170 , length.out=200 )
sigma.list <- seq( from=4 , to=20 , length.out=200 )
post2 <- expand.grid( mu=mu.list , sigma=sigma.list )
post2$LL <- sapply( 1:nrow(post2) , function(i)
  sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] ,
              log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
  dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
                        prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
      col=col.alpha(rangi2,0.1) ,
      xlab="mu" , ylab="sigma" , pch=16 )
dens( sample2.sigma , norm.comp=FALSE )


# quadratic approximation
flist <- alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0, 50)
)
m41 <- map(flist,d2)
precis(m41)
sapply(samples, HPDI)

m42 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 0.1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2 )
precis( m42 )

vcov(m42)
cov2cor(vcov(m42))

# sampling
post <- extract.samples(m41,1e4)
precis(post)


# log variance
m41l <-
  alist(
    height ~ dnorm(mu, exp(log_sigma)),
    mu ~ dnorm(178,20),
    log_sigma ~ dnorm(2,10)
  ) %>% 
  map(d2)
post <- extract.samples(m41l) %>% as.data.table
post[ , sigma := exp(post$log_sigma)]

# 4.4 Adding a predictor ===================================

plot(height ~ weight, d2)

m43 <-
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b * weight,
    a ~ dnorm(156,100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ) %>% 
  map(d2)
precis(m43, corr = TRUE)
cov2cor(vcov(m43))

m43lm <- lm(height ~ weight, d2)
summary(m43lm)

library(brms)
m43brm <- 
  brm(data = d2, family = gaussian,
      formula = height ~ 1 + weight, 
      prior = c(prior(normal(156,100), class = Intercept),
                prior(normal(0,10), class = b),
                prior(normal(0,50),class = sigma)
              ),
      iter = 41000, warmup = 40000, chains = 4, cores = 4, seed = 4 )

summary(m43brm)
plot(m43brm)
# This was another example of how using a uniform prior for σ required we use an
# unusually large number of warmup iterations before the HMC chains converged on
# the posterior. Change the prior to cauchy(0, 1) and the chains converge with 
# no problem, resulting in much better effective samples, too. Here are the 
# trace plots.


# centering
d2[ , weightC := weight - mean(weight)]
m44 <-
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b * weightC,
    a ~ dnorm(156,100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ) %>% 
  map(d2)
precis(m44)
round(cov2cor(vcov(m44)),2)

# plotting
plot( height ~ weight , data=d2 )
abline( a=coef(m43)["a"] , b=coef(m43)["b"] , col = 'red')
post <- extract.samples(m43) %>% as.data.table

# 4.48 only a subset
dataN <- 352
sampleN <- 20
dN <- d2[sample(1:.N,dataN)]
mN <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ dnorm( 178 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=dN )
# extract 20 samples from the posterior
post <- extract.samples( mN , n=sampleN )
# display raw data and sample size
plot( dN$weight , dN$height ,
      xlim=range(d2$weight) , ylim=range(d2$height) ,
      col=rangi2 , xlab="weight" , ylab="height" )
mtext(concat("N = ",dataN))
# plot the lines, with transparency
for ( i in 1:sampleN ){
  abline( a=post$a[i] , b=post$b[i] , col=col.alpha("black",0.3) )
}


muAt50 <- post[,a + b * 50]
dens( muAt50 , col=rangi2 , lwd=2 , xlab="mu|weight=50" )
HPDI(muAt50)
mu <- link(m43)
str(mu)

# define sequence of weights to compute predictions for
# these values will be on the horizontal axis
weight.seq <- seq( from=25 , to=70 , by=1 )
# use link to compute mu
# for each sample from posterior
# and for each weight in weight.seq
mu <- link( m43 , data=data.frame(weight=weight.seq) )
str(mu)

# use type="n" to hide raw data
plot( height ~ weight , d2 , type="n" )
# loop over samples and plot each mu value
for ( i in 1:100 ){
  points( weight.seq , mu[i,] , pch=16 , col=col.alpha(rangi2,0.1) )
}
# summarize the distribution of mu
mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

# plot raw data
# fading out points to make line and interval more visible
plot( height ~ weight , data=d2 , col=col.alpha(rangi2,0.5) )
# plot the MAP line, aka the mean mu for each weight
lines( weight.seq , mu.mean )
# plot a shaded region for 89% HPDI
shade( mu.HPDI , weight.seq )

# prediction samples =======================================

sim.height <- sim(m43, list(weight = weight.seq))
str(sim.height)

height.PI <- apply(sim.height, 2, PI, prob = 0.89)

# plot raw data
plot( height ~ weight , d2 , col=col.alpha(rangi2,0.5) )
# draw MAP line
lines( weight.seq , mu.mean )
# draw HPDI region for line
shade( mu.HPDI , weight.seq )
# draw PI region for simulated heights
shade( height.PI , weight.seq )


# 4.5 polynomial regression ================================
plot(height ~ weight, d)
d[,weightS := scale(weight)]
d[,weightS2 := weightS^2]
plot(height ~ weightS, d)
m45 <- 
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b1*weightS + b2*weightS2,
    a ~ dnorm(178,100),
    b1 ~ dnorm(0,10),
    b2 ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ) %>% 
  map(d)
precis(m45)
cov2cor(vcov(m45))

# parameter uncertainty
weight.seq <- seq( from=-2.2 , to=2 , length.out=30 )
pred_dat <- list( weightS=weight.seq , weightS2=weight.seq^2 )
mu <- link( m45 , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.89 )
# prediction uncertainty
sim.height <- sim( m45 , data=pred_dat )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )


plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) )
lines( weight.seq , mu.mean )
shade( mu.PI , weight.seq )
shade( height.PI , weight.seq )

# converting back to natural scale
plot(height ~ weightS, d, col = col.alpha(rangi2, 0.5), xaxt ='n')
at <- c(-2,-1,0,1,2)
labels <- at*sd(d$weight) + mean(d$weight)
axis( side=1 , at=at , labels=round(labels,1) )


# Practice =====================================================================
# Medium

# 4M1
rnorm(1e4, rnorm(1e4, 0, 10), runif(1e4,0,10))
# 4M2
mm41 <- 
  alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(0, 10),
    sigma ~ dunif(0, 10)
  ) %>% 
  map(d2)

extract.samples(mm41) %>% as.data.table %>% magrittr::extract( , rnorm(1e4, mu,sigma))







# 4H1
m4h1 <-
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b * weight,
    a ~ dnorm(156,100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ) %>% 
  map(d2)

missingPeople <- 
  data.table( individual = 1:5,
              weight = c(46.95, 43.72, 64.78, 32.59, 54.63))

mu <- link(m4h1, missingPeople, 1e4)

missingPeople[ , `expected height` := apply( mu , 2 , mean )]
missingPeople[ , `89% interval` := 
                      apply( mu , 2 , HPDI , prob=0.89 ) %>% 
                      apply(2, round, digits = 2) %>% 
                      apply(2, str_c, collapse = " - ")]

missingPeople

mmm <- lm(height ~ weight, d2)
predict(mmm, missingPeople)


# 4H2
d3 <- d[age < 18]
m4h2 <-
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b * weight,
    a ~ dnorm(108,100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ) %>% 
  map(d3)


weightSeq <- seq(0,45,5)
mu <- link(m4h2, data.frame(weight = weightSeq), 1e4)
muMean <- apply(mu,2,mean)
muHPDI <- apply(mu,2,HPDI)

simHeight <- sim(m4h2, list(weight = weightSeq))
heightPI <- apply(simHeight, 2, PI, prob = 0.89)

plot(height ~ weight, d3, col = col.alpha(rangi2, 0.5))
lines( weightSeq , muMean )
# plot a shaded region for 89% HPDI
shade( muHPDI , weightSeq )
# draw PI region for simulated heights
shade( heightPI , weightSeq )

# 4H3
m4h3 <-
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b * log(weight),
    a ~ dnorm(178,100),
    b ~ dnorm(0,100),
    sigma ~ dunif(0,50)
  ) %>% 
  map(d)
precis(m4h3)


weightSeq <- exp(seq(1.2,4.2,0.25))
mu <- link(m4h3, data.frame(weight = weightSeq), 1e4)
muMean <- apply(mu,2,mean)
muHPDI <- apply(mu,2,HPDI)

simHeight <- sim(m4h3, list(weight = weightSeq))
heightPI <- apply(simHeight, 2, PI, prob = 0.89)

plot( height ~ weight , data=Howell1, col=col.alpha(rangi2,0.4) )
lines( weightSeq , muMean )
# plot a shaded region for 89% HPDI
shade( muHPDI , weightSeq )
# draw PI region for simulated heights
shade( heightPI , weightSeq )





