library(rethinking)

data("Howell1")

kungAdults <- 
  Howell1 %>%
  as.data.table() %>% 
  .[age >= 18]

# 1 ============================================================================
# The weights listed below were recorded in the !Kung census, but heights
# were not recorded for these individuals. Provide predicted heights and 89%
# compatibility intervals for each of these individuals. That is, fill in the table
# below, using model-based predictions.

markedWgts <- c(45,40,65,31,53)
q1Dt <- data.table(weight = 25:75)

meanWgt <- kungAdults[,mean(weight)]

modelSpec <- 
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - meanWgt),
    a ~ dnorm(178,20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  )

model1 <- quap(modelSpec,kungAdults)

post <- extract.samples(model1)

f <- function(wgt,model = model1){
  model@coef['a'] + model@coef['b'] * (wgt - meanWgt)
}

# predictions
q1Dt[ , height_hat := map_dbl(weight,f)]

# compatability intervals
muSamples <- link(model1,q1Dt[,.(weight)])
bands <- 
  apply(muSamples,2,HPDI, prob = 0.89) %>%
  t %>% 
  as.data.table %>%
  .[ , weight := 25:75]

q1Dt <- q1Dt[bands, on = 'weight']
q1Dt[ , Marked := ifelse(weight %in% markedWgts,TRUE,FALSE)]


ggplot(q1Dt, aes(x = weight,y = height_hat)) +
  geom_ribbon(aes(ymin = lower89, ymax = upper89), fill = 'grey70') +
  geom_line() +
  geom_point(data = q1Dt[(Marked)], aes(col = 'red'))


# 2 ============================================================================
# Model the relationship between height (cm) and the natural logarithm of
# weight (log-kg): log(weight). Use the entire Howell1 data frame, all 544
# rows, adults and non-adults. Use any model type from Chapter 4 that you
# think useful: an ordinary linear regression, a polynomial or a spline. Plot
# the posterior predictions against the raw data.

# We'll use splines

modelSpec2 <- alist(
  log(height) ~ dnorm(mu)
  
  
  
)







# 3 ============================================================================
# Plot the prior predictive distribution for the polynomial regression model
# in Chapter 4. You can modify the the code that plots the linear regression
# prior predictive distribution. 20 or 30 parabolas from the prior should suffice
# to show where the prior probability resides. Can you modify the prior
# distributions of alpha, beta_1, and beta_2 so that the prior predictions stay within the
# biologically reasonable outcome space? That is to say: Do not try to fit the
# data by hand. But do try to keep the curves consistent with what you know
# about height and weight, before seeing these exact data.






































