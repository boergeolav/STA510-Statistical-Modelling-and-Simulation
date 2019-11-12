rm(list=ls()) # clear all variables
graphics.off() # clear all figures
cat("\014") # clear console window

###################################################################
# READ ME!

# You should run the whole script once before
# running code snippets separately because some
# may depend on variables and functions defined above it.

###################################################################
##################### Problem 1c ###################################
###################################################################

cat('Problem 1c \n')

# Simulates Nsim individual bet profits.
simPs <- function(Nsim, winprobs, S) {
  X <- matrix(nrow=length(winprobs), ncol=Nsim) # Game outcomes Xi
  for (i in (1:length(winprobs))) {
    p <- winprobs[i]
    X[i,] <- sample(x=c(1, 0), size=Nsim, replace=TRUE, prob=c(p, 1-p))
  }
  TO <- 1/winprobs # Odds TOi = 1/pi
  
  # Ps = (sum(S*TOi*Xi) - n*S) = S*(sum(TOi*Xi - 1))
  profits <- S*colSums(TO*X - 1)
  
  return(profits)
}

Nsim <- 10000
winprobs <- c(0.67, 0.38, 0.4, 0.54)
stake <- 100

Ps <- simPs(Nsim, winprobs, stake)

cat('Estimated expected profit when placing individual bets:', mean(Ps), '\n')

cat('\n')

###################################################################
##################### Problem 1e ###################################
###################################################################
cat('Problem 1e \n')

# Simulates Nsim combined bet profits.
simPc <- function(Nsim, winprobs, S) {
  Xc <- sample(x=c(1, 0), size=Nsim, replace=TRUE,
               prob=c(prod(winprobs), 1-prod(winprobs))) # Combined game outcome Xc
  n <- length(winprobs) # Number of games
  TOc <- prod(1/winprobs) # Combined odds TOc = prod(1/pi)
  
  # Pc = n*S*Xc*TOc - n*S = n*S*(TOc*Xc - 1)
  profits <- n*S*(TOc*Xc - 1)
  return(profits)
}

# NB! Makes use of Nsim, winprobs, and stake defined in the 1c code section!
Pc <- simPc(Nsim, winprobs, stake)

# NB! Makes us of Ps in the 1c code section!
cat('Estimated expected profit when placing individual bets:', mean(Ps), '\n')
cat('Estimated expected profit when combining bets:', mean(Pc), '\n')
cat('Estimated standard deviation when placing individual bets:', sd(Ps), '\n')
cat('Estimated standard deviation when combining bets:', sd(Pc), '\n')

cat('Ps and Pc are both centered around 0, but since Pc has greater variance',
    'more simulations are required to reliably get mean(Pc) close to 0.\n')

cat('\n')

###################################################################
##################### Problem 2b ###################################
###################################################################

cat('Problem 2b \n')

# The first intensity function in problem 2.
lambda <- function(t) {
  return(5 + 50*sin(pi*t/24)^2 + 190*exp(-(t-20)^2/3))
}

# Uses the hit or miss method to estimate the integral of f(t)
# between a and b.
HMMCintegral <- function(n, a, b, c, f) {
  X <- runif(n, a, b)
  Y <- runif(n, 0, c)
  return(c*(b-a)*mean(Y<=f(X)))
}

# Calculates absolute distances between pairs of consecutive elements in
# a vector, and returns the max distance multiplied by the given k.
compensate1 <- function(vec, k) {
  d <- numeric(length(vec)-1)
  for(i in 1:length(d)) {
    d[i] <- abs(vec[i+1]-vec[i])
  }
  return(k*max(d))
}

# Calculates the distance between the max of a vector and its two
# neighbours, and returns the max of these distances multiplied by some
# constant k.
compensate2 <- function(vec, k) {
  i <- which(vec==max(vec))
  if(i<length(vec)) {
    return(k*max(c(abs(vec[i]-vec[i-1])), abs(vec[i+1]-vec[i])))
  } else {
    return(0)
  }
}

### Simulate an estimate of the integral with the number of simulations
### required to have the precision described in point a)
a <- 0
b <- 24
tvec <- seq(from=a, to=b, by = 0.01)
lambdavec <- lambda(tvec)

# Add a small number to the max in case we missed a larger lambda(t)-value.
c <- max(lambdavec) + compensate2(lambdavec, 2)

alpha2b <- 0.05
z <- round(qnorm(p=alpha2b/2, lower.tail=FALSE), digits=2)
e <- 10

n <- ceiling(c^2*(b-a)^2*z^2*0.25/e^2)
cat('Required number of simulations to be ', (1-alpha2b)*100,
    '% certain that the estimate is no more than ', e, 
    ' from the true answer: ', n, '\n', sep='')

lambda.HM <- HMMCintegral(n, a, b, c, lambda)

cat('An estimate of the integral:', lambda.HM, '\n')


### Calculate the probability of having more than 1250 visitors during one day
### (24 hours).
cat('Probability of having more than 1250 visitors during one day',
    ppois(q=1250, lambda=lambda.HM, lower.tail=FALSE), '\n'
)

cat('\n')
###################################################################
##################### Problem 2d ###################################
###################################################################
cat('Problem 2d \n')

# Copied from stochastic_processes_examples.R, and adapted to this problem.
simtNHPP <- function(a,b,lambdamax,lambdafunc){
  # Simple check that a not too small lambdamax is set
  if(max(lambdafunc(seq(a,b,length.out = 100)))>lambdamax)
    stop("lambdamax is smaller than max of the lambdafunction")
  # First simulate HPP with intensity lambdamax on a to b
  expectednumber <- (b-a)*lambdamax  
  Nsim <- 3*expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
  interarrivals <- rexp(Nsim,lambdamax) # Simulate interarrival times
  arrivals <- a+cumsum(interarrivals)   # Calculate arrival times starting at a
  arrivals <- arrivals[arrivals<b] # Dischard the times larger than b
  Nevents <- length(arrivals) # Count the number of events
  # Next do the thinning. Only keep the times where u<lambda(s)/lambdamax
  U <- runif(Nevents)
  arrivals <- arrivals[U<lambdafunc(arrivals)/lambdamax]
  deleted <- 1 - length(arrivals)/Nevents # Proportion of deleted arrival times
  return(list(deleted=deleted, arrivals=arrivals))  # Return the proportion of deleted times, and the remaining times
}

# NB! Makes use of a, b, c, and lambda defined in the 2b code section!
dummy <- simtNHPP(a, b, c, lambda)
deleted2d <- dummy$deleted
arrivals2d <- dummy$arrivals

plot(arrivals2d, 1:length(arrivals2d), type='s', xlab = 'time (hours)', 
     ylab = 'Event number', lwd=1.5,
     main='Data for arrival of visitors during one day')

### Verify by simulations the calculation in point c) of the
### proportion of arrival times being deleted in the thinning
### step of the algorithm.
cat('Proportion of arrival times actually deleted: ',
    round(deleted2d, digits=2)*100.0, '%\n', sep='')

# NB! Makes use of a, b, c, and lambda.HM defined in the 2b code section!
deletedcalc <- 1-lambda.HM/((b-a)*c)
cat('Proportion of arrival times being deleted, calculated as in c): ',
    round(deletedcalc, digits=2)*100.0,
    '%\n', sep='')

if(abs(deleted2d-deletedcalc) < 0.05) {
  cat('Calculated and actual proportion are pretty close!\n')
} else {
  cat('Calculated and actual proportion are not very close!\n')
}

# Copied from stochastic_processes_examples.R.
calculatequeue <- function(arrivaltimes, servicetimes){
  Narrivals <- length(arrivaltimes) # Total number of arrival
  departuretimes <- sort(arrivaltimes+servicetimes) # Calculate and sort the departure times
  eventtimes <- 0 # This will be the vector for event times
  numbersinqueue <- 0 # This will be the vector for numbers in queue, updated at each event time
  currentnumber <- 0  # Keeps track of the current number of customers
  acounter <- 1  # Counter for the arrivals vector
  dcounter <- 1  # Counter for the departures vector
  while(acounter<=Narrivals){
    if(arrivaltimes[acounter]<departuretimes[dcounter]){ # If the next event is an arrival
      currentnumber <- currentnumber+1
      numbersinqueue <- c(numbersinqueue,currentnumber)
      eventtimes <- c(eventtimes,arrivaltimes[acounter])
      acounter <- acounter+1
    }
    else{  # If the next event is an departure
      currentnumber <- currentnumber-1
      numbersinqueue <- c(numbersinqueue,currentnumber)
      eventtimes <- c(eventtimes,departuretimes[dcounter])
      dcounter <- dcounter+1
    }
  }
  return(list(numbers=numbersinqueue,times=eventtimes))
}

### Simulate and plot data for the number of active visitors
### over time for one day.
alpha <- 2
beta <- 3
n <- length(arrivals2d)

servicetimes <- rgamma(n=n, shape=alpha, scale=beta)/60.0 # Convert minutes -> hours
dummy <- calculatequeue(arrivals2d, servicetimes)
times <- dummy$times
numbers <- dummy$numbers

plot(times, numbers, type='s', xlab='time (hours)', 
     ylab='Number of active visitors', lwd=1.5,
     main='Number of active visitors over time for one day')

### Estimate the probability that the maximum number of visitors
### during a day exceeds 30.
### Estimate the median and the 5% and the 95% quantiles of the
### number of visitors at time t = 12.
# Copied from stochastic_processes_examples.R, and adapted to this problem.
t <- 12
Nsim <- 100
maxvalues <- numeric(Nsim)
tvalues <- numeric(Nsim)
for(i in 1:Nsim){
  arrivaltimes <- simtNHPP(a, b, c, lambda)$arrivals
  servicetimes <- rgamma(length(arrivaltimes), shape=alpha, scale=beta)/60.0 # minutes -> hours
  novertime <- calculatequeue(arrivaltimes,servicetimes)
  maxvalues[i] <- max(novertime$numbers)
  # We are not guaranteed that there's an event exactly at t, so let's
  # pick the average of an interval around t.
  # Start with a small interval and make it wider if necessary.
  tvalues[i] <- NaN
  k <- 1
  delta <- 0.01
  while(is.nan(tvalues[i])) {
    index <- abs(novertime$times-t) < k*delta
    # index <- novertime$times>t-k*delta & novertime$times<t+k*delta
    tvalues[i] <- mean(novertime$numbers[index])
    k <- k + 1
  }
}

cat('Probability that the maximum number of visitors during a day ',
    'exceeds 30: ', round(mean(maxvalues>30), digits=3)*100.0, '%\n', sep='')

cat('Median of the number of visitors at time t = ', t , ': ', median(tvalues),
    '\n', sep='')

quantiles <- quantile(tvalues, probs=c(0.05, 0.95))
cat('5% quantile of the number of visitors at time t = ', t , ': ',
    quantiles[1], '\n', sep='')
cat('95% quantile of the number of visitors at time t = ', t , ': ',
    quantiles[2], '\n', sep='')

cat('\n')

###################################################################
##################### Problem 2f ###################################
###################################################################
cat('Problem 2f \n')

# The second intensity function in problem 2, from point e).
lambda2 <- function(t) {
  return(10 + 10*t)
}

# The inverse of the big lambda
lambda2inverse <- function(w) {
  return((-10+sqrt(100+20*w))/10)
}

# Copied from stochastic_processes_examples.R, and adapted to this problem.
# Simulates arrival times for a NHPP between a and b using the
# transformation method.
# NB! Makes use of compensate2 and HMMCintegral defined in the 2b code section!
simtNHPP2 <- function(a, b, lambdafunc, lambdainversefunc) {
  tvec <- seq(a, b, by=0.01)
  c <- max(lambdafunc(tvec)) + compensate2(lambdafunc(tvec), 2)
  # Expected number of events until stop time, calculated using
  # the hit or miss method. Other ways to do this are to let the
  # expeceted value be one of this function's inputs (afterall,
  # the function 10 + 10t is easily integrated), or to give the
  # BigLambda function as input, and calculate the value as BigLambda(b).
  expectednumber <- HMMCintegral(10000, a, b, c, lambdafunc)
  Nsim <- 3*expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
  interarrivals <- rexp(Nsim,1) # Simulate interarrival times from exp(1)
  arrivals <- a+lambdainversefunc(cumsum(interarrivals)) # Calculate arrival times
  arrivals <- arrivals[arrivals<b] # Discard the times larger than b
  return(arrivals)
}

arrivals2f <- simtNHPP2(a, b, lambda2, lambda2inverse)

plot(arrivals2f, 1:length(arrivals2f), type='s', xlab = 'time (hours)', 
     ylab = 'Event number', lwd=1.5,
     main='Using the transformation method')

cat('\n')


###################################################################
##################### Problem 2g ###################################
###################################################################
cat('Problem 2g \n')

### Plot the two intensity functions lambda(t) and lambda2(t)
# NB! Makes use of tvec, and lambda defined in the 2b code section, 
# and lambda2 in the 2f code section!
plot(tvec, lambda(tvec), xlab='t', ylab='f(t)', col='green', type='l',
     main='Plot of the two intensity functions.')
points(tvec, lambda2(tvec), col='red', type='l')
legend('topleft',c('lambda(t)', 'lambda2(t)'), 
       col=c('green', 'red'),lty=1)

# Copied from stochastic_processes_examples.R, and adapted to this problem.
# Simulates arrival times for a NHPP between a and b using the
# alternative thinning algorithm.
# NB! Makes use of HMMCintegral and compensate2 defined in the 2b code section,
# and simtNHPP2 defined in the 2d code section!
simtNHPP3 <- function(a, b, lambdaupper, lambdaupperinverse, lambdafunc) {
  # First simulate NHPP with intensity lambdaupper(t) on a to b
  tvec <- seq(a, b, by=0.01)
  c <- max(lambdafunc(tvec)) + compensate2(lambdafunc(tvec), 2)
  expectednumber <- HMMCintegral(10000, a, b, c, lambdafunc)
  Nsim <- 3*expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
  arrivals <- simtNHPP2(a, b, lambdaupper, lambdaupperinverse) # Simulate an NHPP with intensity lambdaupper(t)
  arrivals <- arrivals[arrivals<b] # Discard the times larger than b
  Nevents <- length(arrivals) # Count the number of events
  # Next do the thinning. Only keep the times where u<lambda(s)/lambdaupper(s)
  U <- runif(Nevents)
  arrivals <- arrivals[U<lambdafunc(arrivals)/lambdaupper(arrivals)]
  deleted <- 1 - length(arrivals)/Nevents # Proportion of deleted arrival times
  return(list(deleted=deleted, arrivals=arrivals))  # Return the proportion of deleted times, and the remaining times
}

# We can verify that the algorithm seems to be correct by comparing
# the plot of this thinning algorithm to the plot of the
# standard thinning algorithm from point 2b), and check if they
# are similar.
for(i in 1:3) {
  # NB! Makes use of a, b, and lambda defined in the 2b code section,
  # and lambda2 and lambda2inverse from the 2f code section.
  dummy <- simtNHPP3(a, b, lambda2, lambda2inverse, lambda)
  deleted2g <- dummy$deleted
  arrivals2g <- dummy$arrivals
  plot(arrivals2g, 1:length(arrivals2g), type='s', xlab = 'time (hours)', 
       ylab = 'Event number', lwd=1.5, col='blue',
       main=paste('Using the alternative thinning algorithm, #', i, sep=''))
  # NB! Makes use of arrivals2d defined the 2d code section!
  points(arrivals2d, 1:length(arrivals2d), type='s',
         lwd=1.5, col='red')
  legend('topleft',c('standard thinning', 'alternative thinning'), 
         col=c('blue', 'red'),lty=1)
}
# The curves for the two thinning methods should be similar in shape.


# NB! Makes use of deleted2d defined in the 2d code section!
cat('Proportion of arrival times actually deleted in point d): ',
    round(deleted2d, digits=2)*100.0, '%\n', sep='')
cat('Proportion of arrival times actually deleted in point g): ',
    round(deleted2g, digits=2)*100.0, '%\n', sep='')
cat('Since lambda2(t) is closer to lambda(t) than lambdamax is,',
    'we don\'t have to delete as many arrival times as in the standard',
    'thinning algorithm.\n')

cat('\n')


