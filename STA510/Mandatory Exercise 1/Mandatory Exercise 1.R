# Mandatory Exercise 1
# Author: Børge Olav Haug

# Remove old variables:
rm(list=ls())

# If a larger chunk of code is very similar to/adapted from a lecture code
# example, I will note which one I've looked at.

###################################################################################
###################################################################################
# Problem 1
# b)

# Using ppois and pgamma to verify the probability calculations in point a):

# P(N >= 4) = 
cat('Calculation by hand gave approximately 0.3527\nppois in R gives:',
    ppois(q=3, lambda=3, lower.tail=FALSE), '\n')

# P(Xi >= 5) =
cat('Calculation by hand gave approximately 0.5037\npgamma in R gives:',
    pgamma(q=5, shape=2, scale=3, lower.tail=FALSE), '\n')

# Implementation of the algorithm from point a):

# A function that generates total visitor times.
gentotalvisitortimes <- function(Nsim, lambda, alpha, beta) {
  totalvisitortimes <- vector(length=Nsim)
  for(i in range(1, Nsim)) {
    # N <- rpois(n=1, lambda=lambda)
    # Xs <- rgamma(n=N, shape=alpha, scale=beta)
    # V <- sum(Xs)
    totalvisitortimes[i] <-
      sum(rgamma(n=rpois(n=1, lambda=lambda), shape=alpha, scale=beta))
  }
  return(totalvisitortimes)
}

# Using the algorithm to estimate the expected total visitor time for
# visitors arriving during one hour (60 minutes), and to estimate the
# probability that the total visitor time exceeds 1000 minutes.
Nsim <- 10000
t <- 60
lambda <- 3
alpha <- 2
beta <- 3

Vt <- gentotalvisitortimes(Nsim,lambda*t,alpha,beta)

# Expected total visitor time:
cat('Expected total visotor time:', mean(Vt))

# Probability that the total visitor time exceeds 1000 minutes:
cat('Probability that the total visitor time exceeds 1000 minutes:', mean(Vt>1000))

###################################################################################
###################################################################################
# Problem 2
# b)

# Implementation of the algorithm from point a):

# A function that generates values from the Weibull-distribution.
rweibull <- function(Nsim, alpha, beta) {
  # Generate values from the uniform(0,1) distribution
  u <- runif(n=Nsim)
  # Apply the invere F^-1(u) to get values from the Weibull-distribution
  return((-log(1-u)/alpha)^(1/beta))
}

# Using the algorithm to generate a large number of observations from the
# Weibull distribution.
Nsim <-10000
alpha <- 1.83 # My height in meters.
betas <- c(2, 0.5, 1)

# A function for calculating the mean of a Weibull-distributed r.v.
mweibull <- function(alpha, beta) {
  return(alpha^(-1/beta)*gamma(1+1/beta))
}

# A function for calulating the Weibull cdf.
dweibull <- function(x, alpha, beta) {
  return(alpha*beta*x^(beta-1)*exp(-alpha*x^beta))
}

# Allow for 3 histograms on the same image.
# Also useful for creating the histograms in problem 2 d)
par(mfrow = c(3, 1))

for(beta in betas) {
  cat('Generating', Nsim, 'values, with alpha =', alpha, 'and beta = ', beta, '\n')
  values <- rweibull(Nsim, alpha, beta)
  hist(values, main=paste('Observations from a Weibull dist.',
                          '\nwith alpha = ', alpha, 'and beta = ', beta),
       probability=TRUE)
  
  # Superimpose the Weibull cdf to see if the inverse transform method did a good job.
  curve(dweibull(x, alpha, beta), from = 0, to = max(values), add=TRUE, col='red')
  
  cat('Average of the data:', mean(values),
      '\nExpectation calculated from the formula:', mweibull(alpha,beta),'\n\n')
}

###################################################################################
###################################################################################
# Problem 2
# d)

# Implementation of the algorithm in c)
# Note: The implementation was done by looking at the simheads2 funtion
# in sums_mixtures_examples.R and adapting it to this problem. 

# Function for generating a matrix containing Npumps columns (each column
# representing a pump) and Nsim rows. Each entry i,j in the matrix
# represents failure time i for pump j. 
genfailuretimes <- function(Nsim, Npumps, alpha, beta) {
  return(matrix(rweibull(Nsim*Npumps, alpha, beta), ncol=Npumps))
}

Nsim <- 10000
Npumps <- 3
alpha <- 0.08
beta <- 1

# Generate failure times.
failuretimes <- genfailuretimes(Nsim, Npumps, alpha, beta)

# Estimating the expected time until failure for the three scenarios:
# i)
# Operational situation 1: It is sufficient that one pump functions
# (i.e. the system fails when all pumps have failed).
# Pick the last pump to fail, i.e. apply the max function to every row.

systemfailuretimes1 <- apply(failuretimes, MARGIN = 1, FUN = max)
cat('Operational situation 1: It is sufficient that one pump functions',
    '\n(i.e. the system fails when all pumps have failed):',
    '\nEstimate of expected time until failure:', mean(systemfailuretimes1),
    '\nProbability that the pump system fails before one year',
    '\nCalculated using simulation:', mean(systemfailuretimes1<12),
    '\nCalculated by hand: 0.2350','\n\n')

# ii)
# Operational situation 2: All pumps have to function (i.e. the
# system fails when the first pump fails):
# Pick the first pump to fail, i.e. apply the min function to every row.

systemfailuretimes2 <- apply(failuretimes, MARGIN = 1, FUN = min)
cat('Operational situation 2: All pumps have to function',
    '\n(i.e. the system fails when the first pump fails)',
    '\nEstimate of expected time until failure:', mean(systemfailuretimes2),
    '\nProbability that the pump system fails before one year',
    '\nCalculated using simulation:', mean(systemfailuretimes2<12),
    '\nCalculated by hand: 0.9439','\n\n')

# iii)
# Operational situation 3: 2 out of 3 pumps have to function (i.e.
# the system fails when at least 2 of the pumps have failed):
# Pick the second pump to fail, i.e. apply the sort function to every row,
# and pick the 2nd element in every row, i.e. the 2nd column of the matrix.
# The result of apply() is "transposed" and the columns have "become" rows.
# Therefore, we pick the 2nd row instead.

systemfailuretimes3 <- apply(failuretimes, MARGIN = 1, FUN = sort)[2,]
cat('Operational situation 3: 2 out of 3 pumps have to function',
    '(i.e. the system fails when at least 2 of the pumps have failed):',
    '\nEstimate of expected time until failure:', mean(systemfailuretimes3),
    '\nProbability that the pump system fails before one year',
    '\nCalculated using simulation:', mean(systemfailuretimes3<12),
    '\nCalculated by hand: 0.6724','\n\n')

###################################################################################
###################################################################################
# Problem 2
# d)

# Implementation of the algorithm in c)

# Function that uses the acceptance-rejection method to generate numbers from
# a triangle distribution. From my solution to problem 2 in exercise set 3.
# Note: The implementation was done by looking at the code in
# acceptance-rejection_examples.R and adapting it to that problem.
rtriangle <- function(Nsim, xmin, xmode, xmax) {
  k <- 0 # counter for accepted
  x <- numeric(Nsim) # vector for accepted
  while(k<Nsim){
    u <- runif(n=1)
    y <- runif(n=1, min=xmin, max=xmax) # proposal distribution g(y) = unif(xmin, xmax) = 1/(xmax-xmin)
    # The peak in triangle distribution occurs at xmode. f(xmode) = 2/(xmax-xmin)
    # g(y) = 1/(xmax-xmin)
    # So, f(y)/g(y) <= (2/(xmax-xmin))/(1/(xmax-xmin)) = 2/(xmax-xmin)*(xmax-xmin) = 2
    # i.e. c = 2
    fy_div_cgy <- 0
    if(y >= xmin & y <= xmode) {
      # f(y) = 2*(y-xmin)/((xmax-xmin)*(xmode-xmin))
      # f(y)/(c*g(y)) = (y-xmin)/(xmode-xmin)
      fy_div_cgy = (y-xmin)/(xmode-xmin) # 
    } 
    else { # if y > xmode & y <= xmax
      # f(y) = 2*(xmax-y)/((xmax-xmin)*(xmax-xmode))
      # f(y)/(c*g(y)) = (xmax-y)/(xmax-xmode)
      fy_div_cgy = (xmax-y)/(xmax-xmode)
    }
    if(u<fy_div_cgy){ # Accept if u < f(y)/(c*g(y))
      k <- k+1
      x[k] <- y
    }
  }
  return(x)
}

Nsim <- 10000
xmin <- 10
xmode <- 20
xmax <- 50
probs <- seq(0.1,0.9,by=0.1) # (Used to calculate quantiles)

# Generate alpha values.
alphas <- 1/rtriangle(Nsim, xmin, xmode, xmax)

# i)
# Operational situation 1: It is sufficient that one pump functions
# (i.e. the system fails when all pumps have failed):
p1 <- (1-exp(-12*alphas))^3
hist(p1, main='P(The last pump to fail\nfails before one year)',
     probability = TRUE)

cat('Quantiles of the distribution in operational situation 1:')
quantile(p1, probs=probs)

# ii)
# Operational situation 2: All pumps have to function (i.e. the
# system fails when the first pump fails):
p2 <- 1-exp(-36*alphas)
hist(p2, main='P(The first pump to fail\nfails before one year)',
     probability = TRUE)

cat('Quantiles of the distribution in operational situation 2:')
quantile(p2, probs=probs)

# iii)
# Operational situation 3: 2 out of 3 pumps have to function (i.e.
# the system fails when at least 2 of the pumps have failed):
p3 <- (1-exp(-12*alphas))^3 + 3*(1-exp(-12*alphas))^2*exp(-12*alphas)
hist(p3, main = 'P(2 or more pumps fail\nbefore one year)',
     probability = TRUE)

cat('Quantiles of the distribution in operational situation 3:')
quantile(p3, probs=probs)

# Are the new pumps better?
cat('If the estimated alpha parameter(s) of the new pumps is lower than the',
    '\nalpha parameter of the previous pumps, the probability of failure before',
    '\none year of a pump will be lower. This in turn will lead to a lower',
    '\nprobability that the system fails before one year in all of the three', 
    '\noperational situations.',
    '\nOld alpha: 0.08',
    '\nExpected new alpha:', mean(alphas),
    '\nNew alpha lower than old alpha?:', mean(alphas)<0.08, '\n'
)
if(mean(alphas)<0.08) {
  cat('We have reason to believe that the new pumps will have have a lower',
      '\nprobability of failure within one year than the old pumps'
  )
} else {
  cat('We have don\'t have any reason to believe that the new pumps will have',
      'have a lower probability of failure within one year than the old pumps'
  )
}

###################################################################################
###################################################################################
# Problem 3
# b)

# A function for generating upper section Yatzy scores.
genscores <- function(Nsim) {
  scores <- vector(length=Nsim)
  # Repeat the Yatzy-game Nsim times.
  for(i in 1:Nsim) {
    score <- 0
    # Simulate round 1 to 6.
    for(round in 1:6) {
      # Number of dice put aside
      naside <- 0
      # Simulate the 3 throws.
      for(throw in 1:3) {
        # Throw 5 dice minus the ones we've put aside.
        dicethrows <- sample(1:6, size=5-naside, replace=TRUE)
        # Update score.
        score <- score + sum(dicethrows == round)*round
        # Put aside dice if any of them matched what we wanted for this round.
        naside <- naside + sum(dicethrows == round)
        # We can end the round if we've put aside all the dice.
        if(naside == 5)
          break
      }
    }
    scores[i] <- score
  }
  return(scores)
}

e <- 0.01
Nsim <- 1/e^2 # 10000

scores <- genscores(Nsim)

# Allow for the plot to cover the whole image.
par(mfrow = c(1, 1))

barplot(table(scores)/Nsim, main='Yatzy upper section score distribution')

# Estimated probability of a score of at least 42
cat('Probability of getting a bonus, i.e. getting a score of at least 42:', mean(scores>=42))
