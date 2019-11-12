rm(list=ls()) # clear all variables
graphics.off() # clear all figures
cat('\014') # clear console window

###################################################################

# Several functions are vectorized with Vectorize() so that they
# also accept a combination of scalars and vectors
# (vectors must be of equal length).

###################################################################
##################### Problem 1b ###################################
###################################################################
cat('Problem 1b \n')

# pdf of the triangle distribution.
dtriangle <- Vectorize(function(x, a, b, c) {
  if((a <= x) && (x <= c)) {
    return(2*(x-a)/((b-a)*(c-a)))
  } else if((c <= x) && (x <= b)) {
    return(2*(b-x)/((b-a)*(b-c)))
  } else {
    return(0.0)
  }
})

# cdf of the triangle distribution.
# NB! Not used anywhere.
cdftriangle <- Vectorize(function(x, a, b, c) {
  if((x <= a)) {
    return(0.0)
  } else if((a < x) && (x <= c)) {
    return((x-a)^2/((b-a)*(c-a)))
  } else if((c < x) && (x <= b)) {
    return(1-(b-x)^2/((b-a)*(b-c)))
  } else { # if(x >= b)
    return(1.0)
  }
})

# Inverse of the cdf of the triangle distribution.
cdfinvtriangle <- Vectorize(function(u, a, b, c) {
  if(u == 0.0) {
    return(a)
  } else if((0 < u) && (u <= (c-a)/(b-a))) {
    return(a+sqrt(u*(b-a)*(c-a)))
  } else if(((c-a)/(b-a) < u) && (u < 1)) {
    return(b-sqrt((1-u)*(b-a)*(b-c)))
  } else { # if(u == 1.0)
    return(b)
  }
})

# Generates numbers from the triangle distribution.
rtriangle <- function(n, a, b, c) {
  u <- runif(n)
  return(cdfinvtriangle(u, a, b, c))
}

# Generate triangle distribution data, and display in
# a histogram with the pdf superimposed.
n.1b <- 10000
a.1b <- 0.7
b.1b <- 2
c.1b <- 1.5

data <- rtriangle(n.1b, a.1b, b.1b, c.1b)

par(mfrow=c(1,1))
hist(x=data, xlab='x', ylim=c(0, dtriangle(c.1b, a.1b, b.1b, c.1b)+0.1),
     probability=TRUE, main='Triangle data with pdf superimposed.')
curve(dtriangle(x, a.1b, b.1b, c.1b), from=a.1b, to=b.1b, add=TRUE,
      col='red', type='l')

cat('\n')

###################################################################
##################### Problem 2b ###################################
###################################################################
cat('Problem 2b \n')

# Run once: install.packages('mvtnorm')
library(mvtnorm)

n.2b <- 10000
muvec <- c(75, 46, 18)
par(mfrow=c(3,1))

# i)
sigmamat1 <- matrix(c(625, -187.5, 0, -187.5, 100, 0, 0, 0, 25), nrow=3)
X1 <- rmvnorm(n.2b, muvec, sigmamat1)
p1 <- mean((X1[,1]>80) & (X1[,2]>50) & (X1[,3]>20))
cat('P(X1 > 80, X2 > 50, X3 > 20) in scenario i):',
    p1, '\n'
    )
cat('P(X1 + X2 + X3 > 150) in scenario i):',
    mean((X1[,1]+X1[,2]+X1[,3])>150), '\n'
)

# ii)
sigmamat2 <- matrix(c(625, 187.5, 0, 187.5, 100, 0, 0, 0, 25), nrow=3)
X2 <- rmvnorm(n.2b, muvec, sigmamat2)
p2 <- mean((X2[,1]>80) & (X2[,2]>50) & (X2[,3]>20))
cat('P(X1 > 80, X2 > 50, X3 > 20) in scenario ii):',
    p2, '\n'
)
cat('P(X1 + X2 + X3 > 150) in scenario ii):',
    mean((X2[,1]+X2[,2]+X2[,3])>150), '\n'
)

# iii)
sigmamat3 <- matrix(c(625, -187.5, 50, -187.5, 100, -25, 50, -25, 25), nrow=3)
X3 <- rmvnorm(n.2b, muvec, sigmamat3)
p3 <- mean((X3[,1]>80) & (X3[,2]>50) & (X3[,3]>20))
cat('P(X1 > 80, X2 > 50, X3 > 20) in scenario iii):',
    p3, '\n'
)
cat('P(X1 + X2 + X3 > 150) in scenario iii):',
    mean((X3[,1]+X3[,2]+X3[,3])>150), '\n'
)

cat('Does the prediction \"ii) more likely to generate a top',
    'tier character than i), which is more likely to generate a top tier',
    'character than iii)\", hold?\n')
if((p2 >= p1) & (p1 >= p3)) {
  cat('Yes, it does!\n')
} else {
  cat('No? Maybe it does, but it didn\'t occur this time.\n')
}

# At first glance, one might think that X1 > 80, X2 > 50, X3 > 20
# and X1 + X2 + X3 > 150 means the same thing, but they don't.
# X1 + X2 + X3 > 150 must be fulfilled for X1 > 80, X2 > 50, X3 > 20
# to be true, but it does not go the other way. For example,
# say that X3 > 20 is false. X1 + X2 + X3 > 150 will still be true
# if X1 + X2 is large enough to push the total sum to over 150.

cat('\n')

###################################################################
##################### Problem 3b ###################################
###################################################################
cat('Problem 3b \n')

# The function we want to integrate.
lambda <- function(t) {
  return(5+50*sin(pi*t/24)^2+190*exp(-(t-20)^2/3))
}

# Calculating the required number of simulations.
alpha.3b <- 0.05
z.3b <- qnorm(alpha.3b/2, lower.tail=FALSE)
a.3b <- 0
b.3b <- 24
e.3b <- 10

nrep.3b <- 10000
var.hat <- var(lambda(runif(nrep.3b, a.3b, b.3b)))

n.3b <- ceiling(z.3b^2*(b.3b-a.3b)^2/e.3b^2*var.hat)

cat('We need', n.3b, 'simulations to be at least', (1-alpha.3b)*100,
    '% certain that the estimated integral is no more than', e.3b,
    'from the true answer.\n')

# Approximating the integral.
theta.CMC <- (b.3b-a.3b)*mean(lambda(runif(n.3b, a.3b, b.3b)))

cat('Using crude Monte Carlo integration to approximate the integral',
    'gives:', theta.CMC, '\n')

cat('\n')
###################################################################
##################### Problem 3d ###################################
###################################################################
cat('Problem 3d \n')

a.3d <- a.3b
b.3d <- b.3b
tvec <- seq(0, 24, by=0.01)
c.3d <- tvec[which.max(lambda(tvec)==max(lambda(tvec)))]
n.3d <- 10000

X <- rtriangle(n.3d, a.3d, b.3d, c.3d)

# Approximating the integral using the importance sampling method.
theta.IS <- mean(lambda(X)/dtriangle(X, a.3d, b.3d, c.3d))
cat('Using importance sampling to approximate the integral',
    'gives:', theta.IS, '\n')

# Estimating and comparing the standard deviations of
# ordinary Monte Carlo and importance sampling.
nrep.3d <- 1000
thetaCMCvec <- numeric(n.3d)
thetaISvec <- numeric(n.3d)
for(i in 1:nrep.3d) {
  thetaCMCvec[i] <- (b.3b-a.3b)*mean(lambda(runif(nrep.3d, a.3b, b.3b)))
  X <- rtriangle(nrep.3d, a.3d, b.3d, c.3d)
  thetaISvec[i] <- mean(lambda(X)/dtriangle(X, a.3d, b.3d, c.3d))
}

cat('Standard deviation when using ordinary Monte Carlo: ', sd(thetaCMCvec), '\n')
cat('Standard deviation when using importance sampling: ', sd(thetaISvec), '\n')

# My choice of importance function didn't really improve anything...


cat('\n')

###################################################################
##################### Problem 4b ###################################
###################################################################
cat('Problem 4b \n')

# Code from bootstrap_examples.R adapted to fit this problem.

# Include the data in 'mandatory3_data.R' 
source('./mandatory3_data.R')

MLEest0 <- 1/mean(ntimesC1)
cat('Nucleation rate for the data, estimated by MLE:', MLEest0, '\n')

# The bootstrap procedure.
B <- 5000
n <- length(ntimesC1)
MLEestB <- numeric(B)
for(b in 1:5000) {
  MLEestB[b] <- 1/mean(sample(x=ntimesC1, size=n, replace=TRUE))
}

bias <- mean(MLEestB)-MLEest0
cat('Estimated bias of the MLE:', bias, '\n')
cat('Bias adjusted MLE estimate:', MLEest0-bias, '\n')
cat('Estimated standard deviation of the MLE:', sd(MLEestB), '\n')


# Bootstrap confidence intervals for the nucleation rate.
alpha.4b <- 0.05

# Percentile bootstrap interval (without using boot).
lowerP <- quantile(MLEestB, alpha.4b/2)
upperP <- quantile(MLEestB, (1-alpha.4b/2))
cat('Percentile bootstrap interval for the nucleation rate:\n',
    '[', lowerP, ' , ', upperP, ']\n', sep='')

# Percentile and BCa interval (using boot).
library(boot)
MLE <- function(data,i) {
  return(1/mean(data[i]))
}
boot.obj <- boot(data=ntimesC1, statistic=MLE, R=5000)
boot.obj
boot.ci(boot.obj,type=c('perc', 'bca'))

cat('\n')

###################################################################
##################### Problem 4d ###################################
###################################################################
cat('Problem 4d \n')

# Code from permutation_examples.R adapted to fit this problem.

# The permutation test.
P <- 10000
nA <- length(ntimesC1)
nB <- length(ntimesC2)
m <- nA+nB

MLEdiff0 <- 1/mean(ntimesC1) - 1/mean(ntimesC2)
MLEdiffP <- numeric(P)

for(p in 1:P) {
  permutation <- sample(x=c(ntimesC1, ntimesC2), size=m, replace=FALSE)
  MLEdiffP[p] <- 1/mean(permutation[1:nA]) - 1/mean(permutation[(nA+1):m])
}

alpha.4d <- 0.05
pval.4d <- mean(abs(MLEdiffP) > abs(MLEdiff0))

if(pval.4d < alpha.4d) {
  cat('Reject the hypothesis that there is NO difference in',
      'rate at the two different temperatures.\n')
} else {
  cat('Keep the hypothesis that there is NO difference in',
      'rate at the two different temperatures.\n')
}

cat('\n')

###################################################################
##################### Problem 4f ###################################
###################################################################
cat('Problem 4f \n')

# The bootstrap procedure.
B <- 5000
n <- length(ntimesC1)
MLEestB <- numeric(B)
for(b in 1:5000) {
  MLEestB[b] <- 1/mean(sample(x=ntimesC1, size=n, replace=TRUE))
}

alpha.4f <- 0.05
pval.4f <- mean(MLEestB<=0.0003)
pval.4f

if(pval.4f < alpha.4f) {
  cat('Reject the hypothesis that the rate',
      'is smaller than or equal to 3\n')
} else {
  cat('Keep the hypothesis that the rate',
      'is smaller than or equal to 3\n')
}


cat('\n')


###################################################################
##################### Problem 5b ###################################
###################################################################
cat('Problem 5b \n')

# Calculates the total production up to time s.
G <- Vectorize(function(s, t0, ts, tp, b, y) {
  if((0 <= s) & (s <= t0)) {
    return(0.0)
  } else if((t0 < s) & (s <= t0 + ts)) {
    return(b/(2*ts)*(s-t0)^2)
  } else if((t0+ts < s) & (s <= t0 + ts + tp)) {
    return(b/2*ts+b*(s-(ts+t0)))
  } else if(s > t0 + ts + tp) {
    return(b/2*ts+b*tp+b/y*(1-exp(-y*(s-(t0+ts+tp)))))
  } else {
    stop('s < 0, invalid input')
  }
})

# Calculating total production up to time s,
# for several different scenarios.
svec <- rep(c(5, 10, 15), 2)
t0vec <- c(rep(1, 3), rep(0.75, 3))
tsvec <- c(rep(1.5, 3), rep(1, 3))
tpvec <- c(rep(6, 3), rep(4, 3))
bvec <- c(rep(8, 3), rep(8.5, 3))
yvec <- c(rep(0.15, 3), rep(0.25, 3))

for(i in 1:6) {
  s <- svec[i]; t0 <- t0vec[i]; ts <- tsvec[i]
  tp <- tpvec[i]; b <- bvec[i]; y <- yvec[i]
  cat('Scenario ', i, ': ', 's=', s, ', t0=', t0, ', ts=', ts,
      ', tp=', tp, ', beta=', b, ', gamma=', yvec[i],
      ', G(s)=', G(s, t0, ts, tp, b, y), '\n\n', sep='')
}

cat('\n')

###################################################################
##################### Problem 5d ###################################
###################################################################
cat('Problem 5d \n')

# Simulates the distribution of total production up to time s.
simTotalProd <- function(n, s) {
  t0vec <- rtriangle(n=n, a=0.85, b=1.5, c=1.1)
  tsvec <- rtriangle(n=n, a=0.7, b=1.7, c=1.0)
  tpvec <- rtriangle(n=n, a=4.0, b=7.0, c=5.0)
  bvec <- rtriangle(n=n, a=7.5, b=8.5, c=8.0)
  yvec <- rtriangle(n=n, a=0.15, b=0.3, c=0.25)
  Gvec <- G(s, t0vec, tsvec, tpvec, bvec, yvec)
  return(Gvec)
}

# Simulating the distribution of 5, 10 and 15 years total production.
n.5d <- 10000
par(mfrow=c(3,1))
svec <- c(5, 10, 15)
for(s in svec) {
  Gsdist <- simTotalProd(n.5d, s)
  summ <- c(summary(Gsdist), Sd=sd(Gsdist), Var=var(Gsdist))
  cat('Summary statistics for the distribution of',
      s, 'years total production:\n')
  print(summ)
  cat('\n')
  hist(x=Gsdist, probability=TRUE, 
       xlab=paste('Total production, G(', s, ')', sep=''),
       main=paste(s, 'y total production distribution', sep='')
       )
}

# Simulating the probability that the 15 years total production exceeds
# 80 with the required precision.
alpha.5d <- 0.05
z.5d <- qnorm(p=alpha.5d/2, lower.tail=FALSE)
e.5d <- 0.02
n.5d <- ceiling(z.5d^2*0.25/(e.5d^2))

G15dist <- simTotalProd(n.5d, 15)
cat('P(G(15)>80) = ', mean(G15dist>80), 'with an error of at most',
    e.5d, 'with', (1-alpha.5d)*100, '% certainty.\n')
cat('\n')

###################################################################
##################### Problem 5f ###################################
###################################################################
cat('Problem 5f \n')

# Simulates times to threshold production level and total volumes
# produced until threshold production level is reached.
simTotalProdThresh <- function(n) {
  t0vec <- rtriangle(n=n, a=0.85, b=1.5, c=1.1)
  tsvec <- rtriangle(n=n, a=0.7, b=1.7, c=1.0)
  tpvec <- rtriangle(n=n, a=4.0, b=7.0, c=5.0)
  bvec <- rtriangle(n=n, a=7.5, b=8.5, c=8.0)
  yvec <- rtriangle(n=n, a=0.15, b=0.3, c=0.25)
  scvec <- log(bvec)/yvec + t0vec + tsvec + tpvec
  Gvec <- G(scvec, t0vec, tsvec, tpvec, bvec, yvec)
  return(list(sc=scvec, Gsc=Gvec))
}

# Calculating the required number of simulations.
nrep.5f <- 10000
alpha.5f <- 0.05
z.5f <- qnorm(p=alpha.5f/2, lower.tail=FALSE)
e.sc <- 1/12
e.Gsc <- 0.25
test <- simTotalProdThresh(nrep.5f)
var.sc <- var(test$sc)
var.Gsc <- var(test$Gsc)
n.sc <- ceiling(z.5f^2*var.sc/e.sc^2)
n.Gsc <- ceiling(z.5f^2*var.Gsc/e.Gsc^2)

# Simulating distributions with the number of simulations needed to
# achieve the wanted precision.
scdist <- simTotalProdThresh(n.sc)$sc
Gscdist <- simTotalProdThresh(n.Gsc)$Gsc

cat('Summary statistics for the distribution of time',
    'until threshold production level:\n')
scquantiles <- quantile(scdist, probs = c(0.1, 0.9))
c(scquantiles, summary(scdist), Sd=sd(scdist), Var=var(scdist))

cat('Summary statistics for the distribution of total volume',
    'produced until threshold production level:\n')
Gscquantiles <- quantile(Gscdist, probs = c(0.1, 0.9))
c(Gscquantiles, summary(Gscdist), Sd=sd(Gscdist), Var=var(Gscdist))

par(mfrow=c(2,1))
hist(x=scdist, probability=TRUE, 
     xlab='sc',
     main='Distribution of time until\nthreshold production level'
)

hist(x=Gscdist, probability=TRUE, 
     xlab='G(sc)',
     main='Distribution of total volume produced\nuntil threshold production level'
)

cat('The quantiles tell us that 80 % of the time, the time until threshold production level will lie between',
    scquantiles[1], 'and', scquantiles[2], 'and the total volume produced (until threshold time) will lie between',
    Gscquantiles[1], 'and', Gscquantiles[2], '\n')

cat('\n')