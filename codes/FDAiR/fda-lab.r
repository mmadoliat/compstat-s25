################################################################################
################################################################################
#                               LAB 1                                          #
################################################################################
################################################################################

# This lab is an introduction to the R package 'fda'. Here we will focus on 
# defining functional objects, smoothing methods, manipulating functional
# data objects and some basis exploratory data analysis. 

### Firstly, let's load up the library and some data. 

library('fda')

# All of our examples will come from the Canadian Weather Data. Not all of these
# examples will therefore be particularly realistic. But they are there in order
# to demonstrate aspects of the code. Realism is retained as much as possible. 

# We can load these data by

data(CanadianWeather)

# Temperature and precipitation are contained in the dailyAV element, and we'll
# extract them. 

temp = CanadianWeather$dailyAv[,,1]
precip = CanadianWeather$dailyAv[,,2]

# We need corresponding time points. I'll put them half-way through a day. This
# is because the period is over 0:365 and we'd like the 365 data points to be
# about symmetric in that period.

daytime = (1:365)-0.5

# This is a bit fine for plotting purposes, so we'll also create a vector of
# points to use for plotting every 5 days

day5 = seq(0,365,5)


### Plot these data

# Note: I am assuming you will use only one plotting window and won't close it
# so I haven't insisted on re-setting the window as being two-pane below. We 
# revert to a single pane graphical window just before the section on smoothing
# functions. 

par(mfrow=c(2,1))
matplot(daytime,temp,type='l')
matplot(daytime,precip,type='l')

# We can also plot by region; Atlantic, Pacific, Central and North. 

matplot(daytime,temp,type='l',col=as.factor(CanadianWeather$region))
matplot(daytime,precip,type='l',col=as.factor(CanadianWeather$region))

# But let's get on to fda. 

################################################################################
###                          DEFINING A BASIS SYSTEM                        ####
################################################################################

# We'll always need the range

dayrng = c(0,365)

#### 1. Fourier Basis with 365 basis functions

fbasis = create.fourier.basis(dayrng,365)

# Plot a basis with just 5 components

plot(create.fourier.basis(dayrng,5))

# etc

# Let's try a simple linear regression of the first temperature record on the 
# first, say, 20 basis functions

fb.values = eval.basis(day5,fbasis)

# the 74 by 365 matrix that results has rows as days, columns as bases

dim(fb.values)
plot(day5,fb.values[,1])
plot(day5,fb.values[,2])

# Extract the first temperature record

ex.temp = temp[,1]

# Run a linear regression on temperature with the first 5 basis functions. In
# this case I need to evaluate the basis at the observation times

Xmat = eval.basis(daytime,fbasis)
Xmat = Xmat[,1:5]   # First 5 basis functions

ex.temp.mod = lm(ex.temp~Xmat)

# Let's plot this; the fitted values returned by lm will represent the smooth
# well enough to plot. 

plot(daytime,ex.temp)
lines(daytime,ex.temp.mod$fitted,col=2,lwd=2)

# We can also look at residuals

plot(daytime,ex.temp.mod$resid)

# There's some clear autocorrelation; more basis functions may be warranted. 

###  EXERCISE: repeat the above with different numbers of basis functions (say
###  20, 50, 100, 200, 365. How many look like they give a reasonable smooth?

#### 2. B-spline bases with knots every 5 days

# First of all define a knot sequence; this will be the same as day5

knots = day5

# We'll use fourth-order B-splines

norder = 4

# this implies the number of basis functions

nbasis = length(knots) + norder - 2

# Now we can define the basis

bbasis = create.bspline.basis(dayrng,nbasis,norder,knots)
  
# If in doubt, we can obtain

bbasis$nbasis    # number of basis functions
bbasis$rangeval   # basis range

# Plotting these is a little crazy

plot(bbasis)

# but we can look at a smaller number

plot(create.bspline.basis(dayrng,nbasis=12,norder))

# We can also look at the inner product of these

in.mat = inprod(bbasis,bbasis)

par(mfrow=c(1,1))
image(in.mat)

# and see that it is zero outside of a diagonal band; this can help computation
# a great deal. 

### EXERCISE: try changing the order of the basis and observe how the width
### of the support of the basis changes and how its smoothness properties change. 

### EXERCISE: obtain a least squares smooth of these data with the Bspline basis
### how does this compare with a least squares smooth using a Fourier basis with
### the same number of basis functions?

################################################################################
###                          SMOOTHING FUNCTIONS                            ####
################################################################################

#### 1. Lfd Objects

# Two common means of generating Lfd objects
# 1. int2Lfd -- just penalize some derivatives. 

curv.Lfd = int2Lfd(2)

# 2. vec2Lfd -- a (constant) linear combination of derivatives; for technical
# reasons this also requires the range of the basis. 

harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0),rangeval=dayrng)

# looking inside these objects is not terribly enlightening. 

#### 2. fdPar objects

# We'll concentrate on B-splines and second-derivative penalties.

# First, a value of lambda  (purposefully large so that we can distinguish a fit
# from data below).

lambda = 1e6

# Now we can define the fdPar object 

curv.fdPar = fdPar(bbasis,curv.Lfd,lambda)

#### 3. Smoothing functions

# We're now in a position to smooth

tempSmooth1 = smooth.basis(daytime,temp,curv.fdPar)

# Let's look at the result

names(tempSmooth1)

# First of all, let's plot it

plot(tempSmooth1$fd)

# There is also a neat utility to go through each curve in turn and look at its
# fit to the data:

plotfit.fd(temp,daytime,tempSmooth1$fd)


# Let's examine some fit statistics

# degrees of freedom

tempSmooth1$df

# Just about equivalent to fitting 5 parameters

# We'll also look at GCV, this is given for each observation

tempSmooth1$gcv

# Let's change to a more realistic value of lambda

lambda = 1e1
curv.fdPar$lambda = lambda

tempSmooth = smooth.basis(daytime,temp,curv.fdPar)

# and repeat the previous steps

plotfit.fd(temp,daytime,tempSmooth$fd)
tempSmooth$df
tempSmooth$gcv

# Here the fit looks a lot better and the gcv values are much smaller. 

#### 4. Choosing smoothing parameters

# We can search through a collection of smoothing parameters to try and find
# an optimal parameter.

# We will record the average gcv and choose lambda to be the minimum of these. 

lambdas = 10^seq(-4,4,by=0.5)    # lambdas to look over

mean.gcv = rep(0,length(lambdas)) # store mean gcv


for(ilam in 1:length(lambdas)){
  # Set lambda
  curv.fdPari = curv.fdPar
  curv.fdPari$lambda = lambdas[ilam]

  # Smooth
  tempSmoothi = smooth.basis(daytime,temp,curv.fdPari)

  # Record average gcv
  mean.gcv[ilam] = mean(tempSmoothi$gcv)
}

# We can plot what we have

plot(lambdas,mean.gcv,type='b',log='x')

# Lets select the lowest of these and smooth

best = which.min(mean.gcv)
lambdabest = lambdas[best]

curv.fdPar$lambda = lambdabest
tempSmooth = smooth.basis(daytime,temp,curv.fdPar)

# And look at the same statistics

plotfit.fd(temp,daytime,tempSmooth$fd)
tempSmooth$df

# We'll also plot these

plot(tempSmooth)

### EXERCISE: try obtaining a smooth of the precipitation data

### EXERCISE: how much does the result change if the basis has a knot every day
### instead of every 5 days?


################################################################################
###       FUNCTIONAL DATA OBJECTS: MANIPULATION AND STATISTICS              ####
################################################################################

## Now that we have a functional data object we can manipulate them in various 
# ways.  First let's extract the fd object

tempfd = tempSmooth$fd

# if we look at what's in this we see

names(tempfd)

# We see a basis, plus coefficient matrix

dim(tempfd$coefs)

# and an array giving names

tempfd$fdnames

# With lists giving names for time points, replicates and dimensions. Each list
# also has a name that can be used in plotting.  Apart from plotting functions, 
# fdnames isn't used and you can generally ignore it. 

# We can also create fd objects by adding a basis and a coefficient array. Let's 
# make a random one, say

newcoefs = matrix(rgamma(nbasis*10,5,2),nbasis,10)
newfd = fd(newcoefs,bbasis)

# Notice that we haven't specified fdnames. 

# The plotting command nicely draws these. 

plot(newfd)

# Not that this looks very nice; we'll stick with the Canadian weather data. 


#### 1. Manipulation

# We can do a number of things with these functions, treating them as data. 
# These operations all result in new functional data objects, but we will plot
# them directly as an illustration. 

# Subset

plot(tempfd[1:10])

# We can add them together; the 'lines' function also works with them

newfd = tempfd[1] + tempfd[2]
plot(newfd)
lines(tempfd[1],col=2)
lines(tempfd[2],col=4)

# We can also multiply

plot(tempfd[1]*tempfd[2])

# And obtain powers

plot(tempfd[1]^2)

# We can also obtain derivatives

plot(deriv.fd(tempfd))

# These are pretty wild because of the roughness of the resulting curves
# instead let's have a look at the over-smoothed data:

plot(deriv.fd(tempSmooth1$fd))

# We can also look at second derivatives

plot(deriv.fd(tempSmooth1$fd,2))

# Note that it is a property of B-splines of order m that the (m-2)th derivative 
# is zero at the end of the interval. 

#### 2. Summary statistics

# The obvious thing to look at is the mean

mtempfd = mean(tempfd)
plot(tempfd,col=4)
lines(mtempfd,lwd=2,col=2)

# We can also examine a variance

temp.varbifd = var.fd(tempfd)

# temp.varbifd is a bivariate functional data object -- meaning it takes values
# on a rectangle. 

# To plot this, we need to evaluate it; here we'll use day5 -- 365 points is
# a bit overkill.

temp.var = eval.bifd(day5,day5,temp.varbifd)
contour(day5,day5,temp.var)

# Mostly high variance in the winter, low in summer. Let's have a look at 
# correlation. In this case, evaluation arguments go in with the function call

temp.cor = cor.fd(day5,tempfd)
filled.contour(day5,day5,temp.cor)

# Here we see high correlation between Summer and Winter temperatures, but 
# much less in the spring and fall (although spring to fall correlation is still
# high). 

### EXERCISE: obtain these for the precipitation data and look at the covariance
### and correlation between temperature and precipitation.  

### EXERCISE: try repeating the above with a Fourier basis and the harmonic 
### acceleration penalty. Does this make much difference?