###========================================================================
### Extracting population projection matrices for animals from COMADRE
###========================================================================
# Download the COMADRE RDATA file here: https://www.compadre-db.org/Data/Comadre
# COMPADRE is plant database
load("~/Downloads/COMADRE_v.2.0.1.RData") # the path may be different for you
sort(unique(comadre$metadata$SpeciesAccepted)) # full list of species in COMADRE
grep("Cebus", comadre$metadata$SpeciesAccepted, ignore.case=T)

comadre$metadata[996,] # Metadata on one of the capuchin monkey entries
comadre$mat[[996]] # comadre$mat is a list of matrices with up to 4 for each observation

###========================================================================
### Overview of stochastic growth rate calculation in 'popbio'
###========================================================================

# install.packages('popbio') # only need to run once
library('popbio')

###========================================================================
### An example using a dataset for a plant included with popbio
###========================================================================
?hudsonia # Golden heather - projection matrices for 1985-1988
data(hudsonia) # stored as a list of matrices (akin to a dict but with R's weird referencing)
hudsonia
hudsonia$A85 # projection matrix for only 1985

### Part 1: Inspect the leading eigenvalue for each individual matrix
lambdas <- lapply(hudsonia, lambda) # "List apply": similar in spirit to python's mapply; applies a function over each element in a list or vector
lambdas # Note that 1985 and 1987 are "bad" years

### Part 2: Calculating stochastic growth rate under different probabilities
    # of observing each transition matrix.
    ## 2.A: Equal probabilities
sgr1 <- stoch.growth.rate(hudsonia, prob=rep(0.25, 4)) # equal probabilities of all years
sgr1 # note that $approx returns Tuljapurkar's approximation of log-lambda and $sim gives you a simulated estimate
    # To extract the true lambdas:
exp(sgr1$approx)
exp(sgr1$sim)
    ## 2.B: More weight on bad years
sgr2 <- stoch.growth.rate(hudsonia, prob=c(4/5,1/5,4/5,1/5), verbose=F) # more weight for "bad" years, and tell it to not print messages
sgr2 # really different results! (big discord between $approx and $sim)
exp(sgr2$approx)
exp(sgr2$sim) # Tuljapurkar's approximation can break when the variation
              # between matrices is large (Morris and Doak Chapter 7)
    ## 2.C: More weight on good years
sgr3 <- stoch.growth.rate(hudsonia, prob=c(1/5,4/5,1/5,4/5), verbose=F) # more weight for "bad" years, and tell it to not print messages
exp(sgr3$sim)

### Useful resources:

# Morris and Doak Chapter 7 (Stage-structured populations; Chapter 6 gives an overview of what you've already seen in Caswell); Chapter 9 pages 3-9 gives a very useful overview of how ecologists think about sensitivity analyses with illustrative examples of geese that are hunted in Alaska.

###========================================================================
### Example sensitivity analysis for non-stage structured populations:
### simulation approach
###========================================================================

# For the birds/mammals, we have a simple function that maps N_t to N_{t+1}.

LogisticGrowthHunt <- function(R, N, K, H) {
	Nnext <- R*N*(1-N/K) - H*N + N
	return(Nnext)
}

VertebratePVA <- function(R, N0, K, HuntMort, HuntProb, t.steps, reps) {
    # R: rmax, N0: starting abundance, K: carrying capacity
    # HuntMort: vector of hunting mortalities (assuming two-state regime, high/low)
    # HuntProb: vector of probabilties for observing each mortality level
    # t.steps: Number of time steps
    # reps: Number of replicates for PVA analysis

    nrun <- function(R, N0, K, HuntMort, HuntProb, t.steps) {
        N <- c(N0)
        for (i in 2:t.steps) {
          Ni <- N[i-1]
          if(runif(1) < HuntProb[1]) {
            # We are in regime 1 of hunting mortality
            N[i] <- LogisticGrowthHunt(R, Ni, K, H=HuntMort[1])
          } else {
            # We are in regime 2
            N[i] <- LogisticGrowthHunt(R, Ni, K, H=HuntMort[2])
          }
        }
        return(N)
    }

    PopMat <- replicate(reps,
                        nrun(R,N0,K, HuntMort, HuntProb, t.steps))



    return(list(Nmat=PopMat, # matrix of population sizes over time
                Nend=PopMat[t.steps,] # final population size
    ))
}

PlotCloud <- function(simdata){
    # From Kevin Shoemaker: http://naes.unr.edu/shoemaker/teaching/NRES-470/LECTURE12.html#step_4:_simulate

    # Visualize population trends
  plot(c(1:nrow(simdata)),simdata[,1],col=gray(0.7),type="l",ylim=c(0,max(simdata)),xlab="Years",ylab="Abundance")

  for(r in 2:ncol(simdata)){
    lines(c(1:nrow(simdata)),simdata[,r],col=gray(0.7),type="l")
  }
}

### Running the functions
PVAruns <- VertebratePVA(1.1, 500, 5000, c(0.05,0.75),c(0.75,0.25),100,50)

### Visualize the different runs
PlotCloud(PVAruns$Nmat)

### What about quasi-extinction risk?

quasiextinction <- function(simdata,Nx){
  # Nx is the quasi-extinction threshold
  Nxvec <- apply(simdata,1,function(t)  length(which(t<=Nx)))
  sum(Nxvec>1)/length(Nxvec)
}

quasiextinction(PVAruns$Nmat, 500) # extinction threshold of 500
quasiextinction(PVAruns$Nmat, 1000) # extinction threshold of 1000

###========================================================================
### Alternative approach for animal sensitivity analysis: ceiling model
###========================================================================
### A different approach assumes geometric population growth (density-independent) yet by including "environmental" stochasticity (the hunting model includes this by definition because the driver of variation for growth rate is higher or lower harvest that is unrelated, at present, to any demographic variables such as birth or death rates).

### The Morris and Doak book that I just uploaded includes this information in
Chapters 2 (pages 3-8, 15-19 in particular); 3 (pages 2-8); 4 describes a "ceiling model" that includes density-dependence.
  ### The Kendall paper describes a caveat for the ceiling model approach (not necessary reading, included for additional context)

### If the ceiling model is interesting, then how would one go about calculating the key parameters, specifically $\mu$ and $\sigma^2$? We do not have observed count data (which the approach assumes as the original type of input data). However, perhaps one could calculate different values of $\lambda$ at low population sizes under different harvest intensities...and use the transition probabilities between high and low harvest pressure to estimate $\sigma^2$ from a set of observations

### Function for the density-independent time to extinction model from Morris and Doak Ch 3
extcdf <- function(mu,sig2,d,tmax=50)
{
	# extcdf(mu,sig2,d,tmax) calculates the unconditional extinction time
	# cumulative distribution function from t=1 to t=tmax for mean and
	# variance parameters mu and sigma2 and log distance from the
	# quasi-extinction threshold d; modified from Lande and Orzack,
	# Proc. Nat. Acad. Sci. USA 85: 7418-7421 (1988), equation 11.
	# modified from the Matlab function of the same name in box 3.3 of
	# Morris and Doak, 2002, Quantitative Conservation Biology

   	G = numeric(tmax)
   	x=1:tmax
   	G = pnorm((-d-mu*x)/sqrt(sig2*x)) +
    	  exp(-2*mu*d/sig2)* pnorm((-d+mu*x)/sqrt(sig2*x))
    return(G)
}

# An example
set.seed(100)
lambdas <- sample(c(1.1, 0.75),100,replace=T, prob=c(0.35,0.65))
occ.probs <- c(length(which(lambdas==1.1))/100, length(which(lambdas==0.9))/100)
geom.mean <- lambdas[1]^occ.probs[1] * lambdas[2]^occ.probs[2] # probably there is a more elegant way to program this, but  ¯\_(ツ)_/¯
extcdf(mu=geom.mean,sig2=var(log(lambdas)), d=log(5000/100))
plot(extcdf(mu=geom.mean,sig2=var(log(lambdas)), d=log(5000/100)), xlab="Time (Years into Future for example)",ylab="Quasi-extinction probability")

### Function for ceiling model
tbar_ceiling_N0K <- function(mu=NULL, sigma2=NULL,K=NULL,Nc=NULL,Nx=NULL) {
    # mu: geometric mean of lambdas
    # sigma2: variance of lambdas
    # K: carrying capacity
    # Nc: Starting population size
    # Nx: Quasi-extinction threshold

    ## Box 4.1 in Morris and Doak with some corrections
    c=mu/sigma2 # top of page 9, convenience parameters
    d=log(Nc/Nx)
    k=log(K/Nx)
    Tbar = (exp(2*c*k)*(1-exp(-2*c*d))-2*c*d)/(2*mu*c) # Equation 4.3 in M&D
    return(Tbar)
}

tbar_ceiling_N0K(mu=geom.mean, sigma2=var(log(lambdas)), K=1000, Nc=1000, Nx=100)
