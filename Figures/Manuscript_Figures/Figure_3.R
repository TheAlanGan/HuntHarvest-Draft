####################################################################################################
########## Modeling Hunting and Harvesting Dynamics between Plant-Disperser Pair ###################
## Kevin De Angeli | kevindeangeli@utk.edu
## Eeman Abbasi    | eabbasi@sas.upenn.edu 
## Alan Gan        | a@utk.edu
## Last updated: Nov 7 2019
##==================================================================================================
### Script description: 
#  Generates a line plot # 1: Sustainable plant harvesting and hunting threshold 
#  x-axis = Constant high hunting (highHunting multiplier)
#  y-axis = Brazil nut stochastic growth rate 

#Generates a line plot # 2: Sustainable disperser harvesting and hunting threshold 
#  x-axis = Constant high hunting (highHunting multiplier)
#  y-axis = Disperser population prportion of the disperser carrying capacity

# Things to note:
# All parameter values are specific to Brazil-nut tree and Agouti 
# Everything is set constant, the only variables that are being varied - representative of a harvesting regime:
# G_t and S_t multipliers - harvest 
# R_hat - offtake multiplier  - hunting
# All of these variables are being varied togther using the same sequence term.  
# Also not sure what is the sustainable level (percentage) for a disperser population a proportion of its carrying capacity.#===================================================================================================
###=================================================================================================
### Packages and Libaries required 
###=================================================================================================
#install.packages('heatmaply')
#install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
#install.packages.2('devtools')
#install.packages("akima")
#install.packages('markovchain')
# make sure you have Rtools installed first! if not, then run:
#install.packages('installr'); install.Rtools()

#devtools::install_github("ropensci/plotly") 
#devtools::install_github('talgalili/heatmaply')
library("heatmaply")

###====================================================================================================
### Part 1: Setting Model Specific Parameters
##=====================================================================================================
highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest
highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees under HIGH harvest
agoutiGrowth <- 1.1 # Growth rate of Agoutis in logistic model
constHunt <- 0.25 # Percentage of agoutis hunted 

# Carrying Capacities for the tree
seedlingCapacity <- 5000 
saplingCapacity <- 500
adultCapacity <- 100
agoutiCapacity <- 5200

ylimit <- 1 # For plotting

perc <- 1.0 #Percent of capacity for which the populations are initialized 
seedlingInit <- perc * seedlingCapacity#5000 # Initial Populations
saplingInit <- perc * saplingCapacity#500
adultInit <- perc * adultCapacity#100
agoutiInit <- perc * agoutiCapacity#5000 * 0.3

m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

time_end <- 1000 # Length of simulation in years
maxt <- 600
high_harv <- matrix(1, nrow = 17, ncol = 17)

##=======The Original 17-Stage Matrix from Zuidema (2001) 
plant_transition_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_transition_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_transition_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_transition_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
##=============high-harvest multiplier
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_t0_low <- plant_transition_mat
plant_t0_high <- plant_transition_mat * high_harv
#===========================================================================

#=============================================================================================================
### Part 2: Setting up FUNCIONS
#=============================================================================================================
#2.1 Sigmoid Function 
# k:  steepness 
# x0: midpoint 
# x:  population size

sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} 

#A different functional form that coudl be implemented.
#It's not used in the paper.
linear <- function(m, x, b)
{
  y <- m*x + b
  return(y)
}

#2.3 Discrete-time  disperser logistic growth function
# R: Growth rate (rmax)
# N: Population size
# K: Disperser carrying capacity 
# H: Offtake rate
# p: Plant affect on disperser (from 0-1)
LogisticGrowthHunt<- function(R, N, K, H, p) 
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- (R*N*(1-N/(K*(p)))+N) *(1-H)
  return(Nnext)
} 


###=======================================================================
### Simulation Function
###=======================================================================
#2.4 Stochastic growth of plant 
# delta_d         : dependence of disperser on the plant (δd,p)
# delta_p         : dependence of plant on disperser (δp,d)
# plant_t0[1:4]   : setting initial population of seedlings 
# plant_t0[5:11]  : setting initial population of saplings
# plant_t0[12:17] : setting initial population of adult trees
# plant_all       : summed plant population at all times 
# agouti_vec      : agouti population at each time step 
# prevN           : storing plant population at the previous time step
# plant_animal_mat: disperser affect on the plant 
# growthRate      : lambda (growth rate) at each time step
stoch_growth <- function(){
  growthRate <- numeric(maxt)
  popAfterTime<- numeric(maxt)
  
  plant_t0 <- matrix(0, nrow = 17)
  plant_t0[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_t0[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_t0[12:17] <- adultInit/6  #Setting initial population of adult trees
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  agouti_vec <- c(agoutiInit)
  
  #Simulation begins here:
  for (i in 1:maxt)
  {
    pmat <- plant_t0_high
    h_off <- constHunt
    
    NPrev<- sum(plant_t0)
    #Given certain dependency and the current population of trees calculate p.
    p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_t0[12:17]))*.5+0.5 # bounded between 0.5 and 1.0
    #p is then used to calculate the population of agouti for the year i+1.
    agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    #popujlation can't be negative:
    if (agouti_vec[j+1] < 0)
    {
      agouti_vec[j+1] <- 0
    }
    #set up the matrix that will be used to caclualte the population of trees for the year i+1:
    plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
    plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i)])
    plant_t0 <- matrix( c((plant_animal_mat * pmat) %*% plant_t0))
    
    #Summing the stages into 3 categories for better plotting
    plant_t0_sum <- c( sum(plant_t0[1:4]), sum(plant_t0[5:11]), sum(plant_t0[12:17])) 
    plant_all <- cbind(plant_all, plant_t0_sum)
    
    N <- sum(plant_t0)
    growthRate[i] <- log(N/NPrev)
    
  }
  
  loglambsim <- mean(growthRate)
  disperser_pop = agouti_vec[length(agouti_vec)] #storing the last value in the animal pop list
  
  
  
  return(list("dispPop" = disperser_pop,"growthRate"= exp(loglambsim)))
}

#=============================================================================================================
### Part 3: Running Simulations 
#=============================================================================================================
#High hunting with high/low harvest and L-H transition matrix rate 
xseq<-seq(0,1,0.001)
growthRate_mat<-matrix(0,length(xseq))
animal_pop_mat <- matrix(0, length(xseq))

num<-1 #Used as an index.
constHunt<- 0.25
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_t0_low <- plant_transition_mat
plant_t0_high <- plant_transition_mat * high_harv  # plant matrix affected by high harvest

for(j in xseq){
  constHunt<-j
  high_harv[1,12:17] <- 1-j 
  high_harv[cbind(12:17,12:17)] <- 1-j
  plant_t0_high <- plant_transition_mat * high_harv
  output<- stoch_growth()
  growth_rate <- output$growthRate
  disp_pop<- output$dispPop
  growthRate_mat[num]<-growth_rate
  animal_pop_mat[num]<- (disp_pop/agoutiCapacity)
  num<-num+1
}

#==============================================================================================================================================
# Part 4: Plotting
#==============================================================================================================================================
plot(xseq,growthRate_mat, xlab="Harvesting Regime (Harvest & Hunting)", ylab=expression(paste("Stochastic Growth Rate (",lambda,")")), cex= 0.1,  xlim= c(0,0.6) )
abline(1,0, lty = 2)
abline(v=0.110,lty = 2)
text(x=0.110, y=0.969, "*", cex=3)
dev.copy(png,"Hunting_PlantGrowthRate.png")
dev.off()

plot(xseq,animal_pop_mat, xlab="Harvesting Regime (Harvest & Hunting)", ylab="Disperser population (Proportion of carrying capacity)", cex= 0.1 , xlim= c(0,0.6))
abline(0.5,0, lty = 2)
abline(v=0.125,lty = 2)
text(x=0.125, y=0, "*", cex=3)
dev.copy(png,"Hunting_AnimalGrowthRate.png")
dev.off()
