library(htmlTable)

###=======================================================================
### Parameters
###=======================================================================
highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees
highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees
lowHunting <- 0.1 # Percentage of agoutis hunted during LOW hunting
highHunting <- 0.25 # Percentage of agoutis hunted during HIGH hunting
agoutiGrowth <- 0.5 # Growth rate of Agoutis in logistic model

seedling_init <- 5000 #Initial Populations
sapling_init <- 500
adult_init <- 100
agouti_init <- 5200

time_end <- 100 # Length of simulation in years

#=========================================================================


###========================================================================
### The Original 17-Stage Matrix from Zuidema and high-harvest multiplier
###========================================================================
plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )

high_harv <- matrix(1, nrow = 17, ncol = 17)
high_harv[1,12:17] <- highHarvestFecundity    # Multiplier for fecundity rate for Adult trees
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
#===========================================================================



###====================================================================
### Plant and animal dynamics under harvest
###====================================================================
# Plant
plant_mat_low <- plant_S_mat
plant_mat_high <- plant_S_mat * high_harv
sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} # plant_animal <- matrix(c(1,1,sigmoid(1, K_animal*0.5, N),1,1,1,1,1,1), nrow=3, byrow=T)
# Animal

LogisticGrowthHunt <- function(R, N, K, H, p) 
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- R*N*(1-N/(K*(p))) - H*N + N
  return(Nnext)
} # some proportion of N_t are harvested
###=====================================================================


library(markovchain)
###====================================================================
### Harvest markov chain
###====================================================================
# Specifying the markov chain
statesNames = c("low","high")
mcHarvest <- new("markovchain", states = statesNames, 
                 transitionMatrix = matrix(data = c(0.2, 0.8, 0.8, 0.2), byrow = TRUE, 
                                           nrow = 2, dimnames=list(statesNames,statesNames)), name="Harvest")
# Simulating a discrete time process for harvest
set.seed(100)
harvest_seq <- rmarkovchain(n=time_end, object = mcHarvest, t0="low")
head(harvest_seq)
###====================================================================



###====================================================================
### For loop simulation
###====================================================================
plant_mat <- matrix(0, nrow = 17)

plant_mat[1:4] <- seedling_init/4   #Setting initial population of seedlings
plant_mat[5:11] <- sapling_init/7   #Setting initial population of saplings
plant_mat[12:17] <- adult_init/6  #Setting initial population of adult trees
agouti_vec <- c(agouti_init) # Initializing the vector containing agouti pop at each timestep

plant_all <- matrix( c(seedling_init, sapling_init, adult_init) ) # This will contain the summed plant populations at ALL timesteps


for (i in 1:time_end) 
{
  h_i <- harvest_seq[i]
  
  if (h_i == "low") 
  {
    pmat <- plant_mat_low
    h_off <- lowHunting
  } 
  
  else 
  {

    pmat <- plant_mat_high
    h_off <- highHunting

  }
  
  p <- sigmoid(0.1, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0
  agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],5200,h_off, p)
  plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
  plant_animal_mat[1,12:17] <- sigmoid(0.0025, 2600, agouti_vec[(i+1)])
  plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
  
  #Summing the stages into 3 categories for better plotting
  plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
  plant_all <- cbind(plant_all, plant_mat_sum)
  
#  print(sum(plant_mat[12:17]))
  
}

par(mar=c(5,4,1,1),oma=c(0,0,0,0))
plot(agouti_vec/5200,xlab="",ylab="Population size/Max size", col="brown", ylim=c(0,1.5), type="l",xlim=c(1,time_end),xaxs="i")
lines(plant_all[1,]/5000, col="forestgreen")
lines(plant_all[2,]/500, col="turquoise3")
lines(plant_all[3,]/100, col="orange")
legend("bottomleft", c("Animal","Adult","Sapling","Seedling"),col=c("brown","orange","turquoise3","forestgreen"),lty=c(1,1,1,1), bty="n",ncol=2)
mtext("Time step", 1, line=1.85, at=25, col="black")
axis(1,1:time_end,labels=toupper( substr(harvest_seq,1,1) ),line=2,col=NA,col.ticks=NA,col.axis="black", cex.axis=0.65)
mtext("Harvest:",1,line=3,at=-2.5,col="black")