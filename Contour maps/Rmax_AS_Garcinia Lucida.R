###====================================================================================================================================
### Parameters
###====================================================================================================================================
# install.packages('heatmaply')
# install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
# install.packages.2('devtools')
# install.packages("akima")
# install.packages('markovchain')
# # make sure you have Rtools installed first! if not, then run:
# #install.packages('installr'); install.Rtools()
# 
# devtools::install_github("ropensci/plotly") 
# devtools::install_github('talgalili/heatmaply')
# library("heatmaply")


highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest
highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees under HIGH harvest
agoutiGrowth <- 1.1 # Growth rate of Agoutis in logistic model

lowHunting <- 0.1 # Percentage of agoutis hunted during LOW hunting
highHunting <- 0.25 # Percentage of agoutis hunted during HIGH hunting

seedlingCapacity <- 5000 # Carrying Capacities. Tree capacities don't matter as much.
saplingCapacity <- 500
adultCapacity <- 100
agoutiCapacity <- 5200
perc <- 0.9
ylimit <- 1 # For plotting

seedlingInit <- perc * seedlingCapacity#5000 # Initial Populations
saplingInit <- perc * saplingCapacity#500
adultInit <- perc * adultCapacity#100
agoutiInit <- perc * agoutiCapacity#5000

m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

time_end <- 1000 # Length of simulation in years

maxt <- 1000
brazilNut <- list(low=plant_mat_low, high=plant_mat_high)


xseq<-seq(0,1,0.25)
low_high_huntseq<- seq(0,0.85,0.05)
gseq<- seq(0.25,1,0.25)
rmaxseq<- seq(0.2,3, 0.25)


##=======The Original 17-Stage Matrix from Zuidema and high-harvest multiplier
plant_S_mat <- matrix( 0, nrow = 8, ncol = 8)
diag(plant_S_mat) <- c(0.90471,0.95385,0.95982,0.95792,00.96938,0.48667,0,0)
plant_S_mat[cbind(2:8,1:7)] <- matrix(c(0.02109,0.03774,0.03361,0.03398,0.01857,0.00067,0))
plant_S_mat[1,5:8] <- matrix( c(13.9346,20.9825,32.3321,33.4871) )
high_harv <- matrix(1, nrow = 8, ncol = 8)
high_harv[1,5:8] <- highHarvestFecundity 
high_harv[cbind(5:8,5:8)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_mat_low <- plant_S_mat
plant_mat_high <- plant_S_mat * high_harv
#===========================================================================

#============FUNCTIONS======================================================

sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} 

linear <- function(m, x, b)
{
  y <- m*x + b
  return(y)
}

LogisticGrowthHunt<- function(R, N, K, H, p) 
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- R*N*(1-N/(K*(p))) - H*N + N
  return(Nnext)
} 

LogisticGrowthHuntRK<- function(R, N, K, H, p,s,m) 
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- R*(s)*N*(1-N/((K*(p))*m)) - H*N + N
  return(Nnext)
} 
# Specifying the markov chain
library('markovchain')

markovChain<- function(){
  
  statesNames = c("low","high")
  mcHarvest <- new("markovchain", states = statesNames, 
                   transitionMatrix = matrix(data = c(0.2, 0.8, 0.8, 0.2), byrow = TRUE, 
                                             nrow = 2, dimnames=list(statesNames,statesNames)), name="Harvest")
  # Simulating a discrete time process for harvest
  set.seed(100)
  harvest_seq <- markovchain::rmarkovchain(n=time_end, object = mcHarvest, t0="low")
  return(harvest_seq)
}
harvest_seq <- markovChain()

stoch_growth <- function(){
  r <- numeric(maxt)
  plant_mat <- matrix(0, nrow = 8)
  plant_mat[1] <- seedlingInit #Setting initial population of seedlings
  plant_mat[2:4] <- saplingInit/3  #Setting initial population of saplings
  plant_mat[5:8] <- adultInit/4  #Setting initial population of adult trees
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  agouti_vec <- c(agoutiInit)
  
  markovChain()
  
  for (i in 1:maxt)
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
    Nprev<- sum(plant_mat)
    p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[5:8]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
    agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    plant_animal_mat <- matrix(1, nrow = 8, ncol = 8)
    if(agouti_vec[i+1]<0)
    {
      agouti_vec[i+1]=0
    }
    plant_animal_mat[1,5:8] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
    #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
    plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
    
    #Summing the stages into 3 categories for better plotting
    plant_mat_sum <- c( sum(plant_mat[1]), sum(plant_mat[2:4]), sum(plant_mat[5:8])) 
    plant_all <- cbind(plant_all, plant_mat_sum)
    
    N <- sum(plant_mat)
    r[i] <- log(N/Nprev)
    
  }
  
  loglambsim <- mean(r)
  
  return(loglambsim)
}
#==========================================================================================================================================================
growth.array<-array(0,dim=c(length(xseq),length(rmaxseq)))
binary_growth.array<-array(0,dim=c(length(xseq),length(rmaxseq)))

agoutiGrowth<-0
num<-1 
num1<-1
high_harv[cbind(5:8,5:8)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees

for(i in xseq)
{
  high_harv[cbind(5:8,5:8)]<-i
  num1<-1
  for(j in rmaxseq){
    
    plant_mat_low <- plant_S_mat
    plant_mat_high <- plant_S_mat * high_harv
    agoutiGrowth<- j     
    growth_rate <- exp(stoch_growth())
    growth.array[num,num1]<-growth_rate
    print(growth.array[num,num1])
    if(growth.array[num,num1]>=1)
    {
      binary_growth.array[num,num1]= 1
    }
    else
    {
      binary_growth.array[num,num1] = 0
    }
    
    
    num1<-num1+1
  }
  
  num<- num+1
  
}

library(akima)
filled.contour(x = xseq,
               y = rmaxseq,
               z = growth.array[,],
               color.palette = colorRampPalette(c("white", "blue")),
               plot.tile= title(xlab = "Adult Survival",
                                ylab = "Rmax", main= "Garcina Lucida"),
               key.title = title(main = "Growth Rate", cex.main = 0.7), nlevels = 20)


#================================================Binary Matrix===========================================================
library(akima)
filled.contour(x = xseq,
               y = rmaxseq,
               z = binary_growth.array[, ],
               color.palette = colorRampPalette(c("red", "blue")),
               plot.tile= title(xlab = "Adult Survival",
                                ylab = "Rmax", main= "Garcina Lucida"),
               key.title = title(main = "Growth Rate", cex.main = 0.7))


#=========================================================================================================================


