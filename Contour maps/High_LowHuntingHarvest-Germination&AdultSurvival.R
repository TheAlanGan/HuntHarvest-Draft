###====================================================================================================================================
##Script Description: 
#  Generates a contour map: Color gradient represents Brazil nut stochastic growth rate 
#
#  We vary:
#  Constant high/low hunting (highHunting/lowHunting multiplier) with H/L harvest 
#  Adult Survival (highHarvestSurvival multiplier)
#  
##========================================================================================================



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

high_harv <- matrix(1, nrow = 17, ncol = 17)

xseq<-seq(0,1,0.05)
gseq<- seq(log(1), log(1000), 0.33)


##=======The Original 17-Stage Matrix from Zuidema and high-harvest multiplier
plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_mat_low <- plant_S_mat
plant_mat_high <- plant_S_mat * high_harv
#==========================================================================================================================================
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
  Nnext <- (R*N*(1-N/(K*(p)))+N) *(1-H)
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

stoch_growth <- function(){
  r <- numeric(maxt)
  plant_mat <- matrix(0, nrow = 17)
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  
  agouti_vec <- c(agoutiInit)
  
  harvest_seq<-markovChain()
  
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
    
    NPrev<-sum(plant_mat)
    p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
    agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
    plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i)]) # k was 0.0025
    #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
    plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
    
    #Summing the stages into 3 categories for better plotting
    plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
    plant_all <- cbind(plant_all, plant_mat_sum)
    
    N <- sum(plant_mat)
    r[i] <- log(N/NPrev)
    
  }
  
  loglambsim <- mean(r)
  
  return(loglambsim)
}
#==============================================================================================================================
#Changing both high and low hunting with respect to germination
low_high_huntseq<- seq(0,0.85,0.05)
growthRate_mat<-matrix(0,length(low_high_huntseq),length(low_high_huntseq))
binary_mat<- matrix(0,length(low_high_huntseq),length(low_high_huntseq))
rownames(growthRate_mat) <- paste(low_high_huntseq)
colnames(growthRate_mat) <- paste(low_high_huntseq)

rownames(binary_mat) <- paste(low_high_huntseq)
colnames(binary_mat) <- paste(low_high_huntseq)

highHunting<-0
lowHunting<-0

num<-1 
num1<-1

high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees

for(i in low_high_huntseq)
{
  high_harv[1,12:17] <- i # Multiplier for fecundity rate for Adult trees
  num1<-1
  for(j in low_high_huntseq){
    
    plant_mat_low <- plant_S_mat
    plant_mat_high <- plant_S_mat * high_harv
    lowHunting<-j
    highHunting<- j+0.15
    growth_rate <- exp(stoch_growth())
    growthRate_mat[num,num1]<-growth_rate
    print(growthRate_mat[num,num1])
    if(growthRate_mat[num,num1]>=1)
    {
      binary_mat[num,num1]<-1
    }
    else{
      binary_mat[num,num1]<-0
    }
    
    num1<-num1+1
  }
  num<- num+1
  
}

library(akima)
filled.contour(x = low_high_huntseq, 
               y = low_high_huntseq,
               z = growthRate_mat,
               color.palette = colorRampPalette(c("red", "blue")),
               xlab = "Germination",
               ylab = "H/L Hunting with H/L Harvest",
               key.title = title(main = "Growth Rate", cex.main = 0.7))


#==============================================================================================================================
low_high_huntseq<- seq(0,0.85,0.05)
growthRate_mat<-matrix(0,length(low_high_huntseq),length(low_high_huntseq))
binary_mat<- matrix(0,length(low_high_huntseq),length(low_high_huntseq))
rownames(growthRate_mat) <- paste(low_high_huntseq)
colnames(growthRate_mat) <- paste(low_high_huntseq)

rownames(binary_mat) <- paste(low_high_huntseq)
colnames(binary_mat) <- paste(low_high_huntseq)

highHunting<-0
lowHunting<-0
num<-1 
num1<-1

high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees

for(i in low_high_huntseq)
{
  high_harv[cbind(12:17,12:17)] <- i # Multiplier for survival rate of Adult trees
  num1<-1
  for(j in low_high_huntseq){
    plant_mat_low <- plant_S_mat
    plant_mat_high <- plant_S_mat * high_harv
    lowHunting<-j
    highHunting<- j+0.15
    growth_rate <- exp(stoch_growth())
    growthRate_mat[num,num1]<-growth_rate
    print(growthRate_mat[num,num1])
    if(growthRate_mat[num,num1]<1)
    {
      binary_mat[num,num1]<-0
    }
    else{
      binary_mat[num,num1]<-1
    }
    
    num1<-num1+1
  }
  num<- num+1
  
}

library(akima)
filled.contour(x = low_high_huntseq, 
               y = low_high_huntseq,
               z = growthRate_mat,
               color.palette = colorRampPalette(c("red", "blue")),
               xlab = "Adult Survival",
               ylab = " H/L Hunting with H/L Harvest",
               key.title = title(main = "Growth Rate", cex.main = 0.7))


