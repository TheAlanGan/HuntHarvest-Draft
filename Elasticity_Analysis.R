#Elasticity Analysis

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
#agoutiInit <- 5000


m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

time_end <- 100 # Length of simulation in years
maxt <- 100 #controls the markovChain

high_harv <- matrix(1, nrow = 17, ncol = 17)

xseq<-seq(0,1,0.05)
low_high_huntseq<- seq(0,0.85,0.05)

##The Original 17-Stage Matrix from Zuidema and high-harvest multiplier---------------------------------
plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_mat_low <- plant_S_mat
plant_mat_high <- plant_S_mat * high_harv


brazilNut <- list(low=plant_mat_low, high=plant_mat_high)

#--------------------------------------------------------------------------------------------------------

#===========================================================================
#============FUNCTIONS======================================================
#===========================================================================

#----------------------------Sigmoid----------------------------
sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} 

#----------------------------Linear----------------------------

linear <- function(m, x, b)
{
  y <- m*x + b
  return(y)
}
#----------------------------LogisticGrowth----------------------------

LogisticGrowthHunt<- function(R, N, K, H, p) 
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- R*N*(1-N/(K*(p))) - H*N + N
  return(Nnext)
} 
#----------------------------LogisticGrowthHunt----------------------------

LogisticGrowthHuntRK<- function(R, N, K, H, p,s,m) 
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- R*(s)*N*(1-N/((K*(p))*m)) - H*N + N
  return(Nnext)
} 
# Specifying the markov chain
#install.packages('markovchain')

#----------------------------MarkovChain----------------------------
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
###===========================================================================
### Sensitivity Matrix
###===========================================================================
sensitivity_matrix <- function(A)
{
  # Computes the Sensitivity and Elasticity Matrices by method specified by Caswell 2001
  #   Chapter 9.
  #
  # Inputs:
  #   A : (n,n) matrix
  #
  # Returns:
  #   A list with 3 items:
  #     sensMat : (n,n) matrix
  #                 i,j values are sensitivity values for each i,j element of A
  #
  #     elasMat : (n,n) matrix
  #                 i,j values are elasticity values for each i,j element of A
  #
  #     stableAge : (n) vector
  #                 The stable Age distribution for A
  
  eigvals <- eigen(A)$values
  eigvecs <- eigen(A)$vectors
  W <- eigvecs
  
  sensMat <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  diag(sensMat) <- c(eigvals)
  
  dominant <- 0
  for (i in eigvals)
  {
    dominant <- max(sqrt(Re(i*Conj(i))), dominant)
  }
  finalCount <- which(eigvals==dominant) # The index of dominant eigenvalue.
  
  V <- Conj(solve(W)) # solve() finds inverse matrix
  w <- matrix(W[,finalCount])
  v <- matrix(Re(V[finalCount,]))
  
  sensMat <- Re(v %*% t(w))
  
  elasMat <- (sensMat * A) / dominant
  
  stableAgeDistr <- ((matrix(Re(eigvecs[,finalCount]) * -1)) / norm(matrix(Re(eigvecs[,finalCount]) * -1))) # Stable age distribution
  
  return(list(sens=sensMat, elas=elasMat, stableAge=stableAgeDistr))
}


#-----------------------------------------------------------------------------
#Bar plot - Elasticity matrix
#-----------------------------------------------------------------------------
#plant_mat_elsaticity <-q
elasticity_demographic_process_bar_plot <- function(plant_mat_elsaticity)
{
elasticity_matrix <- sensitivity_matrix(plant_mat_elsaticity)$elas
growth_elas <- sum(diag(elasticity_matrix[2:17,1:16])) #Elasticity for growth
surv_elas <- sum(diag(elasticity_matrix)) #Elasticity for survival
germination_elas <- sum(elasticity_matrix[1,12:17]) #Elasticity for germination
growth_elas + surv_elas + germination_elas #Must equal 1
df <- data.frame(Demographic_Process=c("Growth", "Survival", "germination_elas"),
                 Elasticity=c(growth_elas, surv_elas, germination_elas))
#head(df)
library(ggplot2)

ggplot(data=df, aes(x=Demographic_Process, y=Elasticity)) +
  geom_bar(stat="identity")

}

#-----------------------------------------------------------------------------
#Contineous graph - Elasticity matrixes. 
#-----------------------------------------------------------------------------
#Pass a vector of values for each demographic process and it will plot it.
#pass a 0 if you want germination to NOT be displayed. 

#It was built to pass the demographic processes of one specific stage 
# (i.e. Adult: Survival, Germination, Growth)
elasticity_contineous_graph <- function(survival,growth,germination)
{
  #Graph Elasticity by stage  ---Sapling
  par(mar=c(5,4,1,1),oma=c(0,0,0,0))
  plot(survival,xlab="",ylab="Elasticity", col="brown", ylim=c(0,.05), type="l",xlim=c(1,time_end),xaxs="i")
  lines(growth, col="orange")
  lines(germination, col="green")
  legend("bottomleft", c("Survival","Growth","Germination"),col=c("brown","orange","green"),lty=c(1,1,1,1), bty="n",ncol=2)
  mtext("Time step", 1, line=1.85, at=25, col="black")
  #axis(1,1:time_end,labels=toupper( substr(harvest_seq,1,1) ),line=2,col=NA,col.ticks=NA,col.axis="black", cex.axis=0.65)
  mtext("Harvest:",1,line=3,at=-2.5,col="black")
}



###===========================================================================
### Running simulation here:
###===========================================================================

#Variables used to buildt the elasticity vectors and create a linear graph.
growth_elas_vec <- c()
surv_elas_vec <- c()
germination_elas_vec <- c()

avg_survival_seedling_elas <- c()
avg_growth_seedling_elas <- c()
avg_survival_sapling_elas <- c()
avg_growth_sapling_elas <- c()
avg_survival_adult_elas <- c()
avg_growth_adult_elas <- c()
avg_germination_adult_elas <-c()

avg_growth_elas <- c()
avg_germ_elas <- c()
#__________________________________________________________________________

harvest_seq <- markovChain()
plant_mat <- matrix(0, nrow = 17)

plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep

plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps

plot_var <- 1 #Controls how many plots will be produced
plot_var2 <- time_end/3 #we can add 1000/2 to plot_var two times to get one plot in the middle and one at the end.
                        #This variable is used below in an "if" statement. 

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
  
  p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
  agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    
  
  plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
  plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
  #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
  plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
  
  #Summing the stages into 3 categories for better plotting
  plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
  plant_all <- cbind(plant_all, plant_mat_sum)
  
#-----------------------------------------------------------------------------
#Elasticity stages
#-----------------------------------------------------------------------------  
  # avg_survival_seedling_elas <- c()
  # avg_growth_seedling_elas <- c()
  # avg_survival_sapling_elas <- c()
  # avg_growth_sapling_elas <- c()
  # avg_survival_adult_elas <- c()
  # avg_growth_adult_elas <- c()
  # avg_germination_adult_elas <-c()
  
dynamic_transition_matrix <- plant_animal_mat * pmat
elasticity_matrix <- sensitivity_matrix(dynamic_transition_matrix)$elas 
#elasticity_matrix <-  plant_S_mat Just used to check if the values that come next match. 
#Seedilging
avg_survival_seedling_elas[i] <- sum(diag(elasticity_matrix[1:4,1:4]))/4  
avg_growth_seedling_elas[i] <- sum(diag(elasticity_matrix[2:5,1:4]))/4 
#Sapling
avg_survival_sapling_elas[i] <- sum(diag(elasticity_matrix[5:11,5:11]))/7
avg_growth_sapling_elas[i] <- sum(diag(elasticity_matrix[6:12,5:11]))/7
#Adult
 avg_survival_adult_elas[i] <- sum(diag(elasticity_matrix[12:17,12:17]))/6
 avg_growth_adult_elas[i] <- sum(diag(elasticity_matrix[13:17,12:16]))/5
 avg_germination_adult_elas[i] <- sum(diag(elasticity_matrix[1,12:17]))/6

#-----------------------------------------------------------------------------
#Elasticity Demographic Proceses 
#-----------------------------------------------------------------------------
#This code takes a matrix as an argument and it descomposes it into three vectors.
#one for each of the three demographic proceses (Surviving,Growing,Germinating)
# dynamic_transition_matrix <- plant_animal_mat * pmat
# elasticity_matrix <- sensitivity_matrix(dynamic_transition_matrix)$elas
# 
growth_elas_vec[i] <- sum(diag(elasticity_matrix[2:17,1:16])) #Elasticity for growth
surv_elas_vec[i] <- sum(diag(elasticity_matrix)) #Elasticity for survival
germination_elas_vec[i] <-sum(elasticity_matrix[1,12:17]) #Elasticity for germination
  

}

##ggplot. 
#Hypercube latinHypercube. 
#elasticity_contineous_graph(avg_survival_adult_elas,avg_growth_adult_elas,avg_germination_adult_elas)

#dynamic_transition_matrix <- plant_animal_mat * pmat
#elasticity_demographic_process_bar_plot(dynamic_transition_matrix)


















