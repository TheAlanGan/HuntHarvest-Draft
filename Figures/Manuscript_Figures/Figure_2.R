####################################################################################################
########## Modeling Hunting and Harvesting Dynamics between Plant-Disperser Pair ###################
## Kevin De Angeli | kevindeangeli@utk.edu
## Eeman Abbasi    | eabbasi@sas.upenn.edu 
## Alan Gan        | a@utk.edu
## Last updated: Nov 4 2019
##==================================================================================================
### Script description: 
# Produces two contour plots using plotly:
# For each individual contour plot we are varying harvest (both Gt and St) and hunting (highHunting multiplier only) levels
# Contour plot x axis: Harvest
# Contour plot y axis: Hunting 

# Things to note:
# We are using Brazil nut and Agouti parameter values 
#===================================================================================================
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

rm(list=ls()) #Removes the data stored in the variables 
library(colorspace)
library(akima)
library(plotly)
library("heatmaply")


###=======================================================================
### Parameters
###=======================================================================

highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest
highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees under HIGH harvest
agoutiGrowth <- 1.1 # Growth rate of Agoutis in logistic model
constHunt <- 0.25# Percentage of agoutis hunted during LOW hunting

# Carrying Capacities for the tree
seedlingCapacity <- 5000 
saplingCapacity <- 500
adultCapacity <- 100
agoutiCapacity <- 5200
ylimit <- 1 # y-axes For plotting
perc <- 0.9 #Proportion of the carrying capacity that is used to initialize the population of trees.

# Initial Population of Tree by stage
seedlingInit <- perc * seedlingCapacity#5000 
saplingInit <- perc * saplingCapacity#500
adultInit <- perc * adultCapacity#100
agoutiInit <- perc * agoutiCapacity#5000

m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).

#The folloing two variables are used in the sigmoid function that quantifies the mutual dependency.  
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

simulation_time <- 500 # Length of simulation in years

##=======The Original 17-Stage Matrix from Zuidema (2001) 
#Note: The original matrix is used as LOW harvesting.
plant_transition_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_transition_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_transition_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_transition_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
##=============high-harvest multiplier
high_harv <- matrix(1, nrow = 17, ncol = 17)
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_t0_high <- plant_transition_mat * high_harv

###=======================================================================
### Functions
###=======================================================================
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
  Nnext <- (R*N*(1-N/(K*(p))) + N) * (1-H)
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
# tree_growthRate : lambda (growth rate) at each time step


simulation <- function(){
  tree_growthRate <- numeric(simulation_time)
  plant_t0 <- matrix(0, nrow = 17)
  plant_t0[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_t0[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_t0[12:17] <- adultInit/6  #Setting initial population of adult trees
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  agouti_vec <- c(agoutiInit) #Population of the agouti/seed disperser
  
  #Simulation begins here:
  for (i in 1:simulation_time)
  {
    pmat <- plant_t0_high
    h_off <- constHunt
    
    NPrev<- sum(plant_t0)
    #Given certain dependency and the current population of trees calculate p.
    p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_t0[12:17]))*.7+0.3 # bounded between 0.3 and 0.7
    #p is then used to calculate the population of agouti for the year i+1.
    agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    #popujlation can't be negative:
    if(agouti_vec[(i+1)]<0){
      agouti_vec[(i+1)]=0
    }
    #set up the matrix that will be used to caclualte the population of trees for the year i+1:
    plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
    plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i)]) 
    #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
    plant_t0 <- matrix( c((plant_animal_mat * pmat) %*% plant_t0)) #This is the new matrix.
    
    #Summing the stages into 3 categories for better plotting.
    plant_t0_sum <- c( sum(plant_t0[1:4]), sum(plant_t0[5:11]), sum(plant_t0[12:17])) 
    plant_all <- cbind(plant_all, plant_t0_sum)
    
    N <- sum(plant_t0)
    tree_growthRate[i] <- log(N/NPrev) #Calculates the stochastic growth during this year.
  }
  
  loglambsim <- mean(tree_growthRate) #output the mean stochastic growth from simulation_time years
  disperser_pop = agouti_vec[length(agouti_vec)] #storing the last value in the animal pop list
  
  return(list("dispPop" = disperser_pop,"growthRate"= exp(loglambsim)))
}

#=============================================================================================================
### Part 3: Running Simulations 
#=============================================================================================================
#In this program, we vary Harvest and Hunting:
harvest_seq<- seq(0,0.4, 0.005)  
hunt_seq<- seq(0,0.8,0.005)

#High hunting with high/low harvest and L-H transition matrix rate 
disp_pop<-matrix(0,length(hunt_seq),length(harvest_seq))
stoch_growthArray<- matrix(0,length(hunt_seq),length(harvest_seq))

#Indexes of the matrix used to store the output from the simulation runs
numRow<-1 
numCol<-1

high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees

#Loop for each possible combination of harvest_seq/hunt_seq and call the Simulation() Function.
for(k in hunt_seq)
{
  numCol<-1
  constHunt<- k  #the offate multiplier (R_hat)
  for(l in harvest_seq)
  {
    high_harv[1,12:17] <- 1-l #Germination Multiplier,inverse of G_t represents the impact of harvest 
    high_harv[cbind(12:17,12:17)] <- 1-l # Multiplier for survival rate of Adult trees, inverse of S_t represents the impact of harvest 
    
    #updates the plant_transition matrix with the new G_t and S_t values under harvest
    plant_t0_high <- plant_transition_mat * high_harv 
    
    output<- simulation()
    disPop <- output$dispPop
    growthRate<- output$growthRate 
    disp_pop[numRow,numCol]<-disPop
    stoch_growthArray[numRow,numCol] <- growthRate 
    
    numCol<-numCol+1
    
  }
  numRow<- numRow+1
  
}

###=======================================================================
### Ploting
###=======================================================================

# These variables are just related to plotly.
f2 <- list(
  family = "Arial, sans-serif",
  size = 24,
  color = toRGB("black"))

f1 <- list(
  family = "Arial, sans-serif",
  size = 24,
  color = toRGB("black"))

abar <- list(
  size=24,
  color = toRGB("black")
)

barFont <- list(
  family = "Arial, sans-serif",
  size = 24,
  color = toRGB("black")
)

batTtile1 <- list(
  text = 'N <sub>disperser',
  font = barFont
)

batTtile2 <- list(
  text = ' Plant \nGrowth \u3bb',
  font = barFont
)




ax <- list(
  tickfont = f2,
  titlefont = f1,
  title = "Hunting",
  #linecolor = toRGB("black"),
  linewidth = 3,
  tickwidth = 3,
  mirror = "ticks"
)
ay <- list(
  tickfont = f2,
  titlefont = f1,
  title = "Harvest",
  #linecolor = toRGB("black"),
  linewidth = 3,
  tickwidth = 3,
  mirror = "ticks"
)


#Colorbar options: https://bl.ocks.org/pstuffa/d5934843ee3a7d2cc8406de64e6e4ea5

# Plotting the plant stochastic growth rate 
p<-plot_ly(
  x = hunt_seq,
  y = harvest_seq, 
  z = t(stoch_growthArray), 
  type = "contour",
  width = 700, height = 600,
  colors = colorRamp(c('darkgreen', "white")))%>% 
  layout(xaxis = ax, yaxis = ay, showlegend = FALSE) %>%
  layout(plot_bgcolor='rgb(254, 247, 234)') %>% 
  layout(paper_bgcolor='transparent') %>%
  #layout(width = 500, height = 500) %>%
  colorbar(title=batTtile2,  x= 1, y=0.52, len =1.09,tickfont = abar)
p


# Plotting the seed disperser population
q<-plot_ly(
  x = hunt_seq,
  y = harvest_seq, 
  z = t(disp_pop), 
  type = "contour",
  width = 700, height = 600,
  colors = colorRamp(c("brown", "white")) )%>%
  layout(xaxis = ax, yaxis = ay, showlegend = FALSE) %>%
  layout(plot_bgcolor='rgb(254, 247, 234)') %>% 
  layout(paper_bgcolor='transparent') %>%
  colorbar(title=batTtile1 ,x= 1, y=0.52, len=1.08, tickfont = abar)
q
