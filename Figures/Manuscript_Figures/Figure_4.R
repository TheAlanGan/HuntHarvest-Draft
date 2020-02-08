####################################################################################################
########## Modeling Hunting and Harvesting Dynamics between Plant-Disperser Pair ###################
## Kevin De Angeli | kevindeangeli@utk.edu
## Eeman Abbasi    | eabbasi@sas.upenn.edu 
## Alan Gan        | a@utk.edu
## Last updated: Nov 4 2019
##==================================================================================================
### Script description: 
# We are varying deltap,d and delta d,p model, all other model parameters are fixed
#===================================================================================================
###=================================================================================================
### Packages and Libaries required 
###=================================================================================================
#install.packages('heatmaply')
#install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
#install.packages.2('devtools')
#install.packages('markovchain')
#install.packages('installr'); install.Rtools()
#devtools::install_github("ropensci/plotly") 

## Clear all
rm(list=ls(all=TRUE)) 
library(plotly)

###====================================================================================================
### Part 1: Setting Model Specific Parameters
##=====================================================================================================

# Multipliers terms 
highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest (G_t)
highHarvestSurvival <- 0.9   # Multiplier for survival rate of Adult trees under HIGH harvest (S_t)
highHunting <- 0.35          # Percentage of agoutis hunted during HIGH hunting (R^)

# Growth rate of Agoutis in logistic model (Rmax)
agoutiGrowth <- 1.1 

#Carrying Capacities 
seedlingCapacity <- 5000 
saplingCapacity <- 500
adultCapacity <- 100
agoutiCapacity <- 5200
perc <- 0.9   

# Initial Populations 
seedlingInit <- perc * seedlingCapacity 
saplingInit <- perc * saplingCapacity
adultInit <- perc * adultCapacity
agoutiInit <- perc * agoutiCapacity

# m is the desired proportion at which sigmoid(m) = m. Ideally it is small (~0.01-0.05).
m <- 0.05   

#Steepness parameters needed for sigmoid(m) = m. See equation (11)
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) 
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  

# Length of simulation in years
maxt <- 500 

##=======The Original 17-Stage Matrix from Zuidema and high-harvest multiplier
plant_transition_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_transition_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_transition_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_transition_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
high_harv <- matrix(1, nrow = 17, ncol = 17)
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival 

####=========Setting up plant transition matrix under high harvest===========================================
plant_t0_high <- plant_transition_mat * high_harv

#=============================================================================================================
### Part 2: Setting up FUNCIONS
#=============================================================================================================

#2.1 Sigmoid Function 
# k:  steepness 
# x0: midpoint 
# x:  population size
sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) 
} 

#2.3 Discrete-time  disperser logistic growth function
# R: Growth rate (rmax)
# N: Population size
# K: Disperser carrying capacity 
# H: Offtake rate
# p: Plant affect on disperser (from 0-1)
LogisticGrowthHunt<- function(R, N, K, H, p) 
{ 
  Nnext <- (R*N*(1-N/(K*(p))) + N) * (1-H)
  return(Nnext)
}

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
# r               : lambda (growth rate) at each time step
stoch_growth <- function(delta_d, delta_p)
{
  r <- numeric(maxt)
  plant_t0 <- matrix(0, nrow = 17)
  plant_t0[1:4] <- seedlingInit/4   
  plant_t0[5:11] <- saplingInit/7  
  plant_t0[12:17] <- adultInit/6  
  
  plant_all <- matrix( c(sum(plant_t0[1:4]), sum(plant_t0[5:11]), sum(plant_t0[12:17])) ) 
  
  agouti_vec <- c(agoutiCapacity)
  
  for (i in 1:maxt)
  {
    pmat <- plant_t0_high
    h_off <- highHunting
    prevN <- sum(plant_t0) 
    p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_t0[12:17]))*delta_d + (1-delta_d) # bounded between 0.9 and 1.0, k was 0.1
    agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    
    if(agouti_vec[(i+1)]<0){
      agouti_vec[(i+1)]=0
    }
    plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
    plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i)])*delta_p + (1-delta_p) # k was 0.0025
    plant_t0 <- matrix( c((plant_animal_mat * pmat) %*% plant_t0))
    
    #Summing the stages into 3 categories for better plotting
    plant_t0_sum <- c( sum(plant_t0[1:4]), sum(plant_t0[5:11]), sum(plant_t0[12:17])) 
    plant_all <- cbind(plant_all, plant_t0_sum)
    
    N <- sum(plant_t0)
    r[i] <- log(N / prevN)
    
  }
  
  loglambsim <- mean(r)
  disperser_pop = agouti_vec[length(agouti_vec)] #storing the last value in the animal pop list
  
  return(list("dispPop" = disperser_pop,"growthRate"= exp(loglambsim)))
}

#=============================================================================================================
### Part 3: Running Simulations 
#=============================================================================================================
delta_p_seq <- seq(0.1, 0.95, length.out=20)
delta_d_seq <- seq(0.1, 0.95, length.out=20)

# initializing arrays to store plant pop and disp pop
stoch_growthArray<- c()
disPop_Array<- c()

#Creating all possible combinations 
multiplier_vals <- expand.grid(delta_d_seq, delta_p_seq) 
set.seed(8)

for(i in 1:nrow(multiplier_vals)){
  
  delta_d_seq <- multiplier_vals[i,1]
  delta_p_seq <- multiplier_vals[i,2]
  output <- stoch_growth(delta_d_seq,delta_p_seq)
  stoch_growthArray[i]<- output$growthRate
  disPop_Array[i] = output$dispPop
}
#==============================================================================================================================================
# Part 4: Plotting Figure 4
#==============================================================================================================================================
delta_d_seq = multiplier_vals[,1]
delta_p_seq = multiplier_vals[,2]

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
  text = 'Plant \nGrowth (\u3bb)',
  font = barFont
)

batTtile2 <- list(
  text ='N<sub>disp',
  font = barFont
)

ax <- list(
  tickfont = f2,
  titlefont = f1,
  title = paste('\u3b4<sub>p→d'),
  linecolor = toRGB("black"),
  linewidth = 3,
  tickwidth = 3,
  mirror = "ticks")

ay <- list(
  tickfont = f2,
  titlefont = f1,
  title = paste('\u3b4<sub>d→p'),
  linecolor = toRGB("black"),
  linewidth = 3,
  tickwidth = 3,
  mirror = "ticks")

p<-plot_ly(
  x = delta_d_seq,
  y = delta_p_seq, 
  z = stoch_growthArray, 
  type = "contour", 
  width = 700, height = 600,
  colors = colorRamp(c('darkgreen', "white")))%>% 
  layout(xaxis = ax, yaxis = ay) %>%
  layout(plot_bgcolor='rgb(0, 75, 0)') %>% 
  layout(paper_bgcolor='transparent') %>%
  colorbar(title =batTtile1, x= 1, y=1.06, len=1.085,tickfont = abar)
p

q<-plot_ly(
  x = delta_d_seq,
  y = delta_p_seq, 
  z = disPop_Array, 
  type = "contour",
  width = 700, height = 600,
  colors = colorRamp(c("brown", "white")) )%>%
  layout(xaxis = ax, yaxis = ay) %>%
  layout(plot_bgcolor='rgb(254, 247, 234)') %>% 
  layout(paper_bgcolor='transparent') %>%
  colorbar(title =batTtile2,x= 1, y=1.035, len=1.059, tickfont = abar)
q