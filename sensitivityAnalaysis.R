###=======================================================================
### Parameters
###=======================================================================
#install.packages('heatmaply')
#install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
#install.packages.2('devtools')
# make sure you have Rtools installed first! if not, then run:
#install.packages('installr'); install.Rtools()

#devtools::install_github("ropensci/plotly") 
#devtools::install_github('talgalili/heatmaply')
#devtools::install_github('spedygiorgio/markovchain')
library("ggplot2")
#library("heatmaply")
#remove.packages("ggplot2")
#install.packages("ggplot2")

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

maxt <- 1000 #controls the markovChain

brazilNut <- list(low=plant_mat_low, high=plant_mat_high)
high_harv <- matrix(1, nrow = 17, ncol = 17)


##=======The Original 17-Stage Matrix from Zuidema and high-harvest multiplier
plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_mat_low <- plant_S_mat
plant_mat_high <- plant_S_mat * high_harv


#===========================================================================
### Additional parameter values
#===========================================================================

### We will use Lamont Cole (1954)'s equation to calculate different ranges of rMax for species

amniote <- read.csv("~/Google Drive/Data/Indicator/Amniote/Amniote_Database_Aug_2015.csv", header=T, stringsAsFactors = FALSE) # global data base on traits for many mammals, birds, and some reptiles
amniote$Latin <- trimws(paste(amniote$genus, amniote$species)) # have to create a Latin name field
rmNAs <- function(x) {ifelse(x==-999, NA, x)} # many databases code missing values as -999 which is annoying
cole <- function(r) { exp(-r) + b*exp(-r*a) - b*exp(-r*(w+1)) -1} # Cole (1954) equation for r-max
find.mass <- function(x) {
  mass.vec <- rmNAs(amniote$adult_body_mass_g[match(x, amniote$Latin)])
  names(mass.vec) <- x
  return(mass.vec)
}

example_species <- c("Tapirus terrestris","Cuniculus paca","Cebus apella","Buceros rhinoceros","Ramphastos toco")

# 1: Indices for extracting values from amniote
inds <- match(example_species, amniote$Latin)

# 2A: Getting values for Cole 1954: Age at first reproduction
  # Amniote AFR is in days (convert to months)
  # Convert -999 to NAs at the same time
a.cole <- rmNAs(amniote$female_maturity_d[inds])/30     # AFR from Amniote

# 2B: Cole 1954: Litter size
b.cole <- rmNAs(amniote$litter_or_clutch_size_n[inds])

# 2C: Cole 1954: age at last reproduction (w) (interpreted as longevity)
w.cole <- rmNAs(amniote$longevity_y[inds])*12

# 3: Body masses
BMs <- find.mass(example_species)

## 4: Using "multiroot" (equivalent to python fsolve) to find Cole's r_max 
rcole <- c(); r.spp <- c(); cmass <- c()
# NB: a and w have to be in years^-1; b = b/2 for females only
for (i in 1:length(a.cole)) {
  a<-a.cole[i]/12; b<-b.cole[i]/2; w<-w.cole[i]/12
  if (!is.na(a) & !is.na(b) & !is.na(w)) {
    rcole <- append(rcole, rootSolve::multiroot(cole,0.5, maxiter=10000, atol=10^-16)$root, after=length(rcole))
    r.spp <- append(r.spp, example_species[i], after=length(r.spp))
    cmass <- append(cmass, BMs[i], after=length(cmass))
  }
}

rmaxdf <- data.frame(Species=r.spp,
                     rMax = rcole,
                     BM_kg = as.numeric(cmass)/1000,
                     stringsAsFactors = FALSE)
rmaxdf

### Another way: allometric scaling
  # Henneman 1983: 4.9*W^-0.2622 # rmax (year-1); rmax (day-1) => 0.0134 * W^-0.2622
  # Hamilton et al 2011, Proc B, Figure 6: rmax (year^-1) = 10^0.63 (4.3) * W (grams)^-1/3 => rmax (day-1) = 0.0117   W^(-1/3)     # probably least reliable - weird marsupial/placental divide
  # Fenchel 1974, end of results right before Discussion, rmax (day^-1) = 10^(-1.4)*X^(-0.275) [0.0398 W ^-0.275 for  day-1]
  # Blueweiss: rmax (day^-1) = 0.025*W^(-0.26) # W is in grams
  # Sibly and West, 2007, PNAS: Results - all mammal production rate: rmax (yr-1): 0.98*M (grams)^-0.275 => 0.00268 W^-0.275 (day-1)
rm.pow <- -0.25 # intercept and exponent for rmax - body mass relationship
rmax_allo <- 0.98*(BMs^rm.pow) # not going to be accurate for birds; need a different intercept
BMlogs <- log10(BMs/1000)

plot(BMlogs, rmax_allo, pch=19, xlab=expression(paste("Body mass (", log[10]," kgs)")), ylab=expression(r[max]), xlim=c(-0.5,2.5),ylim=c(0,0.5), xaxs="i", yaxs="i")
points(log10(rmaxdf$BM_kg),rmaxdf$rMax, col="red", pch=19)
legend("topleft",c("Cole","Allometric"),pch=rep(19,2), col=c("red","black"),bty="n")
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


#----------------------------agouti_Abundance----------------------------

agouti_Abundance<- function(s,m){
  
  plant_mat <- matrix(0, nrow = 17)
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  plant_mat <- plant_mat / sum(plant_mat)
  
  agouti_vec <- c(agoutiInit)
  
  markovChain()
  
  for (i in 1:maxt)
  {
    h_i <- harvest_seq[i]
    
    if (h_i=="low") 
    {
      h_off <- lowHunting
    } 
    
    else 
    {
      h_off <- highHunting
    }
    
    p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
    agouti_vec[(i+1)] <- LogisticGrowthHuntRK(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p,s,m)
    
  }
  return(agouti_vec[length(agouti_vec)])
}

#----------------------------plant_abundance_underHighHunt----------------------------

plant_abundance_underHighHunt<- function(highHunting){
  plant_mat <- matrix(0, nrow = 17)
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  plant_mat <- plant_mat / sum(plant_mat)
  
  agouti_vec <- c(agoutiInit)
  
  for (i in 1:maxt)
  {
    pmat <- plant_mat_high
    h_off <- highHunting 
    
    p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
    agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
    plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
    #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
    plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
    
    #Summing the stages into 3 categories for better plotting
    plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
    plant_all <- cbind(plant_all, plant_mat_sum)
    
  }
  return(plant_mat_sum)
}

harvest_seq <- markovChain() #MaxT controls the marckovchev array.

#----------------------------stoch_growth----------------------------
stoch_growth <- function(){
  r <- numeric(maxt)
  plant_mat <- matrix(0, nrow = 17)
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  
  plant_mat <- plant_mat / sum(plant_mat) #Growth rate!
  agouti_vec <- c(agoutiInit)
  
  #markovChain()
  
  for (i in 1:maxt) #maxt controls markovchain array. 
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
    
    N <- sum(plant_mat)
    r[i] <- log(N)
    plant_mat <- plant_mat / N
  }
  
  loglambsim <- mean(r)
  
  return(loglambsim)
}

#===================================================================================================================
#Binary HeatMap
#===================================================================================================================
growthRate_mat<-matrix(0,21,21)
binary_mat<- matrix(0,21,21)
rownames(growthRate_mat) <- paste(seq(0,1,0.05)) #steps of .05
colnames(growthRate_mat) <- paste(seq(0,1,0.05))

rownames(binary_mat) <- paste(seq(0,1,0.05))
colnames(binary_mat) <- paste(seq(0,1,0.05))

num<-1 
num1<-1
for(i in seq(0,1,0.05)) #steps of .05
{
  high_harv[1,12:17] <- i # Multiplier for fecundity rate for Adult trees
  num1<-1
  for(j in seq(0,1,0.05)){
    
    high_harv[cbind(12:17,12:17)] <- j # Multiplier for survival rate of Adult trees
    plant_mat_low <- plant_S_mat #plat_s_mat is the original 17x17 matrix.
    plant_mat_high <- plant_S_mat * high_harv
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

heatmaply::heatmaply(binary_mat,margins=c(4,4), Rowv=NA, Colv=NA,xlab = "Adult Survival", ylab="Germination")
heatmaply::heatmaply(growthRate_mat,margins = c(4,4), Rowv=NA, Colv=NA,xlab = "Adult Survival", ylab="Germination", scale="none")

#saving these the heatmap as images 
#dir.create("heatMaps")
#library(heatmaply)
#heatmaply(growthRate_mat, file = "heatMaps/heatmaply_plot.png",Rowv=NA, Colv=NA,xlab = "Adult Survival", ylab="Germination")
#browseURL("heatMaps/heatmaply_plot.png")


#=============================================================================================================================
# Blue Heat Map
#=============================================================================================================================
#install.packages("akima")
library(akima)

filled.contour(x = seq(0,1,0.05),
               y = seq(0,1,0.05),
               z = growthRate_mat,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "Adult Survival",
               ylab = "Germination",
               key.title = title(main = "Growth Rate", cex.main = 0.5))


###==========================================================================
### ggplot heatmap example
###==========================================================================

growthRatedf <- reshape2::melt(growthRate_mat)
growthRatedf[,3] <- sample(seq(75,125, by=2.5)/100, replace=T, size=dim(growthRatedf)[1] ) # just to get some values in there without running the simulation

names(growthRatedf) <- c("Germination","AdultSurvival","Lambda") # not sure I got the ordering of these labels right

col.ramp <- c("goldenrod2","#a1dab4","#41b6c4","#2c7fb8","#253494") # create a color ramp; there are also built-in ones in ggplot ?scale_fill_continuous or ?scale_color_gradientn

p <- ggplot(growthRatedf, aes(x=Germination, y=AdultSurvival,fill=Lambda))
p <- p + geom_raster()
p <- p + scale_fill_gradientn(name=expression(lambda), colours=col.ramp, breaks=c(0.75,0.875,1,1.125,1.25)) +
  labs(x="Germination", y=expression(paste("Adult survivorship (", sigma[a], ")"))) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=14),legend.text=element_text(size=11),legend.title=element_text(size=14)) 
p

### Binary example

growthRatedf$LambdaBin <- ifelse(growthRatedf$Lambda > 1, 2, ifelse(growthRatedf$Lambda==1, 1, 0))
growthRatedf$LambdaBin <- factor(growthRatedf$LambdaBin, levels=c("0","1","2"), labels=c("Decline","Constant","Increase"))

p <- ggplot(growthRatedf, aes(x=Germination, y=AdultSurvival,fill=LambdaBin))
p <- p + geom_tile(color="white") # different version of raster, a touch slower, but permits outlines if that is useful
p <- p + scale_fill_manual(name=expression(paste(lambda," category:")),values = c("grey50","slategray3","darkcyan")) +
  labs(x="Germination", y=expression(paste("Adult survivorship (", sigma[a], ")"))) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=14),legend.text=element_text(size=11),legend.title=element_text(size=14)) 
p

# To save a plot to file:
# ggsave("~/LambdaBinomialAdultGermination.eps",plot=p,width=8, dpi=200, units="in") # or something like that. Can use jpeg, tiff, png, etc. I find jpeg and eps best.

#=====================================================================================================================
#Growth Rate vs animal Population
#=====================================================================================================================
agouti_Growth<- matrix(0, 1, 21)
agouti_pop<- matrix(0, 1, 21)
num<-1
for(i in seq(0,1,0.05))
{
  agouti_Growth[num]<- (i)
  agouti_pop[num]<- agouti_Abundance(agouti_Growth[num],1)
  print(c(agouti_Growth[num], agouti_pop[num]))
  num=num+1
}

plot(agouti_Growth,agouti_pop, xlab="Proportion of the Growth Rate", ylab="Animal Population", col="brown", ylim=c(0,4000), type="l",xlim=c(0,1),xaxs="i") 
#=======================================================================================================================
#Carrying Capacity vs Agouti Population
#=======================================================================================================================
agouti_capacity<- matrix(0, 1, 21)
agouti_pop<- matrix(0, 1, 21)
num<-1
for(i in seq(0,1,0.05))
{
  agouti_capacity[num]<-(i)
  agouti_pop[num]<- agouti_Abundance(1,agouti_capacity[num])
  print(c(agouti_capacity[num], agouti_pop[num]))
  num<-num+1
}

plot(agouti_capacity,agouti_pop, xlab="Proportion of the Carrying Capacity", ylab="Animal Population", col="brown", ylim=c(2000,4300), type="l",xlim=c(0.59,1),xaxs="i")


#=========================================================================================================================
#Agouti hunting vs. Plant abundance
#=========================================================================================================================
seedling_hunt_mat<- matrix(0,1,21)
sapling_hunt_mat<- matrix(0,1,21)
adults_hunt_mat<- matrix(0,1,21)

hunting_mat<- matrix(0,1,21)
num<-1
for(i in seq(0,1,0.05))
{
  hunting_mat[num]<-i 
  seedling_hunt_mat[num]<- plant_abundance_underHighHunt(hunting_mat[num])[1]
  sapling_hunt_mat[num]<- plant_abundance_underHighHunt(hunting_mat[num])[2]
  adults_hunt_mat[num]<- plant_abundance_underHighHunt(hunting_mat[num])[3]
  num= num+1
}

par(mar=c(5,4,1,1),oma=c(0,0,0,0))
plot(hunting_mat, log(seedling_hunt_mat), xlab="Proportion of Animal Hunted", ylab="log(Plant Population)", col="brown", ylim=c(0,35), type="l",xlim=c(0,1),xaxs="i")
lines(hunting_mat,log(sapling_hunt_mat), col="black")
lines(hunting_mat,log(adults_hunt_mat), col="green")
legend(0.1, 5, legend=c("Seedling","Sapling","Adults"),col=c("brown", "black","green"), lty=1:2, cex=0.5)








