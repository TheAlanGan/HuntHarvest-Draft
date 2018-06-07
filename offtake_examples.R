# Hill and Padwe
  # Table 5.6
0.57 # high hunting ("sink")
0.037 # low hunting (Mbaracayu Reserve)

# Mena 2000
  # Table 4-3
agouti_mena_kills <- 98*c(7202,2972)/sum(c(7202,2972))*c(1/49, 1/81) # individuals harvested across entire study area (130km2 huntshed) -> kills/km2
  # Table 4-9
agouti_mena_density <- c(30, 14) # individuals per km2, low and high hunting sites
  # Calculating kill rates for "high" and "low"
agouti_mena_kills/agouti_mena_density # % mortality

# Naranjo 2004
  # Table 5
naranjo_kills <- 279 # total agoutis killed in Lancandon forest
naranjo_huntshed <- 11.5^2*pi
naranjo_kills / (14*naranjo_huntshed) # rough estimate of % mortality - very similar to Mena 2000

# Shaffer 2017
    # Table 2
shaffer_paca <- 198 # Agouti paca killed
shaffer_agouti <- 15 # Red-rumped agouti killed
    # Huntshed
shaffer_huntshed <- 6.3^2*pi
    # % Mortality
(shaffer_paca/shaffer_huntshed)/14
(shaffer_agouti/shaffer_huntshed)/16.7 # using Dasyprocta fuliginosa as surrogate from Mena
