# File Name: AdultSurvival_low-HighHunt/Harv_Transition
Description: Generates a contour map we vary: 
1. The probability of transition from high to low, introduced variable p in the markovChain_highTransition function
2. Adult survival multiplier (highHarvestSurvival)


# File Name: AdultSuvival_HighHunting_Germina(constant)
Description: Generates contour maps under varying germination levels (0.25, 0.5, 0.75, 1)
For each contour map we vary:
1. Constant high hunting with high/low harvest 
2. Adult Survival (highHarvestSurvival multiplier)


# File Name: AS_Germination
Description:Generates a contour map, color gradient represents stochastic growth rate of Brazil nut 
For each contour map we vary:
1. Adult Survival (highHarvestSurvival multiplier)
2.Germination (highHarvestFecundity multiplier)

# File Name: constantHunting_AS
Description:Generates a contour map, color gradient represents stochastic growth rate of Brazil nut 
For each contour map we vary:
1.Adult Survival (highHarvestSurvival multiplier)
2.Constant High Hunting (highHunting multiplier) wth H/L harvest 

# File Name: Rmax_AdultSurvival_BrazilNut
Description:Generates a contour map, color gradient represents stochastic growth rate of Brazil nut 
For contour map we vary:
1. Adult Survival (highHarvestSurvival multiplier)
2. Rmax range (0-3) 


# File Name: Delta_AdultSurvival (Growth Rate/Seed disperser pop)
Description: Generates two contour maps
Contour map #1 : Delta and Adult Survival- color gradient represents the agouti population
We vary:
1. Adult Survival (highHarvestSurvival multiplier)
2. Delta animal dependence on disperser-> range (0-1) 
 
Contour map #2 : Delta and Adult Survival- color gradient represents the stochastic growth rate of Brazil Nut
We vary:
1. Adult Survival (highHarvestSurvival multiplier)
2. Delta animal dependence on disperser-> range (0-1) 


# File Name: germination_AdultSurvival_heatmap(changed in matrix)
Generates a contour map
Contour map : Log(Germination) and Adult Survival- color gradient represents the stochastic growth rate of Brazil Nut
We vary:
1. Adult Survival (highHarvestSurvival multiplier)
2. Germination - range(log(1)- log(1000)) values changed within the plant transition matrix


# File Name: High_LowHuntingHarvest-Germination&AdultSurvival
Description:Generates a contour map: Color gradient represents Brazil nut stochastic growth rate 
We vary:
1. High/low hunting (highHunting/lowHunting multiplier) with H/L harvest 
2. Adult Survival (highHarvestSurvival multiplier)


# File Name: HighHunting/Harvest-Germination&AdultSurvival
Description: Generates a contour map: Color gradient represents Brazil nut stochastic growth rate 
We vary:
1. Constant high hunting and harvest (highHunting multiplier and pmat under high harvest) 
2. Adult Survival (highHarvestSurvival multiplier)


# File Name: HighHunting&H/LHarvest-Germination&AdultSurvival
Description:Generates a contour map: Color gradient represents Brazil nut stochastic growth rate 
We vary:
1. Constant high hunting (highHunting multiplier) with H/L harvest 
2. Adult Survival (highHarvestSurvival multiplier)


# File Name: Hunting_EigenValues
Description:Generates a line plot: Sustainable hunting threshold for seed disperser 
-> x-axis = Constant high hunting (highHunting multiplier)
-> y-axis = Brazil nut stochastic growth rate 

# File Name: Rmax_AdultSurvival_BrazilNut
Description:Generates a contour map, color gradient represents stochastic growth rate of Brazil nut 
We vary:
1. Adult Survival (highHarvestSurvival multiplier)
2. Rmax range (0-3) 

# File Name: Rmax_AS_Chamaedorea elegans
Description: 
Uses chamaedorea elegans plant transition matrix
Generates a  contour map: Color gradient represents stochastic growth rate 
We vary:
1. Adult Survival highHarvestSurvival multiplier
2. Rmax range(0-3)


# File Name: Rmax_AS_EutrepeEdulis
Description: 
Uses Eutrepe Edulis plant transition matrix
Generates a  contour map: Color gradient represents stochastic growth rate 
We vary:
1. Adult Survival highHarvestSurvival multiplier
2. Rmax range(0-3)

# File Name: Rmax_AS_PrunusAfricana
Description: 
Uses Prunus Africana plant transition matrix
Generates a  contour map: Color gradient represents stochastic growth rate 
We vary:
1. Adult Survival highHarvestSurvival multiplier
2. Rmax range(0-3)

# File Name: Rmax_AS_Garcinia Lucida
Description: 
Uses Garcinia Lucida plant transition matrix
Generates a  contour map: Color gradient represents stochastic growth rate 
We vary:
1. Adult Survival highHarvestSurvival multiplier
2. Rmax range(0-3)

# File Name: Rmax_AS_Genoma_Deversa
Description: 
Uses Genoma Deversa plant transition matrix
Generates a  contour map: Color gradient represents stochastic growth rate 
We vary:
1. Adult Survival highHarvestSurvival multiplier
2. Rmax range(0-3)


# File Name: Harvest_Hunt_DettaP, DeltaD
Description: Generating 25 contour plots as a single plot
Overall varying delta p and delta d
-> Plot x axis: Delta d
-> Plot y axis: Delta p

For each individual contour plot we are varying harvest (1- Adult Survival multiplier) and hunting (highHunting multiplier) levels
-> Contour plot x axis: Harvest
->Contour plot y axis: Hunting 

We are using Brazil nut and agouti specific parameter values 





