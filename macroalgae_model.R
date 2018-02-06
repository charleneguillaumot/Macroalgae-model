# First draft model 
# without experimental approaches, just using data from the literature 
# Species of interests: 


#######----------------------------------------------------------------------------------------------------
## CHANGES IN BIOMASS, model Duarte et al 1997
#######----------------------------------------------------------------------------------------------------

# Make each variable depends on time (seasonal influence)

P ~ light + Temp + nutrients
R ~ Arrhenius # independent of plant size
PhotoR ~ P # independent of plant size 
Exu ~ P # independent of plant size
FrondBreak ~ current # computed as a function of average class lenght 
M # constant fraction of biomass
# Intraspecific competition neglicted due to high mortality and frond breakage; 
# influence of presence/absence of epiphytes is not significant


# combine this with changes in size-class densities
# definition of size-classes : 5; 5-10 ; 10-15 and >15

dBdt <- P - R - PhotoR - Exu - FrondBreak - M

#######----------------------------------------------------------------------------------------------------
## Model equation
#######----------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------
############ ENVIRONMENTAL DATA to force the model (data series) ###########
#---------------------------------------------------------------------------
# TPM: total suspended particulate materials  (mg.L-1)
TPM <- c(10,13,12,14,15,16,10)
# DIP: dissolved inorganic phosphorus (micromol.L-1)
DIP <- c(100,123,112,14,135,126,110)
# DIN: dissolved inorganic nitrogen (micromol.L-1)
DIN <- c(300,423,212,144,435,226,210)
# temperature
temp <- c(5,1,3,4,6,5,2)
# time
times <- seq(1,7,1)

# initial data 
Nint <- 100
Pint <- 200
B <- 5
Ggrowth <- 0.3


# what defines the time step : precision of the data time series 
for (i in 1:length(temp)){ # work maybe in runge kutta or euler? 
  
  
#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
###### Nutrients concentration #######
#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-

# Update of the remaining nutrients concentration within the plant tissues 
# dynamics in uptake of nutrients 
Next= DIN[i] # external nutrient concentration in water 
Kn = 2 # half saturation constant for the uptake of nutrients  (micromol.L-1)
Nimax= 3000 # maximum internal concentration (micromol.gDW-1)
Nimin= 500 # minimum internal concentration (micromol.gDW-1)
VmaxN= 90 #maximum uptake rate of the nutrient (micromol.ngDW-1. d-1)

Pext=DIP[i] # external nutrient concentration in water 
Kp = 0.1 # half saturation constant for the uptake of nutrients (micromol.L-1)
Pimax= 250 # maximum internal concentration (micromol.gDW-1)
Pimin= 31 # minimum internal concentration (micromol.gDW-1)
VmaxP= 7 #maximum uptake rate of the nutrient (micromol.ngDW-1. d-1)

# nutrients assimilation 
PsiN= ((Nimax-Nint)/(Nimax-Nimin))*VmaxN*(Next/(Kn+Next)) 
PsiP= ((Pimax-Pint)/(Pimax-Pimin))*VmaxP*(Pext/(Kp+Pext))
# Xext/(Kx+Xext) Michaelis Menten kinetics
# (Ximax-Xint)/(Ximax-Ximin) accounts for the internal nutrient concentration 

LambdaN = Nint * Ggrowth # use of internal nutrient (N) concentration
LambdaP = Pint * Ggrowth # use of internal nutrient (P) concentration

# remaining internal nutrient concentration 
Nint= PsiN-LambdaN
Pint= PsiP-LambdaP

# the model considers that the nutrients are not limited by the external content 
# growth is dependent of internal concentration of nutrients 
# NP:internal concentration of nutrients 
fN= 1-(Nimin-Nint)
fP= 1-(Pimin-Pint)

# according to time, we measure the amount of nutrients that are available within the algae by substracting the amount of nutrients consumed (lambda) from the uptake of nutrients (psi)
fNP= min (fN, fP)

# if Pimin ou Nimin < Pint ou Nint, we have f(N) or f(P)=0, du coup f(NP)=0, the growth stops


#___________________
###### LIGHT #######
#___________________
# Case 1: Zhang 2016, seabed light
k= 0.0484*TPM[i]+0.0243 # light coefficient extinction of the station, in Zhang2016, they use an empirical relationship (k=0.0484*TPM+0.0243), where TPM= total particulate materials 
z= 0.2 # depth (m)
Is= 200.38-116.47*cos(2*pi*(times[i]-1)/365)# light intensity at the surface
I=Is*exp(-k*z)
I0=180 # W.m2, optimum light intensity for growth
fI= (I/I0)*exp(1-(I/I0))

# Case 2 : Clark 2013, seabed light
#ln(i) = sin(d/365.2*pi)+cos(d/365.2*pi)+d
fI =exp(sin(times[i]/365.2*pi)+cos(times[i]/365.2*pi)+times[i])
# light at the seabed per day
# i: irradiance at the seabed (mol photons/m2/d)
# times= date=> we should accordate the time with the date of the year because the model is working like this 

#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
###### INFLUENCE OF TEMPERATURE #######
#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
#temperature-optimum curve
b= 0.5 #adjustment parameter (°C) # how can we find it? 
Tw=temp[i] # sea water T° (°C)
Topt= 13 # otpimum T° for growth (°C)
Tmax= 23 # upper temperature limit above which growth ceases (°C)
Xt= (Tw-Topt)/(Tw-Tmax)
fT=(2*(1+b)*Xt)/(Xt^2 +2*b*Xt + 1)


#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
###### Respiration rates #######
#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
Rmax20= 0.015 #max respiration rate at 20°C (d-1)
r= 1.07 # empirical coefficient
Resp = Rmax20 * r^(temp[i]-20)

#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
###### Growth rate equation #######
#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
mumax=0.135 # max growth rate (d-1)
Ggrowth = mumax*fT*fNP*fI
NGR = Ggrowth - Resp

#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
###### Erosion rate #######
#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
ER=0.01/100 #individual erosion rate per kelp (d-1)


#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-
###### GLOBAL EQUATION BIOMASS CHANGE #######
#-_-_-_-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-

dBdt = (NGR-ER)*B }




#### Equation LEO
# integrate in Zhang 2016 equation the velocity (Ren 2014)
B(t+1) = R(B(t))* G(~light,depth, sediments, resp) - IScour(~ depth, distance from glacier) - M + epsilon

