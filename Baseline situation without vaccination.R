# We are simulating the behaviour of a baseline situation where no vaccination is in place. 
# We use the same model as the one for the vaccination and assign zeros to all states and parameters associated with the vaccination.

library(deSolve) 
library(tidyverse)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - - - Initial Values- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

# We divide the total population in two age groups, over and under 65. 
# For brevity we call the two populations older and younger repsectively.

# Initial state values 
Y_Pop <- 4000000    #Total younger susceptible population
O_Pop <- 900000     #Total older susceptible population
young_exp <- 2000   #Number of younger exposed
old_exp <- 200      #Number of older exposed
young_inf <- 2000   #Number of younger infected
old_inf <- 200      #Number of older infected 
young_rec <- 200000 #Number of recovered younger people 
old_rec <- 100000   #Number of recovered older people 

# Initial parameter values
young_nv <- 0.23    #Percentage of younger people unwilling to get vaccinated
old_nv <- 0.07      #Percentage of older people unwilling to get vaccinated
mean_h.t_inf <- 7.4 #Mean holding time in the infectious state
mean_h.t_exp <- 6.6 #Mean holding time in the exposed state
R0_oo <- 1.2        #Mean number of older people an older person can infect before recovering
R0_oy <- 0.9        #Mean number of younger people an older person can infect before recovering
R0_yo <- 0.9        #Mean number of older people an younger person can infect before recovering
R0_yy <- 1.2        #Mean number of younger people an younger person can infect before recovering

time_horizon <- 300 #Number of days for simulation

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - Creating state and parameter vectors- - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

State <- c(OS = O_Pop-(old_exp+old_inf+old_rec),         #Initial number of susceptinble people in the older age group
           OV =0, ON =0, OU =0,                          #Initial values for states associated with vaccination in the older age group set to 0
           OE =old_exp, OI =old_inf, OR =old_rec,        #Initial values for exposed, infectious and recovered states in the older age group
           OP =0,                                        #Initial value for protected people from the vaccine set to 0 in the older age group
           YS = Y_Pop-(young_exp+young_inf+young_rec),   #Initial number of susceptinble people in the younger age group
           YV =0, YN =0, YU =0,                          #Initial values for states associated with vaccination in the younger age group set to 0
           YE =young_exp, YI =young_inf, YR =young_rec,  #Initial values for exposed, infectious and recovered states in the younger age group
           YP =0)                                        #Initial value for protected people from the vaccine set to 0 in the younger age group

Parameters <- c(betaoo = R0_oo/(mean_h.t_exp+mean_h.t_inf), #Rate at which an older person becomes infected by a older person       
                betayo = R0_yo/(mean_h.t_exp+mean_h.t_inf), #Rate at which an older person becomes infected by a younger person 
                betaoy = R0_oy/(mean_h.t_exp+mean_h.t_inf), #Rate at which an younger person becomes infected by a older person 
                betayy = R0_yy/(mean_h.t_exp+mean_h.t_inf), #Rate at which an younger person becomes infected by a younger person 
                gammaV = 0, alphaV = 0,                     #Vaccine-related parameters, rate at which vaccine becomes effective, vaccine effectiveness
                gammaE = 1/mean_h.t_exp,                    #Rate at which an exposed person becomes infectious
                gammaI = 1/mean_h.t_inf,                    #Rate at which an infectious person recovers
                WO = 10^11, WY = 10^11)                     #Age specific weight constants (not used in this part of the simulation as no control is applied)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - Time sequence and control - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

# The time sequence for the state system. These will be used in the ODE solver. 
times_state <- seq(from=0,to=time_horizon,by=1)

# Initialisation of control function to all zeros. This function will not be updated in this simulation as we are looking into the baseline model where no 
# vaccination strategy is in place. 

# We use the R function approxfun to turn the sequence of zeros to a function of time so that we can use it in the ODE solver
signal_uO = as.data.frame(list(times = times_state, import = rep(0,length(times_state))))
signal_uY = signal_uO
uY = approxfun(signal_uY, rule = 2)
uO = approxfun(signal_uO, rule = 2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - -Define the state ODEs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

state_system <- function(t, state, parms){
  #Extract state variables
  OS = state[1]
  OV = state[2]
  ON = state[3]
  OU = state[4] 
  OE = state[5] 
  OI = state[6] 
  OR = state[7] 
  OP = state[8]
  YS = state[9] 
  YV = state[10] 
  YN = state[11] 
  YU = state[12] 
  YE = state[13] 
  YI = state[14] 
  YR = state[15] 
  YP = state[16]
  
  #Extract parameter values
  betaoo = parms[1]
  betayo = parms[2]
  betaoy = parms[3]
  betayy = parms[4]
  gammaV = parms[5]
  alphaV = parms[6]
  gammaE = parms[7]
  gammaI = parms[8]
  
  # Total number of individuals in each age group
  OT <- OS + OV + ON + OU + OE + OI + OR + OP
  YT <- YS + YV + YN + YU + YE + YI + YR + YP
  
  # Exctracting the control functions 
  uOt <- uO(t)
  uYt <- uY(t)
  
  #State equations
  dOS <- - (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*OS-uOt*OS
  dOV <- uOt*OS - (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*OV - (1-alphaV)*gammaV*OV - alphaV*gammaV*OV
  dON <- (1-alphaV)*gammaV*OV - (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*ON
  dOU <- -(betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*OU
  dOE <- (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*(OS + OV + ON + OU) - gammaE*OE
  dOI <- gammaE*OE - gammaI*OI
  dOR <- gammaI*OI
  dOP <- gammaV*alphaV*OV
  dYS <- - (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*YS-uYt*YS
  dYV <- uYt*YS - (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*YV - (1-alphaV)*gammaV*YV - alphaV*gammaV*YV
  dYN <- (1-alphaV)*gammaV*YV - (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*YN
  dYU <- -(betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*YU
  dYE <- (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*(YS + YV + YN + YU) - gammaE*YE
  dYI <- gammaE*YE - gammaI*YI
  dYR <- gammaI*YI
  dYP <- gammaV*alphaV*YV
  
  dxdt <- c(dOS, dOV, dON, dOU, dOE, dOI, dOR, dOP, dYS, dYV, dYN, dYU, dYE, dYI, dYR, dYP)
  list(dxdt)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - Solving the system of ODEs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

ode(
  func=state_system,
  y=State,
  times=times_state,
  parms=Parameters
) %>%
  as.data.frame() -> out_state

# Getting the total percentage of the population that got infected in each age group and in total
anv_or<- out_state[,c("OR")] # All the older people that recovered (they are all that have been infected)
last(anv_or)/O_Pop # Obtain percentage
anv_yr<- out_state[,c("YR")] # All the younger people that recovered (they are all that have been infected)
last(anv_yr)/Y_Pop # Obtain percentage

(last(anv_or)+last(anv_yr))/(O_Pop+Y_Pop) #Total percentage of the population that got infected

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - Plotting the trajectories - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

#Choose the states we want to plot
state_to_plot <- as.data.frame(out_state[,c("time","OS","YS","OI","YI","OR","YR")])

#Produce the plot
state_to_plot%>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  theme(legend.position="none")+ 
  #scale_y_log10() +  # Logscale can be used depending on the specific simulation and specifically the transmission rates
  labs(x='time (days)',y='number of individuals')

#Save plot in pdf file
ggsave("restrictions, all state, logscale.pdf")

# Saving results to use when plotting comparison plots in the code with the vaccine in place 
nv_os<- state_to_plot[,c("time","OS")]
nv_ys<- state_to_plot[,c("time","YS")]
nv_oi<- state_to_plot[,c("time","OI")]
nv_yi<- state_to_plot[,c("time","YI")]
nv_or<- state_to_plot[,c("time","OR")]
nv_yr<- state_to_plot[,c("time","YR")]
