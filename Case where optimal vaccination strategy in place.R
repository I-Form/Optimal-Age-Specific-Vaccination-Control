# We are simulating the behaviour of a situation where the optimal vaccination strategy is in place. 
# Our code determines the optimal control and applies it to get the optimal trajectory.

library(deSolve)
library(tidyverse)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - - - Initial Values- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

#We divide the total population in two age groups, over and under 65. 
#For brevity we call the two populations older and younger repsectively.

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
young_nv <- 0.23      #Percentage of younger people unwilling to get vaccinated
old_nv <- 0.07        #Percentage of older people unwilling to get vaccinated
mean_h.t_inf <- 7.4   #Mean holding time in the infectious state
mean_h.t_exp <- 6.6   #Mean holding time in the exposed state
R0_oo <- 1.2          #Mean number of older people an older person can infect before recovering
R0_oy <- 0.9          #Mean number of younger people an older person can infect before recovering
R0_yo <- 0.9          #Mean number of older people an younger person can infect before recovering
R0_yy <- 1.2          #Mean number of younger people an younger person can infect before recovering
mean_h.t_vacc <- 14   #Mean holding time before the vaccine becomes effective
effectiveness <- 0.9  #Vaccine effectiveness

time_horizon <- 300 #Number of days for simulation

d <- 0.0005     #Threshold for convergence to optimal trajectory
error_term <- 1 #Random initialisation to enter the while loop

upp.control <- 0.3 #Upper bound for control

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - Creating state and parameter vectors- - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

State <- c(OS = (1-old_nv)*(O_Pop-(old_exp+old_inf+old_rec)),             #Initial number of susceptinble people in the older age group
           OV =0, ON =0,                                                  #Initial values for states associated with vaccination in the older age group set to 0
           OU =old_nv*(O_Pop-(old_exp+old_inf+old_rec)),                  #Initial number of people refusing the vaccine in the older age group               
           OE =old_exp, OI =old_inf, OR =old_rec,                         #Initial values for exposed, infectious and recovered states in the older age group
           OP =0,                                                         #Initial value for protected people from the vaccine set to 0 in the older age group
           YS = (1-young_nv)*(Y_Pop-(young_exp+young_inf+young_rec)),     #Initial number of susceptinble people in the younger age group
           YV =0, YN =0,                                                  #Initial values for states associated with vaccination in the younger age group set to 0
           YU =young_nv*(Y_Pop-(young_exp+young_inf+young_rec)),          #Initial number of people refusing the vaccine in the younger age group
           YE =young_exp, YI =young_inf, YR =young_rec,                   #Initial values for exposed, infectious and recovered states in the younger age group
           YP =0,                                                         #Initial value for protected people from the vaccine set to 0 in the younger age group
           LOS = 0, LOV = 0, LON = 0, LOU = 0, LOE = 0, LOI = 0,          #All adjoint states are intialised to 0 and will be solved backwards in time
           LYS = 0, LYV = 0, LYN = 0, LYU = 0, LYE = 0, LYI = 0) 
           
           
Parameters <- c(betaoo = R0_oo/(mean_h.t_exp+mean_h.t_inf+mean_h.t_vacc), #Rate at which an older person becomes infected by a older person       
                betayo = R0_yo/(mean_h.t_exp+mean_h.t_inf+mean_h.t_vacc), #Rate at which an older person becomes infected by a younger person 
                betaoy = R0_oy/(mean_h.t_exp+mean_h.t_inf+mean_h.t_vacc), #Rate at which an younger person becomes infected by a older person 
                betayy = R0_yy/(mean_h.t_exp+mean_h.t_inf+mean_h.t_vacc), #Rate at which an younger person becomes infected by a younger person 
                gammaV = 1/mean_h.t_vacc,                                 #Rate at which vaccine becomes effective
                alphaV = effectiveness,                                   #Vaccine effectiveness
                gammaE = 1/mean_h.t_exp,                                  #Rate at which an exposed person becomes infectious
                gammaI = 1/mean_h.t_inf,                                  #Rate at which an infectious person recovers
                WO = 10^11, WY = 10^11)                                   #Age specific weight constants 


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - Time sequence and intial control- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

# The time sequence for the state and adjoint systems. These will be used in the ODE solver. 
times_state <- seq(from=0,to=time_horizon,by=1)
times_adjoint <- seq(from=time_horizon,to=0,by=-1) #This is defined backwards because the adjoint equations are solved backwards in time

# Initialisation of control function to all zeros. This function will be updated until convergence.

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
  LOS = state [17]
  LOV = state[18]
  LON = state[19]
  LOU = state[20]
  LOE = state[21]
  LOI = state[22]
  LYS = state[23]
  LYV = state[24]
  LYN = state[25]
  LYU = state[26]
  LYE = state[27]
  LYI = state[28]
  
  #Extract parameter values
  betaoo = parms[1]
  betayo = parms[2]
  betaoy = parms[3]
  betayy = parms[4]
  gammaV = parms[5]
  alphaV = parms[6]
  gammaE = parms[7]
  gammaI = parms[8]
  WO = parms[9]
  WY = parms[10]
  
  # Total number of individuals in each age group
  OT <- OS + OV + ON + OU + OE + OI + OR + OP
  YT <- YS + YV + YN+ YU + YE + YI + YR + YP
  
  #Read in the control function 
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
  dLOS <- 0
  dLYS <- 0
  dLOV <- 0
  dLYV <- 0
  dLON <- 0
  dLYN <- 0
  dLOU <- 0
  dLYU <- 0
  dLOE <- 0
  dLYE <- 0
  dLOI <- 0
  dLYI <- 0
  
  dxdt <- c(dOS, dOV, dON, dOU, dOE, dOI, dOR, dOP, dYS, dYV, dYN, dYU, dYE, dYI, dYR, dYP,
            dLOS, dLOV, dLON, dLOU, dLOE, dLOI, dLYS, dLYV, dLYN, dLYU, dLYE, dLYI )
  list(dxdt)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - -Solve the state ODEs forwards in time- - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

ode(
  func=state_system,
  y=State,
  times=times_state,
  parms=Parameters
) %>%
  as.data.frame() -> out_state

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - -Define the adjoint ODEs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

adjoint_system <- function(t, state, parms){
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
  LOS = state [17]
  LOV = state[18]
  LON = state[19]
  LOU = state[20]
  LOE = state[21]
  LOI = state[22]
  LYS = state[23]
  LYV = state[24]
  LYN = state[25]
  LYU = state[26]
  LYE = state[27]
  LYI = state[28]
  
  #Extract parameter values
  betaoo = parms[1]
  betayo = parms[2]
  betaoy = parms[3]
  betayy = parms[4]
  gammaV = parms[5]
  alphaV = parms[6]
  gammaE = parms[7]
  gammaI = parms[8]
  WO = parms[9]
  WY = parms[10]
  
  # Total number of individuals in each age group
  OT <- OS + OV + ON + OU + OE + OI + OR + OP
  YT <- YS + YV + YN+ YU + YE + YI + YR + YP
  
  #Read in the control function 
  uOt <- uO(t)
  uYt <- uY(t)
  
  #Adjoint equations
  dOS <- 0
  dOV <- 0
  dON <- 0
  dOU <- 0
  dOE <- 0
  dOI <- 0
  dOR <- 0
  dOP <- 0
  dYS <- 0
  dYV <- 0
  dYN <- 0
  dYU <- 0
  dYE <- 0
  dYI <- 0
  dYR <- 0
  dYP <- 0
  dLOS <- (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT + uOt)*LOS - uOt*LOV - (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*LOE 
  dLYS <- (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT + uYt)*LYS - uYt*LYV - (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*LYE 
  dLOV <- (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT + gammaV)*LOV - (1-alphaV)*gammaV*LON - (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*LOE
  dLYV <- (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT + gammaV)*LYV - (1-alphaV)*gammaV*LYN - (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*LYE
  dLON <- (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*(LON-LOE)
  dLYN <- (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*(LYN-LYE)
  dLOU <- (betaoo*(OE+OI)/OT + betayo*(YE+YI)/YT)*(LOU-LOE)
  dLYU <- (betaoy*(OE+OI)/OT + betayy*(YE+YI)/YT)*(LYU-LYE)
  dLOE <- betaoo*OS/OT*(LOS-LOE) + betaoo*OV/OT*(LOV-LOE) + betaoo*ON/OT*(LON-LOE) + betaoo*OU/OT*(LOU-LOE) + 
    betaoy*YS/OT*(LYS-LYE) + betaoy*YV/OT*(LYV-LYE) + betaoy*YN/OT*(LYN-LYE) + betaoy*YU/OT*(LYU-LYE) + gammaE*LOE - gammaE*LOI
  dLYE <- betayo*OS/YT*(LOS-LOE) + betayo*OV/YT*(LOV-LOE) + betayo*ON/YT*(LON-LOE) + betayo*OU/YT*(LOU-LOE) + 
    betayy*YS/YT*(LYS-LYE) + betayy*YV/YT*(LYV-LYE) + betayy*YN/YT*(LYN-LYE) + betayy*YU/YT*(LYU-LYE) + gammaE*LYE - gammaE*LYI
  dLOI <- betaoo*OS/OT*(LOS-LOE) + betaoo*OV/OT*(LOV-LOE) + betaoo*ON/OT*(LON-LOE) + betaoo*OU/OT*(LOU-LOE) + 
    betaoy*YS/OT*(LYS-LYE) + betaoy*YV/OT*(LYV-LYE) + betaoy*YN/OT*(LYN-LYE) + betaoy*YU/OT*(LYU-LYE) + gammaI*LOI -1
  dLYI <- betayo*OS/YT*(LOS-LOE) + betayo*OV/YT*(LOV-LOE) + betayo*ON/YT*(LON-LOE) + betayo*OU/YT*(LOU-LOE) + 
    betayy*YS/YT*(LYS-LYE) + betayy*YV/YT*(LYV-LYE) + betayy*YN/YT*(LYN-LYE) + betayy*YU/YT*(LYU-LYE) + gammaI*LYI - 1
  
  dxdt <- c(dOS, dOV, dON, dOU, dOE, dOI, dOR, dOP, dYS, dYV, dYN, dYU, dYE, dYI, dYR, dYP,
            dLOS, dLOV, dLON, dLOU, dLOE, dLOI, dLYS, dLYV, dLYN, dLYU, dLYE, dLYI )
  list(dxdt)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - -Solve the adjoint ODEs backwards in time - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

ode(
  func=adjoint_system,
  y=as.numeric(out_state[length(out_state),-1]), #The input is the output from the state ODEs
  times=times_adjoint,
  parms=Parameters
) %>%
  as.data.frame() -> out_adjoint


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - -Repeat process until convergence - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

# We solve the system of ODEs for the state and adjoint systems, updating the control every time, until convergence
# Convergence is checked by a randomly initialised error term that we want to be under a predetermined threshold 

while (error_term > d) {
  
  # Control update obtain by solving the optimal control problem
  new_uO <- out_state[,"OS"]/Parameters["WO"]*(out_adjoint[,17]- out_adjoint[,18])
  new_uY <- out_state[,"YS"]/Parameters["WY"]*(out_adjoint[,23]- out_adjoint[,24])
 
  # Making sure the control functions are under the upper bound 
  if(any(new_uO>upp.control)){ 
    new_uO <- rep(upp.control,length(times_state))
  }else if(any(new_uO<0)){
    new_uO <- 0
  }
  
  if(any(new_uY>upp.control)){ 
    new_uY <- rep(upp.control,length(times_state))
  }else if(any(new_uY<0)){
    new_uY <- 0
  }
  
  # Turning the updated controls into functions to be called in the ODE solver
  signal_uO = as.data.frame(list(times = times_state, import = new_uO))
  signal_uY = as.data.frame(list(times = times_state, import = new_uY))
  uY = approxfun(signal_uY, rule = 2)
  uO = approxfun(signal_uO, rule = 2)
  
  # Saving the previous output of the solver to use for the error term
  state_previous <- as.matrix(out_state[,-1]) 
  dimnames(state_previous)<- NULL 
  
  # Solving the state and adjoint equations with the updated control
  ode(
    func=state_system,
    y=State,
    times=times_state,
    parms=Parameters
  ) %>%
    as.data.frame() -> out_state
  
  ode(
    func=adjoint_system,
    y=as.numeric(out_state[length(out_state),-1]),
    times=times_adjoint,
    parms=Parameters
  ) %>%
    as.data.frame() -> out_adjoint
  
  state_current <- as.matrix(out_state[,-1]) # Needed for the error term
  dimnames(state_current)<- NULL
  
  # Error term used for the while loop
  error_term <- norm(state_previous-state_current,type = "1")/norm(state_previous,type = "1")
}  


# Getting the total percentage of the population that got infected in each age group and in total
av_or<- out_state[,c("OR")] # All the older people that recovered (they are all that have been infected)
last(av_or)/O_Pop # Obtain percentage
av_yr<- out_state[,c("YR")] # All the younger people that recovered (they are all that have been infected)
last(av_yr)/Y_Pop # Obtain percentage

(last(av_or)+last(av_yr))/(O_Pop+Y_Pop) # Total percentage of the population that got infected

# Getting the percentages of the populations that were succesfully protected from the vaccine.
av_op<- out_state[,c("OP")] # All the older people that are protected
last(av_op)/O_Pop # Obtain percentage
av_yp<- out_state[,c("YP")] # All the younger people that are protected
last(av_yp)/Y_Pop # Obtain percentage

(last(av_op)+last(av_yp))/(O_Pop+Y_Pop)  #Total percentage of the population that got succesfully vaccinated

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#- - - - - - - - - - - - - - - - - - Plotting the optimal trajectories - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -# 

# Extracting the states we want to plot 
state_to_plot_vaccine<- as.data.frame(out_state[,c("time","OS","YS","OI","YI","OR","YR")])
colnames(state_to_plot_vaccine)<- c("time","O65 Susceptible","Y65 Susceptible","O65 Infectious",
                            "Y65 Infectious","O65 Recovered","Y65 Recovered")  

# Plotting the optimal trajectories of the states after the optimal vaccination strategy is applied 
state_to_plot_vaccine %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  scale_y_log10()+ # This is to be commented out in the case of high transmission rates.
  scale_colour_hue(name="Compartment")+
  labs(x='time (days)',y='number of individuals')
ggsave("With vaccine all states.pdf")

# Extracting the optimal controls to plot. They are in the form of percentages
uO<-signal_uO[,"import"]
uY<-signal_uY

# Plot both optimal controls together in logscale
cbind(uY,uO)%>%
  gather(variable,value,-times) %>%
  ggplot(aes(x=times,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  scale_y_log10()+ 
  labs(x='time (days)',y='Percentage to be vaccinated (log scale)')+
  scale_colour_hue(name="Control functions",
                   labels=c("Under 65", "Over 65"))
ggsave("restrictions, inputs, log scale.pdf")

# Plotting the percentages of older people to get vaccinated at each time point
signal_uO %>%
  ggplot( aes(x=times, y=import)) +
  geom_line() +
  geom_point()+
  theme_classic()+
  labs(x='time (days)',y='uO')
ggsave("restrictions uO.pdf")

# If we want to plot the numbers of older people to vaccinate as opposed to the percentages we use the following
new<-data.frame(signal_uO$times,signal_uO$import*O_Pop)
colnames(new)<-c("times","import")
new %>%
  ggplot( aes(x=times, y=import)) +
  geom_line() +
  geom_point()+
  theme_classic()+
  labs(x='time (days)',y='uO')

# Plotting the percentages of younger people to get vaccinated at each time point
signal_uY %>%
  ggplot( aes(x=times, y=import)) +
  geom_line() +
  geom_point() +
  theme_classic()+
  labs(x='time (days)',y='uY')
ggsave("restrictions uY.pdf")

# Plotting the number of younger people to vaccinate as opposed to percentage
newY<-data.frame(signal_uY$times,signal_uY$import*O_Pop)
colnames(newY)<-c("times","import")
newY %>%
  ggplot( aes(x=times, y=import)) +
  geom_line() +
  geom_point()+
  theme_classic()+
  labs(x='time (days)',y='uY')

# Extracting states to use in comparison plots with the baseline situation without vaccine
v_osu<- new_out_state[,"O65 Susceptible"]
v_ysu<- new_out_state[,"Y65 Susceptible"]
v_oi<- new_out_state[,"O65 Infectious"]
v_yi<- new_out_state[,"Y65 Infectious"]
v_or<- new_out_state[,"O65 Recovered"]
v_yr<- new_out_state[,"Y65 Recovered"]

# The following comparison plots can only be run after running the baseline situation code first

# Comparison plot of the over 65 infectious with and without vaccine
cbind(nv_oi,v_oi)%>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  theme(legend.position="none")+
  scale_y_log10()+
  scale_colour_hue(name="Infectious individuals",
                   labels=c("Without vaccine", "With vaccine"),direction = 1)+
  labs(x='time (days)',y='number of individuals ')
ggsave("restrictions, infected, over 65.pdf")

# Comparison plot of the under 65 infectious with and without vaccine
cbind(nv_yi,v_yi)%>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  scale_y_log10()+
  scale_colour_hue(name="Infectious individuals",
                   labels=c("With vaccine", "Without vaccine"),direction = -1)+
  labs(x='time (days)',y='number of individuals ')
ggsave("restrictions, infected, under 65.pdf")

# We can produce more comparison plots in a similar way. 
