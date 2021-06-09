# Optimal-Age-Specific-Vaccination-Control

This code can be used to simulate the evolution of the COVID-19 pandemic with and without optimal vaccination control in place using our proposed compartmental model. A preprint containing the model can be found on ArXiv and Researchgate under the title "Optimal age-specific vaccination control for COVID-19".

The code was written in R. 

# Baseline situation without vaccination.R

This code can be used to simulate the evolution of the pandemic if there is no vaccination strategy being applied. We use initial conditions representative of the situation in the republic of Ireland at the beginning of 2021, but any sort of initialisation can be used. The model used is the same as in the case of vaccination being present and all states and parameters related to the vaccination are being set to zero. At the end of the code we can save some information from this simulation to use in comparative plots within the next code. 

# Case where optimal vaccination strategy in place.R

This code can be used to obtain the optimal vaccination strategy given specific intial conditions (same as the case without vaccination) and simulate the evolution of the pandemic with the vaccination strategy in place. We define the state and adjoint equations and then solve them repeatedly, updating the optimal control until convergence. At the end of the code we can produce some comparative plots to see the effect of the optimal vaccination, using the information saved from the baseline code.
