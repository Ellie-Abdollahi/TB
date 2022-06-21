# Evaluation of TB elimination strategies in Canadian Inuit population 
Elaheh Abdollahi, 2022  
Agent-Based Modeling Laboratory
York University, Toronto, Canada


### How to use the codes

The model is run by 1) instantiating a `ModelParameters` object with the desired parameter values and 2) running the `main` function. 

The function `main` returns a population size x model time matrix where each element is one of the following model states: `(SUSC LAT LATT ACT ACTT DS)` coded from `1 ... 6`.  This matrix can be used to calculate the incidence and prevalence of the different model states. 


Since the model is stochastic, many realizations are required to get an accurate picture of the results. We recommend running this in parallel manner. This essentially means running the `main` function repeatedly, saving the results for each replicate. This can be done very easily using Julia's Distributed library. Simply `addprocs` or using `ClusterManagers` to connect to a cluster to launch `n` number of parallel workers. Then use `pmap` to run the main function on each worker, saving the results to be processed later. Here, we have defined run_and_save() function for running different scenarios in parallel mode. The arguments of this function are model parameters, LTBI treatment regimen, active treatment regimen, vaccine coverage (this is set to zero when there is no vaccine in the model), testing method, contact tracing and time-to-identification. To evaluate different scenarios set the desired arguments.
(See Final_Run.jl)

For the scenario in which average household size is reduced, “TB_AHS.jl” should be used instead of “TB_main.jl”, as the population structure changes after increasing the number of households.
