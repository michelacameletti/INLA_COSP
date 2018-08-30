# INLA_COSP
Code for the simulation study presented in the paper "Bayesian modeling for spatially misaligned health and air pollution
data through the INLA-SPDE approach" by Michela Cameletti, Virgilio Gómez-Rubio and Marta Blangiardo.

The code is composed by the following files:
1) **Run_simulation.R**
  With this file all the libraries and additional source files are loaded. 
  Moreover, all the parameter values and settings for the simulation study are chosen.
  Finally, simulation function (contained in the **Simulation_singlestep.R** file, see 2)) is run in parallel
  (i.e. scenarios are run indipendently at the same time) or sequentially (i.e. one scenario at a time).
  
  
2) **Simulation_singlestep.R**
  This file contains the *Simulation_singlestep* function that for a given number of iterations (given by n_simulations):
  a) simulates data (*the true exposure values*) at the point and area level (using the function contained in the
    **Data_simulation.R** file); it also sample a given number of monitoring stations (whose locations remain fixed
    across simulations).
  b) uses the data available at the monitoring stations to run the spatial model (using the INLA-SPDE approach);
  c) estimates exposure at the area level using two linear combination methods for uspcaling: linear combination with
  neighbourhood intersections (Method1) and simple mean (Method2);
  d) runs the Poisson model considering three approaches for linking count (health) data with exposure estimates: 
    the plug-in (PI), the feed-forward (FF) and prior-exposure (PE) approach, each run under both Method1 
    and Method2.
  e) computes performance indexes (bias and RMSE) and some plots for the eposure model parameters, for the area 
     exposure prediction and for the Poisson parameters.
    

3) **Data_simulation.R**
This file contains the *Data_simulation* function that is run for a number of iterations (given by
*n_simulations*) for each scenario.
In particular this function:
    a) simulates the exposure values for a high-resolution grid of points;
    b) computes the true area exposure by averaging the exposure values observed for all the points inside each area;
    d) it defines a prediction grid (which is a sub-grid of the one defined at point a);
    c) computes the intersections + weights between the shapefile and the regular grid of prediction points
       (just for the first iteration); this is required by Method1;
    d) computes the prediction points inside each area (and the corresponding weights) or the closest one if 
       there are no points inside (just for the first iteration); this is required by Method2;
    e) draws a sample of monitoring stations from the set of points defined at point a) (their locations are 
       fixed across all the iterations of a given scenario) and compute the exposure values (given by the true exposure 
       value + a measurement error);
    f) simulates the Poisson count data by using the true area exposure.   
    
    

4) **Extra_functions.R**
This file contains some additional functions that are used by other functions. The most important one is *fit.inla*
which is used to implement the PE approach when running the Poisson model.


Required data files:
1) **Scenarios_bhicv.csv**
  This file contains the parameter values and the simulation setting according to the considered scenarios.
  

To run the code you must organize your working directory with the following sub-folders:
- **R_Code_LATEST** (this must contain all the *.R* and data files specified above);
- **Data** (here the simulated exposure data + weights matrices (for Method1 and Method2) + some plots will be saved for each scenario);
- **Figures** (here all the posterior analysis plots will be saved for each scenario);
- **Output** (here all the outputs and summary statistics will be saved for each scenario).


  
  
