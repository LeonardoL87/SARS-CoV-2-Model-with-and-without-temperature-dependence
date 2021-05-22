# SARS-CoV-2-Model-with-and-without-temperature-dependence
SARS-CoV-2 Model 

Stochastic model COVID-19
stochastic version of the NHB and RINP compartmental model
Odds are deducted from rates

Populations:

	S: Susceptible
	
	E: Exposed
	
	I: Infected
	
	Q: Infected Reported
	
	R: Recovered
	
	D: Dead
	
	P: Protected Population
	
	Y: Total cases
	
Parameters:

	-alpha: confination rate
	
	-beta: infection rate
	
	-gamma: incubation rate
	
	-delta: detected infected
	
	-Lambda: recovery rate
	
	-kappa: death rate
	
	-tau: deconfination or realese rate 

The models are separated into three different directories.
  - No Temperature: This directory contains the stochastic model that models the dynamics of COVID-19 without taking into account the effect of temperature.

  - Temperature: This directory contains the stochastic model that models the dynamics of COVID-19 taking into account the effect of temperature.

  - Seasonal: This directory contains the stochastic model that models the dynamics of COVID-19 without taking into account the effect of temperature but modeling the infection rate as a function of time t through a sinusoidal function of the annual period.

Each of the models can be run to 4 different locations: Barcelona, Catalonia, Lombardia and Thüringen. These locations can be selected by uncommenting the corresponding part of the code.



