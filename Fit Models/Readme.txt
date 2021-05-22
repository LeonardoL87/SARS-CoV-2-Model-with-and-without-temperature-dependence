Stochastic model COVID-19 stochastic version of the NHB and RINP compartmental model Odds are deducted from rates

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
	
Each script can run the fittings for the model with temperature (SEIQRDFitTemperature.py), No temprature  (SEIQRDFitNoTemperature.py) and sesonal (SEIQRDFitSeasonal.py)
By uncommening the lines corresponing to location X and wave Y it's possible to perform the titting for each case.
X={Barcelona, Catalonia, Lombardia, Thüringen}
Y={firts wave, second wave, third wave (if X=Barcelona or Catalonia )}
1
