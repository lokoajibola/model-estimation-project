Sequential Parameter and Model Estimation for Financial Timeseries MATLAB codes.

The Main Code is divide into 6 sections:

A. SIMULATION TYPE SELECTION
This section requires setting model and filter to suit the required simulation. Commenting out the unrequired condition simply gets the code running. Choices to be made are
1. MODEL: 	'md_RBF' - for radial basis function network 
		'BS' - for Black-Scholes model
2. FILTER:	'SOEKF' - for Second Order EKF
		'EKF' - for Extended Kalman Filter
3. DATA TYPE:	'real' - for generating the FTSE100 data
		'synthetic' - for synthesizing data with steady states
4. OPTIONS CONTRACT TYPE:	'call' - for European call option
				'put' - for European put option
5. DATA YEAR:	real DATA TYPE requires a choice of a year. Years of data 		available are 2012 2013 2014
6. STRIKE PRICE(X):	choice of the strike price for the contract picked eg 					6700. The available prices for the various years are:
		2014 - 6800 or 6700; 2013 - 6100 or 6000; 2012 - 5500 or 5600


B. PARTAMETER INITIALIZATION

Here the value for the prior state(X), covariance(P) and the noise(Q) and (R) are set respectively


C.  STEADY STATE PARAMETER INITIALIZATION

Under synthetic DATA TYPE, it required that you set the states to constant true state vector values depending on the MODEL choice.

For Black-Scholes(BS) model choice, set: 
r - risk-free interest rates divided by 100
sig - volatility

For RBF model choice, set:
mn - representing the mean vector to a row vector

D.  MODEL AND DATA PREPARATION
This section gerates the models and the data.
This section of the code needs no tweeking. It simply runs the functions and computes the results based on the simulation type selection previously made.

E.  KALMAN FILTER ITERATION PROCESS
This section is the loop where the Kalman filter runs and estimations are made.
This section of the code needs no tweeking. It simply runs the functions and computes the results based on the simulation type selection previously made.

F. PERFORMANCE ANALYSIS
This section computes the resulting error and plots for visualization of the performance based on the choice of simulation type made.
This section of the code needs no tweeking. It simply runs the functions and computes the results based on the simulation type selection previously made.
