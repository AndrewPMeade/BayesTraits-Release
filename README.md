ReadMe for the Fabric-regression in BayesTraits V5.0 Beta – a prerelease version of the program. Not to be distributed.

The manual for BayesTraits which details, how to build and run the program, set parameters and analyse data, can be found in the manual section of the repository.  Note that the manual is not fully up-to-date but along with this README the Fabric-regression model as used in the paper can be run.

The commands below were used to run the Fabric-regression model for the brain and body size dataset used in the paper. 

The first two lines select the regression model and MCMC analysis. 

The next three lines set burn-in to 10M iterations, the chain to run for 20M iterations and to sample the chain every 100K iterations. The sampling frequency is wide to return a low autocorrelation among successive samples.

Line 6 selects the Fabric model, and the next four lines set the priors on the betas, evolvability (variance) scalars and regression coefficients. A threshold value of -2.25 log units is applied to the Fabric betas as a partial prior on entry into the model. Priors were chosen on the basis of our previous work (Pagel, O’Donovan and Meade, Nat Comms 2022).

The next line sets the stepping stone sampler to run using 1000 stones, each for 50K iterations. 

The last two lines specify that the reverse jump component of the model is turned off during the estimation of the marginal likelihood and the model is run. For a more detailed description of the settings, please see the manual. 
<br/>
<br/>
9<br/>
2<br/>
burnin 100000000<br/>
iterations 200000000<br/>
sample 100000<br/>
Fabric<br/>
Prior FabricBeta gamma 3.18 0.20<br/>
Prior VRNode gamma 1.2 5<br/>
Prior Beta-1 uniform 0 2<br/>
Prior Beta-2 normal 0 1<br/>
RJThreshold FabricBeta -2.25<br/>
Stones 1000 50000<br/>
RJLockModel<br/>
run<br/>
