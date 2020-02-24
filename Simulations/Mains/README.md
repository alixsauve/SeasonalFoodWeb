# Simulating the dynamics of a seasonal food web

`ParamFW_MainSimu_*.m` simulates the dynamics of a seasonal food web, parameterised with the data collected on the Bia&#322;owie&#380;a forest. Both type I and type II functional responses are considered in each script which only differ by the estimating method for the predator density-dependent mortality rate (`*_DDModel.m` or `*_DDData.m`).

# Testing sensitivity to initial conditions

The following scripts simulate the dynamics of the same seasonal food web model as `ParamFW_MainSimu_*.m` but increasing or decreasing the initial densities of a couple of species by 10 or 50%.

* `VarIC_SmX_GiY_IntroUncert.m` modifies the densities of fish, reptiles, and amphibians separately and altogether.
It simulates the dynamics with a type `X` functional response (either 1 or 2), and with the predator density-dependent mortality either estimated based on predator maximal intake of prey (`VarIC_SmX_GiModel_IntroUncert.m`) or on data-based estimates of predator intrinsic growth rate (`VarIC_SmX_GiData_IntroUncert.m`).


* `VarIC_SmX_GiY_DegraInform.m` modifies the densities of birds and mammals separately, and both at the same time.
It simulates the dynamics with a type `X` functional response (either 1 or 2), and with the predator density-dependent mortality either estimated based on predator maximal intake of prey (`VarIC_SmX_GiModel_DegraInform.m`) or on data-based estimates of predator intrinsic growth rate (`VarIC_SmX_GiData_DegraInform.m`).

# Investigating to what extent model outcomes are sensitive to the quality of empirical abundance data 

The following scripts simulate the dynamics of the same seasonal food web model as `ParamFW_MainSimu_*.m` but altering predator discovery rates assuming some level of sampling error for some prey groups. We consider two scenarios:

* *Introducing uncertainty*: We modify the densities for summer and winter for lesser–known taxa
(i.e., fish, reptiles, and amphibians) by increasing or decreasing their estimates by 10% or 50%
and then recalculate the seasonal discovery rates. We check whether this affects the persistence of
the whole community, as well as the predictions of the dynamical model. This exercise is meant
to better inform on the dynamics and persistence of the actual food web, given uncertainties on
these abundance data.
The simulations for this scenario are produced by MATLAB scripts named as `VarDiscRates_SmX_GiY_IntroUncert.m`.

* *Degrading information*: This time, we alter the densities of the best–studied taxonomic groups
(namely mammals and birds). This is meant to test, in a more general manner, whether degrading
the accuracy of the information we have on their abundances alters discovery rates sufficiently to
affect the dynamical outcomes of the model.
The simulations for this scenario are produced by MATLAB scripts named as `VarDiscRates_SmX_GiY_DegraInform.m`.
