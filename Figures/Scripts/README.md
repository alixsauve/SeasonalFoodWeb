# Scripts for figures

The following scripts produces the figures displayed in our manuscript on the parameterisation of seasonal food webs. Note that the script to calculate and draw food web metrics can be found in `../../FoodWebMetrics`.

* `SeasonalFoodWeb.R` produces a two-panel figure with the food web of the Bia&#322;owie&#380;a forest during summer and winter.
It uses R packages `bipartite` and `RColorBrewer` (the later not being mandatory). 

* `IntraSpParamFigures.R` creates a multi-panel figure displaying all parameter estimates, but discovery rates and handling times, against species body masses.

* `PhasePortraits.R` draws a four-panel figures with phase portraits and time series for different types of parameterisation (for the predator density-dependent mortality rate) and types of predator functional response.

* `SimuVSObsSeasonDens.R` compares simulated with observed densities for different seasons. It also draws insets to display the distribution of the simulated to observed density ratio.

* `AfterVSBeforeBreedingDens.R` compares species densities after and before the breeding period, both with simulated dynamics and observed data.

* `SeasonalSSEAndSAD.R` compares squared errors between parameterisation methods of the predator density-dependent mortality, as well as species abundance distributions across simulations and seasons.
