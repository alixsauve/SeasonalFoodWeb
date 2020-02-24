The following tables contains estimates for the seasonal food web model of the Bia&#322;owie&#380;a forest:

* `SpDensBiomass.csv` is a table listing species biomass densities during spring (column `InitDensity_Nha`) and their body mass (column `BodyMass_g`). These information are used for the initial conditions of the simulations described in `../../Simulations`.

* `YearPopParam.csv` is a table describing parameter of population growth due to intra-specific processes: `R` the intrinsic growth rate of prey species; `BETA` the intra-specific competition rate of prey species; `M` the baseline mortality rate of predator species; `DDModel` and `DDData` two estimates of density-dependent mortality rate of predator species. Column `IsCst` is `TRUE` if the resource is assumed to have same biomass for all time *t*.

* `YearIntParam_TypeI*.csv` are tables describing the mean discovery rate (column `G`) of prey `LowerTaxon` by predator `UpperTaxon`. Column `E` gives the strength of the seasonal forcing exerted on the discovery rate, and `SeasonDom` is the season during which the discovery rate is the highest (either summer `S` or winter `W`). This type of table has also been created for food web subsets corresponding to interactions occurring at least daily (`YearIntParam_TypeIDaily.csv`), weekly (`YearIntParam_TypeIMonthly.csv`), monthly (`YearIntParam_TypeIMonthly.csv`) or seasonally (`YearIntParam_TypeISeason.csv`).
These estimates are calculated under the working hypothesis that predation is ruled by a **type I functional response**.

* `YearIntParam_TypeII.csv` is a table similar to `YearIntParam_TypeI.csv` but is generated under the working hypothesis that predation is ruled by a **type II functional response**.
