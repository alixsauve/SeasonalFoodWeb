To produce input tables for simulating predator-prey dynamics in the Bia&#322;owie&#380;a forest with the seasonal food web model described in our manuscript, the following R scripts are run with the following order:

* To prepare collected data:

1. `SeasonalDensity.R` converts collected density data in numbers of individuals per square kilometers (N/km<sup>2</sup>). Output density table is written in `../../Data`.

2. `ConvertRawDietPref.R` reads predator dietary preferences as described in J&J1998, and recalculate preferences when considering only the trophic interactions described in SEtAl2016's work.

* **Intra-specific parameters**

6. `IntrinsicGrowthRates.R` estimates the intrinsic growth rates of prey species from various groups (voles and mices, amphibians, shrews, moles, birds, reptiles, ungulates) assuming an allometric scaling of this parameter.

7. `IntraSpCompRates.R` estimates the intra-specific competition rate of prey species from various groups (voles and mices, amphibians, shrews, moles, birds, reptiles, ungulates) based on the greatest recorded densities and their intrinsic growth rates.

8. `PredBaseMortRates.R` estimates the baseline mortality rates of predator species based on their maximum life expectancies.

9. `PredBirthRates.R` calculates the birth rate of each predator species based on the population increase during the breeding season.

10. `DDMortRates.R` estimates the density-dependent mortality rates of predator species with two methods: one is based on the population growth during the breeding season while the other is based on the maximum intake of food one individual predator have each year.

* **Inter-specific parameters**

2. `ImpactOnPrey.R` calculates the expected amount of prey consummed by each predator per season (*G<sub>ki</sub><sup>S,W</sup>* in our manuscript).

3. `ImpPredImpact.R` explores the magnitude of each trophic link. This is not directly used in our manuscript.

4. `EstimatingDiscoveryRates_TypeI.R` estimates the predator discovery rates across seasons, under the hypothesis of a type I functional response of the predators.

5. `EstimatingDiscoveryRates_TypeII.R` estimates the predator discovery rates across seasons, under the hypothesis of a type II functional response of the predators.


* **Final verification**

11. `VerifyInputData.R` reads the tables describing the food web parameterisation and verify that the modelled species and the interaction lists match.

