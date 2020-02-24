# Data from J&#281;drzejewska & J&#281;drzejewski (1998)

We provide here data extracted from [J&#281;drzejewska & J&#281;drzejewski's book](https://www.springer.com/gp/book/9783540641384) (1998). Please cite their monograph when using these tables or the resulting parameterisation described in our manuscript.

* `Table5-1.csv` (pp. 329-330) lists the predator species we model and their activity in the Bia&#322;owie&#380;a forest. Column names have the following meanings: `Taxon` for the scientific names; `Clade` is `Raptor` for birds of prey and family name for mammals; `DailyFoodIntake` the amount of food consummed per day by adults (in g/day, c.f., column `Unit_DailyFoodIntake`); `DailyFoodIntake_Juv` is the daily food intake of juveniles; `NdayPresent_W` and `NdayPresent_S` the number of days the species is present on site during winter and summer respectively; `NdayPresent_Juv` is the number of days juveniles are present on site; `TableDiet` is the name of the table describing species diet; `BodyMass` is the mass of adult individuals (in g, c.f., column `Unit_BodyMass`) and comes from Table 3.3 (p. 164).

* Other tables describe the dietary preferences of each predator species:
	- `Table4-3.csv` for *Lynx lynx*;
	- `Table4-8.csv` for *Canis lupus*;
	- `Table4-11.csv` for *Vulpes vulpes*;
	- `Table4-14.csv` for *Nyctereutes procyonoides*;
	- `Table4-15.csv` for *Meles meles*;
	- `Table4-17.csv` for *Martes martes*;
	- `Table4-19.csv` for *Mustela erminea*;
	- `Table4-21.csv` for *Mustela nivalis*;
	- `Table4-23.csv` for *Mustela putorius*;
	- `Table4-25.csv` for *Lutra lutra*;
	- `Table4-27.csv` for *Mustela vison*;
	- `Table4-30.csv` for *Accipiter gentilis*;
	- `Table4-31.csv` for *Accipiter nisus*;
	- `Table4-32.csv` for *Buteo buteo*;
	- `Table4-35.csv` for *Aquila pomarina*;
	- `Table4-40.csv` for *Falco subbuteo*;
	- `Table4-42.csv` for *Strix aluco*;
	- `Table4-45.csv` for *Glaucidium passerinum*;
	- `Table4-46.csv` for *Aegolius funereus*;
	- `Table4-47.csv` for *Asio otus*;
	- `p268.csv` for *Pernis apivorus* which is not originally a table but a paragraph on the share of vertebrate prey in this predator's diet.
Each table gives the common name of the prey species (`Common_name`), their scientific name (`Scientific_name`), and the percentage biomass of the prey species in the predator's diet for both  summer (`FracBioInDiet_W`) and winter (`FracBioInDiet_S`) if the predator is active on site during both seasons.

* `PredPostBreedingBiomass.csv` (partial extraction of Table 3.3, p. 164) describes the density of predator populations after the breeding-season (`PostBreedingDensity_N10km2` in N/10km<sup>2</sup>), their biomass (`PostBreedingBiomass_gha` in g/ha = kg/10km<sup>2</sup>) assuming juveniles weight half the mass of the adults. This table also provides the mean density of each predator population (`MeanDensity_N10km2` in N/10km<sup></sup>) which can be interprated as the density before the breeding season.

* `HypAverageIntrinsicGR.csv` lists the maximum number of juveniles produced per breeding pair in spring-summer (column `NMaxJuv`) for various groups of prey (column `Clade`). These numbers are extracted from Table 2.14 (p. 95). These numbers are used to calculate the mean intrinsic growth rate for these prey types (`MeanR`) with average body mass (`MeanBodyMass` in g) with the R script `../../Parameterisation/Scripts/IntrinsicGrowthRates.R`.

* `IntListSubsets.csv` is produced by `../../Parameterisation/Scripts/ImpPredImpacts.R` and describes the interactions happening in our food web model. Each line corresponds to one interaction with `LowerTaxon` being the prey name (from prey group `LowerClade`) and `UpperTaxon` being the predator name (from predator group `UpperClade`). Columns `SeasonX`, `MonthlyX`, `WeeklyX` and `DailyW` are `TRUE` if one individual prey is killed at least seasonally, monthly, weekly or daily during season `X` (`X` being `S` for summer or `W` for winter).
