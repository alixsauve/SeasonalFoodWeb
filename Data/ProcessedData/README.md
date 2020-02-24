The following files correspond to files either resulting from the processing of "raw" data:


* `PerCapitaIntakePreyBiomass.csv` lists the amount of prey biomass consummed by each individual predator according to their daily food intakes and their dietary preferences.
Each line corresponds to the biomass (in g) of prey species `LowerTaxon` (from prey group `LowerClade`) consumed by an individual of the predator population `UpperTaxon` (from predator group `UpperClade`) during winter (`IntakePrey_W`) and summer (`IntakePrey_S`) respectively.
This table is generated with the R script `../../Parameterisation/Scripts/ImpactOnPrey.R`.

* `PredatorsBirthRates.csv` is a table of estimates for predator birth rates (column `B`) and their intrinsic growth rates (column `R`). Columns `SourceB` and `SourceR` provide the source of the input values for these estimates (c.f., Supporting Information B).

* `PredMortality.csv` is a table of predator life expectancy in the wild or/and in captivity (`LifeExpectWild` and `LifeExpectCapti`) according to the data base [*AnAge*](https://genomics.senescence.info/species/), and the resulting estimates of baseline mortality rates (`MbasedonLifeExpectWild` and `MbasedonLifeExpectCapti`).

* `PredPref.csv` is a compilation of the tables listed in `../BPF` and describing predator dietary preferences, with correction on percentage biomasses to have their sum equal to 1 for each predator/season. Only interactions described in our food web model are kept in this list. This table is generated with `../../Parameterisation/Scripts/ConvertRawDietPref.R`.

* `SpDensBiomass.csv` is a table listing species biomass densities during spring (column `InitDensity_Nha`) and their body mass (column `BodyMass_g`). These information are used for the initial conditions of the simulations described in `../../Simulations`.

* `SpDensitySeasons*.csv` are tables of species densities across seasons. The several versions (1, 2 and *Final*) correspond to small edits such as conversion of density units, addition of seasonal densities under specific working hypotheses and subsetting of the table. Details are given in the R scripts executed following the workflow in `../../Parameterisation/Scripts`.
Each line correspond to one species named with its scientific name (`Taxon`). For each, we specified its species group with `Clade` and `CladeBis` which is a simplified version of column `Clade` (e.g., bird family is no longer specified). The `Trophic_level` (*Predator* or *Prey*) is specified, as well as the migratory status of birds (`IsMigratoryBird`). Each column `*_density` gives the species density for a given season with units as in columns `Unit_*_density`.
