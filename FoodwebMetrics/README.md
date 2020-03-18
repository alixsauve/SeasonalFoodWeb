In this folder, all food web metrics discussed in our manuscript are calculated by `FWMetrics.R`.

`FWMetrics.R` requires the R package `bipartite` (by [Dormann et al.](https://cran.r-project.org/web/packages/bipartite/bipartite.pdf)) and calls other functions:

* `shannonSeason.R` returns our index of predation seasonality which is a Shannon index applied on species diets across seasons. It takes a data frame `EdgeList` describing each interaction (`LowerTaxon` and `UpperTaxon`) and their weight. Other input variables are the columns names (`SummerWeightName` and `WinterWeightName`) corresponding to the weight of each interaction for summer and winter respectively.

* `edgelist2Mat.R` returns a matrix with lower taxa in rows and upper taxa in columns based on a data frame listing each interaction (variable `EdgeList`) with at least two columns `LowerTaxon` and `UpperTaxon`. `ColWeightName` is an optional argument, a character string providing the column name for the interaction weight.

* `dietBrayCurtis.R` returns the Bray-Curtis dissimilarity between two diets. This function takes two vectors of the same size as inputs.

* `jaccardDistDiet.R` returns the distance between two diets based on the Jaccard index.
