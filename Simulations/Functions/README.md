The following functions describe population growth rates (temporal derivates of population densities):

* `SeasonLV_TypeI.m` when predation is ruled by a **type I functional response**.

* `SeasonLV_TypeII.m` when predation is ruled by a **type II functional response**.

These functions take the following arguments as inputs: `tt` for time (a numeric variable), `xt` the species biomass densities (a numeric vector), `SpParam` the table describing parameters of population growth due to intra-specific processes, `IntParam` the table describing each interaction with the mean discovery rate and its seasonal forcing strength, `IsStep` a logical variable that is `TRUE` if the forcing signal is rectangular, and `IsReproSeasonal` another logical variable that is `TRUE` if the prey reproduction is seasonal.
Like most ODE functions, they return the population growth (numeric vector `xdot`) at time `tt`.
