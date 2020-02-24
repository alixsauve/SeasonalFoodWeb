# In this script, we estimate the intrinsic growth rates of prey species from various groups (voles and mices, amphibians, shrews, moles, birds, reptiles, ungulates).

rm(list=ls(all=TRUE))

library(RColorBrewer)

ALPHA <- 0.75 # scaling parameter (3/4 power law)

# define working directory
DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"
PAR_DIR <- "../../Parameterisation/OutputTables"

# load data
setwd(WK_DIR)
SpQuantities <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)

PreyList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"])
PredList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Predator"])

# verify order of species
SpQuantities <- SpQuantities[match(c(PredList, PreyList), SpQuantities$Taxon),]

# creating and initiating a data frame compilling parameters tuning species growth and mortality
PopParam <- subset(SpQuantities, select = c(Taxon, Clade, Trophic_level, Body_mass))
PopParam$R <- NA; PopParam$BETA <- NA  # i.e. intrinsic growth rate (prey), intraspecific competition (prey)
PopParam$M <- NA; PopParam$DD <- NA # i.e. baseline mortality rate (predator), density-dependent mortality (predator)

# data for prey groups
IntrinsicGR_Clades <- data.frame(Clade = c("VolesMice", "Amphibians", "Shrews", "Birds", "Reptiles", "Fish"),
                                 NMaxJuv = c(250, 160, 30, 15, 30, NA)) # NMaxJuv is the maximum number of juveniles produced per breeding pair in spring-summer
IntrinsicGR_Clades$MeanR <- log((IntrinsicGR_Clades$NMaxJuv + 2) / 2)
IntrinsicGR_Clades$Slope <- NA

# Fish mean intrinsic growth rate
IntrinsicGR_Clades$MeanR[IntrinsicGR_Clades$Clade == "Fish"] <- 0.27*15

# a list of clades, for each we provide a data frame with taxa and index in PopParam
Clades <- vector("list", nrow(IntrinsicGR_Clades)); names(Clades) <- IntrinsicGR_Clades$Clade
# Voles/Mice
VolesMice <- which((PopParam$Clade == "Rodent") & (PopParam$Trophic_level == "Prey"))
VolesMice <- data.frame(Taxon = PopParam$Taxon[VolesMice], Index = VolesMice, stringsAsFactors = FALSE)
Clades$VolesMice <- VolesMice
# Amphibians
Amphibians <- which((PopParam$Clade == "Amphibian") & (PopParam$Trophic_level == "Prey"))
Amphibians <- data.frame(Taxon = PopParam$Taxon[Amphibians], Index = Amphibians, stringsAsFactors = FALSE)
Clades$Amphibians <- Amphibians
# Shrews
Shrews <- which((PopParam$Clade == "Shrew/Mole") & (PopParam$Taxon != "Talpa_europaea") & (PopParam$Trophic_level == "Prey"))
Shrews <- data.frame(Taxon = PopParam$Taxon[Shrews], Index = Shrews, stringsAsFactors = FALSE)
Clades$Shrews <- Shrews
# Birds
Birds <- c(grep("forme", PopParam$Clade), 
           which(((PopParam$Clade == "Bird") | (PopParam$Clade == "Passerine")) & (PopParam$Trophic_level == "Prey")))
Birds <- data.frame(Taxon = PopParam$Taxon[Birds], Index = Birds, stringsAsFactors = FALSE)
Clades$Birds <- Birds
# Reptiles
Reptiles <- which((PopParam$Clade == "Reptile") & (PopParam$Trophic_level == "Prey"))
Reptiles <- data.frame(Taxon = PopParam$Taxon[Reptiles], Index = Reptiles, stringsAsFactors = FALSE)
Clades$Reptiles <- Reptiles
# Fish
Fish <- which((PopParam$Clade == "Fish") & (PopParam$Trophic_level == "Prey"))
Fish <- data.frame(Taxon = PopParam$Taxon[Fish], Index = Fish, stringsAsFactors = FALSE)
Clades$Fish <- Fish

# for each prey group listed above, we estimate the parameter of body mass scaling, assuming the quater power law
# and then estimate r_k values for taxa on this basis
for (SpGroup in names(Clades)){
  MeanBodyMass <- mean(SpQuantities$Body_mass[!is.na(match(SpQuantities$Taxon, Clades[[SpGroup]]$Taxon))], na.rm = TRUE)
  if (SpGroup == "Fish"){
    MeanBodyMass <- 7.9
  }
  Slope <- log(IntrinsicGR_Clades$MeanR[IntrinsicGR_Clades$Clade == SpGroup]) - (ALPHA - 1) * log(MeanBodyMass)
  PopParam$R[Clades[[SpGroup]]$Index] <- round(exp(Slope) * SpQuantities$Body_mass[Clades[[SpGroup]]$Index] ^ (ALPHA - 1), 5)
  
  # save the reference values for plotting
  IntrinsicGR_Clades$MeanBodyMass[IntrinsicGR_Clades$Clade == SpGroup] <- MeanBodyMass
  IntrinsicGR_Clades$Slope[IntrinsicGR_Clades$Clade == SpGroup] <- Slope
}

# for ungulates, cf. Table 2.1 in Jedrzejewska and Jedrzejewski, 1998
PopParam$R[PopParam$Taxon == "Bison_bonasus"] <- 0.17
PopParam$R[PopParam$Taxon == "Capreolus_capreolus"] <- 0.66
PopParam$R[PopParam$Taxon == "Cervus_elaphus"] <- 0.43
PopParam$R[PopParam$Taxon == "Sus_scrofa"] <- 0.93

# for moles, cf. Table 2.14 in Jedrzejewska and Jedrzejewski, 1998
PopParam$R[PopParam$Taxon == "Talpa_europaea"] <- round(log((10+2)/2), 5) # (N_Juv)_max = 10/litter; 1 litter/female/y; N0 = 2; N1 = (N_Juv)_max; N1 = N0 e^r_k

# for bats, we take half the birth mass for Vespertilionidae family (Kunz and Fenton, 2005)
# In the calculation, we assume a 1:1 sex ratio in bats populations.
PopParam$R[PopParam$Taxon == "Undet._bats"] <- log((2 + 1.2*1.2)/2)
PopParam$R[PopParam$Taxon == "Undet._bats"] <- round(PopParam$R[PopParam$Taxon == "Undet._bats"], 5)

# for crayfish, we use observations from a study on eggs incubation (Policar et al, 2004).
# In outdoor conditions, they observed up to 137 stage 2 juveniles per female, and a survival of 47.2% at the end of the incubation period.
# In the calculation, we assume a 1:1 sex ratio in Astacus population.
PopParam$R[PopParam$Taxon == "Astacus_sp."] <- round(log((47.2*137/100 + 2) / 2), 5)

# for the brown hare (Lepus europaeus), we use the greatest number of juveniles (up to 14.5) produced per female per year
PopParam$R[PopParam$Taxon == "Lepus_europaeus"] <- round(log((14.5 + 2) / 2), 5)

# draw a plot
# pick a colour palette for the different clades
CladeCol <- brewer.pal(nrow(IntrinsicGR_Clades) + 3, "Set1"); names(CladeCol) <- c(as.character(IntrinsicGR_Clades$Clade), "Other mammals", "Crustaceans", "Ungulates")
names(CladeCol)[1] <- "Voles & mice"
# assign a colour to each species
PopParam$SpCol <- NA
for (SpGroup in 1:length(Clades)){
  PopParam$SpCol[Clades[[SpGroup]]$Index] <- CladeCol[SpGroup]
}
PopParam$SpCol[(PopParam$Taxon == "Lepus_europaeus") | (PopParam$Taxon == "Undet._bats") | (PopParam$Taxon == "Talpa_europaea")] <- CladeCol[SpGroup + 1]
PopParam$SpCol[PopParam$Taxon == "Astacus_sp."] <- CladeCol[SpGroup + 2]
PopParam$SpCol[PopParam$Clade == "Ungulate"] <- CladeCol[SpGroup + 3]

layout(matrix(c(1, 1, 1, 1, 2), 1, 5))
# the actual plot
par(mar = c(5, 4, 1, 1) + 0.1)
plot(log(PopParam$Body_mass), log(PopParam$R), pch = 21, bg = PopParam$SpCol,
     xlab = "Body mass (g)", ylab = "Intrinsic growth rate (/year)",
     xaxt = "n", yaxt = "n")
axis(1, at = log(10^seq(1, 6)), labels = 10^seq(1, 6))
axis(2, at = log(c(0.10, 1, 5, 10, 15)), labels = c(0.10, 1, 5, 10, 15))
# now add segments to highlight body mass scaling and reference mean species
for (SpGroup in names(Clades)){
  # we draw AB segments
  RangeBodyMass <- range(SpQuantities$Body_mass[Clades[[SpGroup]]$Index], na.rm = TRUE)
  A <- log(RangeBodyMass) + c(-1, 1) * 0.5 # coordinates of point A
  RangeIGR <- range(PopParam$R[Clades[[SpGroup]]$Index], na.rm = TRUE)
  B <- rev(log(RangeIGR)) + c(-1, 1) * (ALPHA-1) * 0.5 # coordinates of point B
  segments(A[1], B[1], A[2], B[2], lty = 2, col = CladeCol[which(names(Clades) == SpGroup)])
  points(log(IntrinsicGR_Clades$MeanBodyMass[IntrinsicGR_Clades$Clade == SpGroup]),
         log(IntrinsicGR_Clades$MeanR[IntrinsicGR_Clades$Clade == SpGroup]),
         pch = 24, col = CladeCol[which(names(Clades) == SpGroup)], bg = "black", cex = 1.5)
}

# the legend
par(mar = c(5, 1, 1, 1) + 0.1)
plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA)
legend("topleft", legend = names(CladeCol), pch = 21, pt.bg = CladeCol, bty = "n")
legend("left", legend = c("Hypothetical average", "Estimates"), pch = c(24, 1), pt.bg = "black", bty = "n", pt.cex = c(1.5, 1))

setwd(PAR_DIR)
colnames(PopParam) <- c("Taxon", "Clade", "TrophicLevel", "BodyMass", "R", "BETA", "M", "DD", "SpCol")
write.csv(PopParam, file = "YearPopParam.csv", quote = FALSE, row.names = FALSE)
setwd(DAT_DIR)
write.csv(IntrinsicGR_Clades, file = "HypAverageIntrinsicGR.csv", quote = FALSE, row.names = FALSE)
