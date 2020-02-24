# In this script, we calculate the birth rate of each predator species

rm(list=ls(all=TRUE))

library(RColorBrewer)

ALPHA <- 0.75 # scaling parameter (3/4 power law)

# define working directory
DAT_DIR <- "../../Data/BPF_Data"
WK_DIR <- "../../Data/ProcessedData"

# load data
setwd(DAT_DIR)
PostBreedingDens <- read.csv("PredPostBreedingBiomass.csv", header = TRUE, stringsAsFactors = FALSE)
PredInfo <- read.csv("Table5-1.csv", header = TRUE, dec = ",", stringsAsFactors = FALSE)

setwd(WK_DIR)
SpQuantities <- read.csv("SpDensitySeasonsFinal.csv", header = TRUE, stringsAsFactors = FALSE)
PreyList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Prey"])
PredList <- sort(SpQuantities$Taxon[SpQuantities$Trophic_level == "Predator"])

# verify order of species
SpQuantities <- SpQuantities[match(c(PredList, PreyList), SpQuantities$Taxon),]
PredQuantities <- subset(SpQuantities, Trophic_level == "Predator")

# create and initialise a dataset for predators birth rates
BirthRatesPred <- subset(PredQuantities, select = c(Taxon, Clade, Trophic_level))
BirthRatesPred$Clade[BirthRatesPred$Clade == "Raptor"] <- "Bird"
BirthRatesPred$Clade[BirthRatesPred$Clade != "Bird"] <- "Mammal"
BirthRatesPred$B <- NA; BirthRatesPred$SourceB <- NA # birth rates, i.e. number of births per individual per year

# predators' intrinsic growth rates (ie. bi-mi)
BirthRatesPred$R[match(BirthRatesPred$Taxon, PostBreedingDens$Taxon)] <- 365/168*log(PostBreedingDens$PostBreedingBiomass_gha[match(BirthRatesPred$Taxon, PostBreedingDens$Taxon)]/SpQuantities$Mean_density[match(BirthRatesPred$Taxon, SpQuantities$Taxon)])
BirthRatesPred$SourceR <- "Table3.3"

# birth rates
# Accipiter gentilis
# breeding pairs raised 1.1 young per year on average
BirthRatesPred$B[BirthRatesPred$Taxon == "Accipiter_gentilis"] <- 1.1/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Accipiter_gentilis"] <- "p131"

# Accipiter nisus
# up to 77% of breeding pairs raised 1.2 young per year
BirthRatesPred$B[BirthRatesPred$Taxon == "Accipiter_nisus"] <- 1.2/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Accipiter_nisus"] <- "p134"

# Aegolius funereus
# p. 155 "In seven observed broods, from three to five owlets (on average 4.3) were recorded (Domaszewicz 1993)."
BirthRatesPred$B[BirthRatesPred$Taxon == "Aegolius_funereus"] <- 4.3/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Aegolius_funereus"] <- "p155"

# Aquila pomarina
# p. 141: "In 1985-1994, from 27 to 75% of pairs (on average 52%) ended their breeding season
# with a success, but fratricide meant that, except for one case, all successful
# pairs produced only one young (Fig. 3.48)."
BirthRatesPred$B[BirthRatesPred$Taxon == "Aquila_pomarina"] <- 0.52/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Aquila_pomarina"] <- "p141"

# Asio_otus
# p. 157: "In seven broods, three to four juveniles (on average 3.3) were seen (Domaszewicz 1993)."
BirthRatesPred$B[BirthRatesPred$Taxon == "Asio_otus"] <- 3.3/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Asio_otus"] <- "p157"

# Buteo buteo
# p. 137: "Of the mean number of 2.4 eggs laid per pair, only 0.7 young fledged"
BirthRatesPred$B[BirthRatesPred$Taxon == "Buteo_buteo"] <- 0.7/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Buteo_buteo"] <- "p137"

# Canis lupus
# each pack (4-5 individuals) produces 0.9 pups per year on average
BirthRatesPred$B[BirthRatesPred$Taxon == "Canis_lupus"] <- 0.9/4.5
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Canis_lupus"] <- "p108"

# Falco subbuteo
# p. 147: "On average, 0.9 juveniles fledged per breeding pair in the Polish part (Pugacewicz 1996)"
BirthRatesPred$B[BirthRatesPred$Taxon == "Falco_subbuteo"] <-  0.9/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Falco_subbuteo"] <- "p147"

# Glaucidium_passerinum
# p. 154: "A few of our own observations of parents feeding
# their fledged owlets indicated that, on average, three young per successful brood were reared."
BirthRatesPred$B[BirthRatesPred$Taxon == "Glaucidium_passerinum"] <- 3/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Glaucidium_passerinum"] <- "p154"

# Lutra lutra
# p. 129: "Fifteen observations of litters showed that
# Bialowieza's otters had from one to three pups (on average two)"
# we assume half the population reproduces.
BirthRatesPred$B[BirthRatesPred$Taxon == "Lutra_lutra"] <- 2/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Lutra_lutra"] <- "p129"

# Lynx lynx
# 23% of females reproduce, and each has 1.6 kitten per year which reach independency
# we assume a 1:1 sex ratio
BirthRatesPred$B[BirthRatesPred$Taxon == "Lynx_lynx"] <- 0.23*1.6/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Lynx_lynx"] <- "p103"

# Martes_martes
# p. 118: "From one to four kittens in a litter (on average 2.6) were seen when they leave nests"
BirthRatesPred$B[BirthRatesPred$Taxon == "Martes_martes"] <- 2.6/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Martes_martes"] <- "p118"

# Meles meles
# p. 117: "Sliwinski (1987), who studied badger population in the deciduous and mixed
# forests of Rog6w (central Poland), found that from one to ten badgers lived
# in one sett in summer, on average 3.2 individuals per active sett.
# This number included a mean of 0.9 cubs/sett (range from 0 to 6)."
BirthRatesPred$B[BirthRatesPred$Taxon == "Meles_meles"] <- 3.2/(3.2-0.9)
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Meles_meles"] <- "p117"

# Mustela_nivalis
# p. 122: "Weasels reproduce from spring until early autumn. After 35 days of pregnancy, females
# give birth to four to eight young (on average 5.2), and may become pregnant
# again while still nursing their first litters (Jedrzejewska 1987). Thus, an adult
# female may bear two litters in one summer."
BirthRatesPred$B[BirthRatesPred$Taxon == "Mustela_nivalis"] <- 2*5.2/2
# BirthRatesPred$B[BirthRatesPred$Taxon == "Mustela_nivalis"] <- 1.6 # Turchin & Hanski 1997
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Mustela_nivalis"] <- "p122"

# Mustela putorius
# p. 128: "Scarce information on polecat litters indicates that the most common litter
# size is four young per female."
BirthRatesPred$B[BirthRatesPred$Taxon == "Mustela_putorius"] <- 4/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Mustela_putorius"] <- "p128"

# Mustela vison
# p. 131: "No data on mink reproduction in Bialowieza are available. In the forests
# of Belarus, Sidorovich {1993) reported that, in summer, on average 3.3 young followed the female.
BirthRatesPred$B[BirthRatesPred$Taxon == "Mustela_vison"] <- 3.3/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Mustela_vison"] <- "p131"

# Nyctereutes_procyonoides
# p. 115: "The reproductive potential of raccoon dogs is very high; from four to
# eight pups (on average 6.2) were seen in March through July (Bunevich 1983 a; unpubl. data)"
BirthRatesPred$B[BirthRatesPred$Taxon == "Nyctereutes_procyonoides"] <- 6.2/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Nyctereutes_procyonoides"] <- "p115"

# Pernis apivorus
# p. 141: "0.5 juveniles fledged per breeding pair."
BirthRatesPred$B[BirthRatesPred$Taxon == "Pernis_apivorus"] <- 0.5/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Pernis_apivorus"] <- "p141"

# Strix aluco
# p. 151: "If we accept that estimate as typical for all years with moderate numbers of rodents, then in
# 1988-1989, the numbers of owls in late summer were 1.7 times higher than in winter (Fig. 3.56)."
BirthRatesPred$B[BirthRatesPred$Taxon == "Strix_aluco"] <- 1.7-1
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Strix_aluco"] <- "p151"

# Vulpes vulpes
# p. 112: "In most cases (30%) a litter consisted of 3 pups, and the mean was 3.2 pups/litter."
BirthRatesPred$B[BirthRatesPred$Taxon == "Vulpes_vulpes"] <- 3.2/2
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Vulpes_vulpes"] <- "p112"

# # Mustela erminea
BirthRatesPred$B[BirthRatesPred$Taxon == "Mustela_erminea"] <- mean(BirthRatesPred$B[grep("Mustela", BirthRatesPred$Taxon)], na.rm = TRUE)
BirthRatesPred$SourceB[BirthRatesPred$Taxon == "Mustela_erminea"] <- "MeanMustela"

 BirthRatesPred$B <- round(BirthRatesPred$B, 2); BirthRatesPred$R <- round(BirthRatesPred$R, 2)
write.csv(BirthRatesPred, file = "PredatorsBirthRates.csv", row.names = FALSE, quote = FALSE)