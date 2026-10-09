
library(rotl)
library(treeplyr)
library(ape)
library(dplyr)
library(tibble)

#get the Pantheria dataset
getwd()
setwd("~/Desktop/Current Research/Batsfile")

#Here I am importing bats file with their lifehistory data that  I will match/join 
file1 <- read.csv("AllPanTHERIA.txt", head = TRUE, sep ="\t",stringsAsFactors=F)
#View(file1)

#subset for the chiropteras
columnames <- c(
  "MSW05_Binomial","X5.1_AdultBodyMass_g",
                 "C","X5.2_BasalMetRateMass_g","X17.1_MaxLongevity_m","X21.1_PopulationDensity_n.km2")

file2<-dplyr::select(file1,MSW05_Order,
  MSW05_Binomial,X5.1_AdultBodyMass_g,X5.2_BasalMetRateMass_g,X17.1_MaxLongevity_m,X21.1_PopulationDensity_n.km2,X18.1_BasalMetRate_mLO2hr)

#View(file2)


#getting the columnames
colnames(file2)
columnames<-colnames(file2)
#View(columnames)



#dplyr way
chiroptera1 <- filter(file2, MSW05_Order=="Chiroptera") %>% dplyr::select(., columnames)
#View(chiroptera1)

#rename first the colnames
rename(chiroptera1, species=MSW05_Binomial)


Chiroptera3<-rename(chiroptera1, species=MSW05_Binomial,mass=X5.1_AdultBodyMass_g,BMR.R=X18.1_BasalMetRate_mLO2hr)
#View(Chiroptera3)

Chiroptera4<-rename(Chiroptera3,lifespan=X17.1_MaxLongevity_m,Popdensity=X21.1_PopulationDensity_n.km2,BMR.Mass=X5.2_BasalMetRateMass_g)
#View(Chiroptera4)


#replace -999 with NA
#na_if(x, y) where x is the vector and y is the value to be replaced

Chiroptera5 <- Chiroptera4 %>%
  mutate_if(is.numeric, ~ na_if(., -999))

write.csv(Chiroptera5, "Chiroptera5.csv", row.names = FALSE)
Chiroptera5<- read.csv("Chiroptera5.csv")
#View(Chiroptera5)

#fungal load data from Dr. Hoyt 
#getwd()
#setwd("~/Desktop/Current Research/Batsfile")
batsHoyt = read.csv("2mergedbatsetSheet1.csv")
#View(batsHoyt)

#match the two datasets by "species"

Joinedbats1<-left_join(batsHoyt,Chiroptera5,by="species")
#View(Joinedbats1)

#combine the 2 mass and (by loop! ) log transform the values
Joinedbats1$logmerge_mass <- NA
for(i in 1:nrow(Joinedbats1)){
  Joinedbats1$logmerge_mass[i] <- log(mean(c(Joinedbats1$mean_mass[i], Joinedbats1$mass[i]), na.rm=TRUE))
}

#log transform values first before joining
Joinedbats1$logcluster=log( Joinedbats1$mean_cluster)
Joinedbats1$logforearm=log( Joinedbats1$forearm_length..mm.)
Joinedbats1$logpopsize=log( Joinedbats1$mean_popsize)
Joinedbats1$logtouching=log( Joinedbats1$mean_touching+0.01)
#Joinedbats1$logtemp=log( Joinedbats1$mean_temp + 273.15)
Joinedbats1$logfungalloads=(Joinedbats1$fungal_loads)
colnames(Joinedbats1)
#View(Joinedbats1)
#has the additional log transformed columns
saveRDS(Joinedbats1, file = "batdata.RDS")
batdata <- readRDS("batdata.RDS")

batdata$log_load_sd=log(batdata$loads_sd)

batdata$logmerge_mass
Joinedbats1$mean_mass
Joinedbats1$mass
#tree
battree <- rotl::get_study(study_id='ot_254')

#match the treetips to your data names then 'join' tree with data
battree$tip.label <- gsub("_", " ", battree$tip.label)

#Visualizinng the tree
#plot(battree,no.margin=TRUE,edge.width=2)
#plot(battree,no.margin=TRUE,edge.width=1)

#save as csv (to add location) and then RDS
fungals71<-
  batdata%>% 
  dplyr::select(species,disease._state,fungal_loads, loads_sd, log_load_sd, log_c_loads,logcluster,logmerge_mass,logforearm,logpopsize,logtouching,mean_temp)
write.csv(fungals71, file = "fungals71.csv")
fungals71 <- read.csv("fungals71.csv")
View(fungals71)

fungals71<-
  batdata%>% 
  dplyr::select(species,disease._state,fungal_loads, loads_sd, log_load_sd,log_c_loads,logcluster,logmerge_mass,logforearm,logpopsize,logtouching,mean_temp)
saveRDS(fungals71, file = "fungals71.RDS")
fungal71<-readRDS("fungals71.RDS")
View(fungal71)

#make tree data and then add location

#td71
td63 <- make.treedata(battree, fungals71)
#View(td63$dat)

#check which species where dropped
fungals71$species == td63$dat$species
not_in_td63 <- fungals71[!fungals71$species %in% td63$dat$species, ]
View(not_in_td63)


pred_ranges <- read.csv("predictions_ranges.csv",header=TRUE)
View(pred_ranges)

# Combine trait data with location
td63$dat <- data.frame(td63$dat, Location = pred_ranges$Location, group=pred_ranges$Group)

# Set species names to match the phylogeny tips
td63$dat <- data.frame("species"=td63$phy$tip.label, td63$dat)
#row.names(td63$dat) <- td63$phy$tip.label

td63$dat <- td63$dat[,-c(2, 12)]
View(td63$dat)

saveRDS(td63, file = "td63.RDS")
td63 <- readRDS("td63.RDS")
View(td63$dat)
View(td63$phy$tip.label)

plot(td63$phy, show.tip.label = FALSE)
?plot


#### for species with log_c_loads data (27)
df_clean <- td63$dat[!is.na(td63$dat$log_c_loads), ]
species_to_drop <- setdiff(td63$phy$tip.label, df_clean$species) #gives you species in the tree that are not in df_clean
tree_clean <- drop.tip(td63$phy, species_to_drop)
td27 <- make.treedata(tree_clean, df_clean)

td27$dat <- data.frame("species"=td27$phy$tip.label, td27$dat)
View(td27$dat)

saveRDS(td27,file="td27.RDS")
td <- readRDS("td27.RDS") #or
td27 <- readRDS("td27.RDS")
View(td$dat) #27


td$dat <- data.frame("species"=td$phy$tip.label, td$dat)
#check 
all(td$dat$species == td$phy$tip.label)  # Should return TRUE

