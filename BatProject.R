#get the Pantheria dataset
getwd()
setwd("~/Desktop/Current Research/Batsfile")
install.packages("rotl")
install.packages("remotes")
install.packages("devtools")
devtools::install_github("uyedaj/treeplyr")
install.packages("ape")
library(remotes)
library(devtools)
library(dplyr)
library(tidyr)
library(readxl)
library(dplyr)
library(rotl)
library(treeplyr)
library(ape)

#Here I am importing bats file with their lifehistory data that  I will match/join 
#fungal load data from Dr. Hoyt 
Allpantheria<-read.csv("AllPanTHERIA.txt")
#View(Allpantheria)
#mind you this is all the species. not only bats
file1 <- read.csv("AllPanTHERIA.txt", head = TRUE, sep ="\t",stringsAsFactors=F)
#View(file1)

#this is the data from Zukal 
zukalbats<-read_excel("zukalbats.xlsx")
#View(zukalbats)

#subset for the chiropteras

columnames <- c(
  "MSW05_O5","MSW05_Binomial","X5.1_AdultBodyMass_g","X8.1_AdultForearmLen_mm",
  "C","X5.2_BasalMetRateMass_g","X17.1_MaxLongevity_m","X21.1_PopulationDensity_n.km2")

file2<-file1 %>% dplyr::select(c("MSW05_Order",
                                 "MSW05_Binomial","X5.1_AdultBodyMass_g","X5.2_BasalMetRateMass_g","X8.1_AdultForearmLen_mm","X17.1_MaxLongevity_m","X21.1_PopulationDensity_n.km2","X18.1_BasalMetRate_mLO2hr"))

#View(file2)


#getting the columnames
colnames(file2)
columnames<-colnames(file2)
#View(columnames)

#dplyr way
chiroptera1<-file2[file2$MSW05_Order == 'Chiroptera',]
#View(chiroptera1)


#rename first the colnames
rename(chiroptera1, species=MSW05_Binomial)


Chiroptera3<-rename(chiroptera1, species=MSW05_Binomial,mass=X5.1_AdultBodyMass_g,BMR.R=X18.1_BasalMetRate_mLO2hr)
#View(Chiroptera3)

Chiroptera4<-rename(Chiroptera3,lifespan=X17.1_MaxLongevity_m,Forearm=X8.1_AdultForearmLen_mm,Popdensity=X21.1_PopulationDensity_n.km2,BMR.Mass=X5.2_BasalMetRateMass_g)
#View(Chiroptera4)


#replace -999 with NA
#na_if(x, y) where x is the vector and y is the value to be replaced
Chiroptera5<-replace(Chiroptera4, Chiroptera4==-999,NA)
#View(Chiroptera5)


#bring in the older data set to later match it
batsHoyt = read.csv("2mergedbatsetSheet1.csv")
#View(batsHoyt)

#match the two datasets by "species"

Joinedbats1<-left_join(batsHoyt,Chiroptera5,by="species")
#View(Joinedbats1)

#combine the 2 mass and (by loop! ) log transform the values
Joinedbats1$logmerge_mass <- NA
for(i in 1:nrow(Joinedbats1)){
  Joinedbats1$logmerge_mass[i] <- log(mean(c(Joinedbats1$logmean_mass[i], Joinedbats1$mass[i]), na.rm=TRUE))
}


#View(Joinedbats1)

###

#Match with bat tree
battree <- rotl::get_study(study_id='ot_254')

#match the treetips to your data names then 'join' tree with data
battree$tip.label <- gsub("_", " ", battree$tip.label)

#Visualizinng the tree
plot(battree,no.margin=TRUE,edge.width=2)
plot(battree,no.margin=TRUE,edge.width=1)

plot(battree,no.margin=TRUE,edge.width=0.5)


#log transform values first before joining
Joinedbats1$logcluster=log( Joinedbats1$mean_cluster)
Joinedbats1$logforearm=log( Joinedbats1$forearm_length..mm.)
Joinedbats1$logpopsize=log( Joinedbats1$mean_popsize)
Joinedbats1$logtouching=log( Joinedbats1$mean_touching+0.01)
#Joinedbats1$logtemp=log( Joinedbats1$mean_temp + 273.15)
Joinedbats1$logfungalloads=(Joinedbats1$fungal_loads)
colnames(Joinedbats1)

write.csv(Joinedbats1, "Joinedbats1.csv")

fungals71<-
  Joinedbats1%>% 
  dplyr::select(species,disease._state,log_c_loads,logcluster,logmerge_mass,logforearm,logpopsize,logtouching,mean_temp)
#View(fungals71)

#td71
td63 <- make.treedata(battree, fungals71)
View(fungals71)
View(td63$dat)

#For later use 30 species with available fungal load data
fungals30<-
  fungals71%>% drop_na(log_c_loads) 
View(fungals30) #30

#td30
td27 <- make.treedata(battree, fungals30)
View(fungals30)
View(td27$dat) #27

#####
fungals71$disease._state[64:71] <- NA
COLS = ifelse(fungals71$disease._state=="epiizootic","blue","orange")
#COLS = ifelse(fungals71$disease._state=="endemic","orange","orange")
names(COLS) <- fungals71$species
#plot(td27$phy,tip.color=COLS[td27$phy$tip.label])
plot(td63$phy,tip.color=COLS[td63$phy$tip.label])
contMap(td63$phy, )
#####
View(td27$dat)

#NON PHYLOGENETIC MODEL LME4 WITH ALL THE VARIABLES FROM THE BEST FITTING PHYLOGENETIC MODEL BELOW
library(lme4)
lmmodel = lmer(log_c_loads ~ logmerge_mass +
                          logforearm + (1|logcluster) +
                          (1|logtouching) +(1|mean_temp), data=td27$dat, na.action = na.omit)

#find phylogenetic signal
library(phytools)
#NOT WORKKING
?phylosig
phylosig(td27$phy, td27[['logmerge_mass']], method="K", nsim=1000, test=TRUE)
phylosig(td27$phy, td27[['logmerge_mass']], method="lambda", nsim=1000, test=TRUE)
phylosig(td27$phy, td27[['log_c_loads']], method="lambda", nsim=1000, test=TRUE)

phylosig(td27$phy, log(td27[['logmean_popsize']]), method="lambda", nsim=1000, test=TRUE)
#all seem to have a high phylogenetic signal <0.05

#NOT WORKING?
plot(td27$logmerge_mass, td27$log_c_loads)


#####################
library(Rphylopars)
td27$dat <- data.frame("species"=td27$phy$tip.label, td27$dat)

plm_BM1 <- phylopars.lm(log_c_loads ~ logmerge_mass+mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM1)
.dat1 <- dplyr::select(td27$dat, species, log_c_loads, logmerge_mass, temp)
.dat1[1,2] <- NA

# if you want to predict the values
#?DONT WORK
pred_BM1 <- phylopars(trait_data=.dat1, tree = td27$phy, model="BM")
recon1 <- cbind(.dat1, pred_BM1$anc_recon[1:63,])
recon1

#for every trait
td27$dat <- data.frame("species"=td27$phy$tip.label, td27$dat)

plm_BM1 <- phylopars.lm(log_c_loads ~ logmerge_mass,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM1)

plm_BM2 <- phylopars.lm(log_c_loads ~ logforearm,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM2)

plm_BM3 <- phylopars.lm(log_c_loads ~ logcluster,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM3)

plm_BM4 <- phylopars.lm(log_c_loads ~ logtouching,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM4)

plm_BM5 <- phylopars.lm(log_c_loads ~ logpopsize,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM5)

plm_BM6 <- phylopars.lm(log_c_loads ~ mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM6)


#DOUBLE
plm_BM7 <- phylopars.lm(log_c_loads ~ mean_temp+logcluster,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM7)

plm_BM8 <- phylopars.lm(log_c_loads ~ logmerge_mass+logcluster,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM8)

plm_BM9 <- phylopars.lm(log_c_loads ~ logmerge_mass+mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM9)

plm_BM10 <- phylopars.lm(log_c_loads ~ logforearm+logcluster,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM10)


plm_BM11 <- phylopars.lm(log_c_loads ~ logpopsize+mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM11)

#with 3 variables
plm_BM12 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM12)

plm_BM13 <- phylopars.lm(log_c_loads ~  logmerge_mass + logcluster + mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM13)

#without mass and forearm_length
plm_BM14 <- phylopars.lm(log_c_loads ~  mean_temp + logcluster + logtouching,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM14)

#without touching
plm_BM15 <- phylopars.lm(log_c_loads ~ logmerge_mass + logcluster + logforearm + mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM15)

#without forearm
plm_BM16 <- phylopars.lm(log_c_loads ~  logmerge_mass + logcluster + logtouching + mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM16)

#with all variables
plm_BM17 <- phylopars.lm(log_c_loads ~ logmerge_mass + logcluster + logforearm + logpopsize + logtouching + mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM17)

plm_BM18 <- phylopars.lm(log_c_loads ~ logmerge_mass + logcluster +logtouching+ logforearm + mean_temp,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM18)

#phylogeny only model
plm_BM19 <- phylopars.lm(log_c_loads ~ 1,  tree = td27$phy, trait_data=td27$dat, model = "BM")
summary(plm_BM19)

pp_BM18 <- phylopars(td27$dat[,-c(2,3)], tree = td27$phy)
pp_BM18$anc_recon
pc_BM18 <- princomp(pp_BM18$anc_recon[1:27,])

par(mfrow=c(2,3))
par(mfrow=c(1,1))
sm <- make.simmap(td27$phy, td27[["disease._state"]])
phytools::phylomorphospace(sm, pc_BM18$scores[,1:2], xlim=c(-20, 10), ylim=c(-10,10), xlab='PC 1',ylab='PC 2',node.by.map=TRUE, colors=setNames(c("red", "black"), c("endemic", "epiizootic")),ftype='off')
phytools::phylomorphospace(sm, pc_BM18$scores[,c(1,3)], xlim=c(-15, 10), ylim=c(-10,10), node.by.map=TRUE, colors=setNames(c("blue", "red"), c("endemic", "epiizootic")))
phytools::phylomorphospace(sm, pc_BM18$scores[,2:3], xlim=c(-15, 10), ylim=c(-10,10), node.by.map=TRUE, colors=setNames(c("blue", "red"), c("endemic", "epiizootic")))

?phytools::phylomorphospace
phytools::fancyTree(td27$phy, type="scattergram", X=pc_BM18$scores[,1:4], ftype="off")

biplot(pc_BM18, choices = c(2,3))

autoplot(data.frame(pc_BM18$scores), colour=td27$dat$disease._state, loadings=TRUE)
summary(pc_BM18)

?phylopars
#TO NOTE THAT, ORIGINALLY WE HAD RUN THESE ANALYSIS WITH 71 SPECIES, BUT MOST OF HAD MISSING DATA EXCEPT FOR MASS AND 
#FOREARM. SO I RERUN WITH ONLY THE 27 SPECIES. 

str(plm_BM18$PPE)
asr <- plm_BM18$PPE$anc_recon
#how many spp do you have?
nrow(plm_BM18$original_data) #27 

#this gives you tiny garphs of the variables against each other to show correlation
plot(asr[1:27,'logmerge_mass'], asr[1:27,'log_c_loads'], pch=21, col=1+!is.na(plm_BM18$original_data$log_c_loads))
pairs(asr)
?pch

#A tree with the bar plots to show the actual fungals-re , and the predicted-black
barcols <- 1+!is.na(plm_BM19$original_data$log_c_loads)
plotTree.wBars(td27$phy,asr[1:63,"log_c_loads"]-min(asr[1:63,"log_c_loads"]), col=barcols)

#,args.plotTree=list(ftype="off"),
args.barplot=list(col=barcols,border=barcols,xlab="", xlim=c(0,70)),args.axis=list(at=seq(0,1,by=0.2)))
plot(tree, ftype="off")
legend(x="topright",legend=names(colors),pch=22,pt.cex=2,pt.bg=colors,
       box.col="transparent")

#subset
View(newfungals)



#TRYING FOR A BEST CROSS VALIDATION MODEL

species_with_data <- which(!is.na(td27$dat$log_c_loads))
PRED <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_BM18 <- phylopars.lm(log_c_loads ~  mean_temp+logforearm+logmerge_mass,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  summary(.plm_BM18)
  
  predicted <- .plm_BM18$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM18$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED <- rbind(PRED, preds)
}

View(species_with_data)

###### Non phylogenetic
ltrns <- geiger::rescale(td27$phy, "lambda")
PRED_nonphylo <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_BM18 <- phylopars.lm(log_c_loads ~  logtouching + mean_temp+logforearm+logmerge_mass+logcluster,  tree = ltrns(0), trait_data=.td27$dat, model = "BM")
  summary(.plm_BM18)
  
  predicted <- .plm_BM18$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM18$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_nonphylo <- rbind(PRED_nonphylo, preds)
}

###### Phylo-only
ltrns <- geiger::rescale(td27$phy, "lambda")
PRED_phyloonly <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_BM18 <- phylopars.lm(log_c_loads ~  1,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  summary(.plm_BM18)
  
  predicted <- .plm_BM18$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM18$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_phyloonly <- rbind(PRED_phyloonly, preds)
}


LOOCV <- function(model, lambda=1){
  ltrns <- geiger::rescale(td27$phy, "lambda")
  PRED <- NULL
  for(i in species_with_data){
    .td27 <- td27
    
    refvalue <- .td27$dat$log_c_loads[i]
    .td27$dat$log_c_loads[i]  <- NA
    if(lambda==0){
      .plm_BM18 <- phylopars.lm(model,  tree = ltrns(0), trait_data=.td27$dat, model = "BM")
      
    } else {
      .plm_BM18 <- phylopars.lm(model,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
      
    }
    summary(.plm_BM18)
    
    predicted <- .plm_BM18$PPE$anc_recon[i,'log_c_loads']
    predicted_var <- .plm_BM18$PPE$anc_var[i,'log_c_loads']
    
    preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
    PRED <- rbind(PRED, preds)
  }
  return(list(PRED=PRED, cor=cor(PRED$truevalue, PRED$pred), rmse=mean((PRED$truevalue - PRED$pred)^2)))
}
##Fit all the models that you want to
models <- list(log_c_loads~1, 
               log_c_loads~logmerge_mass,
               log_c_loads~mean_temp, 
               log_c_loads~logmerge_mass + mean_temp,
               log_c_loads~logmerge_mass + logforearm + mean_temp, 
               log_c_loads~logmerge_mass + logcluster + mean_temp,
               log_c_loads~logmerge_mass + logforearm+logcluster+logtouching + mean_temp)
fits_loocv <- lapply(models, LOOCV)
mse <- sapply(fits_loocv, function(x) x$rmse)
cor <- sapply(fits_loocv, function(x) x$cor)
best_model <- which(mse==min(mse))
best_model

nonphylo_loocv <- LOOCV(log_c_loads~logmerge_mass + logforearm + mean_temp, lambda=0)



plot(fits_loocv[[best_model]]$PRED[,2:3], xlab="True Log C-loads", "Predicted Log C-Loads", xlim=c(-10,-3), ylim=c(-10,-3))
lm1 <- lm(fits_loocv[[best_model]]$PRED[,3]~ fits_loocv[[best_model]]$PRED[,2])
abline(lm1, lty=2)

res <- LOOCV(log_c_loads~logforearm, lambda=1)

##Residual Mean Square Error

mean((PRED$truevalue - PRED$pred)^2)
mean((PRED_nonphylo$truevalue - PRED_nonphylo$pred)^2)
mean((PRED_phyloonly$truevalue - PRED_phyloonly$pred)^2)

##  Phenogram of All species predicted under best model
td63$dat <- data.frame("species"=td63$phy$tip.label, td63$dat)
plm_AllBest <- phylopars.lm(log_c_loads~logmerge_mass+logforearm+mean_temp,  tree = td63$phy, trait_data=td63$dat, model = "BM")
phenogram(td63$phy, plm_AllBest$PPE$anc_recon[,'log_c_loads'], fcols=is.na(td63$dat$log_c_loads), fsize=0.3)
nH <- max(branching.times(td63$phy))
upper <- sqrt(plm_AllBest$PPE$anc_var[,'log_c_loads']) + plm_AllBest$PPE$anc_recon[,'log_c_loads']
lower <- -1*sqrt(plm_AllBest$PPE$anc_var[,'log_c_loads']) + plm_AllBest$PPE$anc_recon[,'log_c_loads']


COLS = ifelse(td63$dat$disease._state=="epiizootic","blue","orange")
names(COLS) <- td63$phy$tip.label
#plot(td63$phy,tip.color=COLS[td63$phy$tip.label])
object <- contMap(td63$phy,  plm_AllBest$PPE$anc_recon[1:63,'log_c_loads'], outline = FALSE)

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
## for fun, let's change our contMap gradient
object<-setMap(object,c("purple","blue","green","yellow", "white"))
plot(object,ftype="off",xlim=lastPP$x.lim,fsize=c(0,0.6),
     leg.txt="Log fungal loads",sig=1, outline=FALSE)
## make sure our discrete character is in the order of our tree
status <-setNames(is.na(td63$dat$log_c_loads)+1, td63$phy$tip.label)
cols<-setNames(c("red", "black"),levels(status))
for(i in 1:length(status))
  text(lastPP$xx[i],lastPP$yy[i],object$tree$tip.label[i],
       pos=4,cex=0.6,col=cols[status[i]],font=3)

add.simmap.legend(colors=setNames(cols, c("Known", "Predicted")),prompt=FALSE,
                  x=60,y=-5.4,vertical=FALSE,fsize=0.8)
#tiplabels(td63$phy$tip.label, col=COLS[td63$phy$tip.label])


#gives you a graph to show the correlation between the predicted values and the actual values
par(mar=c(5,4,2,2))
plot(PRED$truevalue, PRED$pred, xlim=c(-10,2), ylim=c(-10,2))
abline(0,1)

corval <- cor.test(PRED$pred, PRED$truevalue)
View(corval)
#corval[["estimate"]]
0.5605851 
#corval[["p.value"]]
0.00235435


# When you remake the plot, you can set 
plot(...., xlab = "True Value", ylab="Predicted Value", xlim=c(-10, 0), ylim=c(-10, 0)) 
#and it will give prettier axis labels than PRED$pred  and cut off the extra 2 units that aren't needed.



#Use the same phylopars lm code, but replace the tree with ltrns(0)
#urooting the tree to 0 to 'make' it a normal linear modal
ltrns <- geiger::rescale(tree30, "lambda")
fit <- phylopars(log_c_loads ~  logtouching + mean_temp+logforearm+logmerge_mass,trait_data=batsnew, tree=ltrns(0.9), pheno_error=FALSE, phylo_correlated=TRUE, pheno_correlated=FALSE, model="lambda")
Summary(fit)

ltrns <- geiger::rescale(td27$phy, "lambda")
l0fit <- phylopars.lm(log_c_loads ~  logmerge_mass + logcluster +logtouching+ logforearm + mean_temp,  tree = ltrns(0), trait_data=td27$dat, model = "BM")
summary(l0fit)
#R-squared: 0.1965868 
#p-value:  0.08450395 


summary(plm_BM18)
#R-squared: 0.3917246 
# p-value: 1.069627e-06 

fit <- phylopars(trait_data=td27$dat, tree=ltrns(0), pheno_error=FALSE, phylo_correlated=TRUE, pheno_correlated=FALSE, model="lambda")
Summary(fit)

## add phylopars command with just variables from PM19 and species name in dat






#Known Vs To be predicted
COLS1 = ifelse(Joinedbats1$log_c_loads=="TRUE","blue","orange")
names(COLS1) <- Joinedbats1$species
plot(pruned.tree,tip.color=COLS1[pruned.tree$tip.label],no.margin=TRUE,edge.width=0.5)




When you remake the plot, you can set 
plot(...., xlab = "True Value", ylab="Predicted Value", xlim=c(-10, 0), ylim=c(-10, 0)) and it will give prettier axis labels than PRED$pred  and cut off the extra 2 units that aren't needed.
 
 