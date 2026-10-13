get the dataset
```{r}

getwd()
setwd("~/Desktop/Current Research/Batsfile")

td <- readRDS("td27.RDS") #or
td27 <- readRDS("td27.RDS")
View(td27$dat)

```


with best model with the outgroup(model 27) no phylogeny
```{r}

#td27$dat$Group <- ifelse(td$dat$group == "ingroup", 1, 0)
ltrns <- function(lambda) geiger::lambdaTree(td$phy, lambda)
plm_group26 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + Group,  tree = ltrns(0), trait_data=td$dat, model = "BM")
summary(plm_group26)


#run from here
species_with_data <- which(!is.na(td27$dat$log_c_loads))
PREDgroup <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_group26 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + Group,  tree = ltrns(0), trait_data=.td27$dat, model = "BM")
  summary(.plm_group26)
  
  predicted <- .plm_group26$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_group26$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PREDgroup <- rbind(PREDgroup, preds)
}

write.csv(PREDgroup ,"PREDgroup.csv",row.names = FALSE)
PREDgroup <- read.csv("PREDgroup.csv")
View(PREDgroup)

#comparison....A lower RMSE means better predictive accuracy.
rmsegroup <- sqrt(mean((PREDgroup$truevalue - PREDgroup$pred)^2))
#1.068426

#plots
plot(PREDgroup$truevalue,PREDgroup$pred, main="group model non phy", xlab="True values", ylab="Predicted values")
model <- lm( PREDgroup$pred~ PREDgroup$truevalue)
abline(model, col = "blue", lwd = 2)

##
sst <- sum((PREDgroup$truevalue - mean(PREDgroup$truevalue))^2)
sse <- sum((PREDgroup$truevalue - PREDgroup$pred)^2)
r2 <- 1 - (sse / sst)
# 0.3470824

```
with the outgroup(model 28) plus phylogeny
```{r}

td27$dat$Group <- ifelse(td$dat$group == "ingroup", 1, 0)
plm_group27 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + Group,  tree = .td27$phy, trait_data=td$dat, model = "BM")
summary(plm_group27)


#here
species_with_data <- which(!is.na(td27$dat$log_c_loads))
PREDgroupphy <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_group27 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + Group,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  summary(.plm_group27)
  
  predicted <- .plm_group27$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_group27$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PREDgroupphy <- rbind(PREDgroupphy, preds)
}

write.csv(PREDgroupphy ,"PREDgroupphy.csv",row.names = FALSE)
PREDgroupphy <- read.csv("PREDgroupphy.csv")
View(PREDgroupphy)


#comparison
rmsegroupphy <- sqrt(mean((PREDgroupphy$truevalue - PREDgroupphy$pred)^2))
# 0.9826101

#plots
plot(PREDgroupphy$truevalue,PREDgroupphy$pred, main="group model phy", xlab="True values", ylab="Predicted values")
model <- lm( PREDgroupphy$pred~ PREDgroupphy$truevalue)
abline(model, col = "red", lwd = 2)

##
sst <- sum((PREDgroupphy$truevalue - mean(PREDgroupphy$truevalue))^2)
sse <- sum((PREDgroupphy$truevalue - PREDgroupphy$pred)^2)
r2 <- 1 - (sse / sst)
# 0.4477546

```

 model 13- mean_temp+logforearm+logmerge_mass

```{r}
species_with_data <- which(!is.na(td27$dat$log_c_loads))
PRED <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_BM12 <- phylopars.lm(log_c_loads ~  mean_temp+logforearm+logmerge_mass,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  summary(.plm_BM12)
  
  predicted <- .plm_BM12$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM12$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED <- rbind(PRED, preds)
}

write.csv(PRED ,"PRED1.csv",row.names = FALSE)
PRED1 <- read.csv("PRED1.csv")
View(PRED1)
##Quantiles



quants <- pnorm(q = PRED1[,2], mean = PRED1[,3], sd = sqrt(PRED1[,4]))
hist(quants, breaks=10)
abline(h=length(quants)/10, lty=2)
ks.test(quants, punif, 0, 1) #Kolmogorov-Smirnov test. Tests whether a distribution is significantly different from another distribution. Testing whether our true data quantiles fall uniformly within our prediction interval by seeing if their quantile distribution is signficantly different from a uniform distribution between 0 and 1. 

#for(i in 1:nrow(PRED1)){
#  quants[i] <- pnorm()
#}

#comparison
rmse <- sqrt(mean((PRED1$truevalue - PRED1$pred)^2))
# 0.9520708

#plots

pdf("./LOOCV_fig.pdf", height=5, width=5)
plot(PRED1$truevalue,PRED1$pred, main="True Values and Predicted Values", xlab="True values", ylab="Predicted values")
model <- lm( PRED1$pred~ PRED1$truevalue)
abline(model, col = "maroon", lwd = 2)

# Add ideal y = x line (perfect predictions)
abline(a = 0, b = 1, col = "blue", lwd = 2, lty = 2)

dev.off()
system ("open LOOCV_fig.pdf")

##
sst <- sum((PRED1$truevalue - mean(PRED1$truevalue))^2)
sse <- sum((PRED1$truevalue - PRED1$pred)^2)
r2 <- 1 - (sse / sst)
#0.4815484

```

LOOCV using 4 traits Model 18..Full with Cluster

```{r}
species_with_data <- which(!is.na(td27$dat$log_c_loads))
PRED18 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_BM18 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + logcluster,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  summary(.plm_BM18)
  
  predicted <- .plm_BM18$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM18$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED18 <- rbind(PRED18, preds)
}

write.csv(PRED18, "PRED18.csv",row.names = FALSE)
pred18<- read.csv("PRED18.csv")
View(pred18)

rmse18 <- sqrt(mean((PRED18$truevalue - PRED18$pred)^2)) #1.055708

```


Model 21

```{r}
ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_21 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM21 <- phylopars.lm(log_c_loads ~ logmerge_mass + logforearm + mean_temp + logcluster + logtouching,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
summary(.plm_BM21)
  
  predicted <- .plm_BM21$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM21$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_21 <- rbind(PRED_21, preds)
}
View(PRED_21)

write.csv(PRED_21, "PRED_21.csv",row.names = FALSE)
PRED_21<- read.csv("PRED_21.csv")
View(PRED_21)

rmse21 <- sqrt(mean((PRED_21$truevalue - PRED_21$pred)^2)) # 1.237356


```

model 7
```{r}
ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_7 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM7 <- phylopars.lm(log_c_loads ~ logmerge_mass+logforearm,  tree = .td27$phy, trait_data=.td27$dat, model= "BM")
  
  predicted <- .plm_BM7$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM7$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_7 <- rbind(PRED_7, preds)
}
View(PRED_7)

write.csv(PRED_7, "PRED_7.csv",row.names = FALSE)
PRED_7<- read.csv("PRED_7.csv")
View(PRED_7)

rmse7 <- sqrt(mean((PRED_7$truevalue - PRED_7$pred)^2)) #0.9661093
```

model 26 Endemic model
```{r}
ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_26 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM26 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp+ endemism,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  
  predicted <- .plm_BM26$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM26$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_26 <- rbind(PRED_26, preds)
}
View(PRED_26)

write.csv(PRED_26, "PRED_26.csv",row.names = FALSE)
PRED_26<- read.csv("PRED_26.csv")
View(PRED_26)

rmse26 <- sqrt(mean((PRED_26$truevalue - PRED_26$pred)^2)) #1.031838
```
Model 23-5 traits
```{r}

ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_23 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM23 <- phylopars.lm(log_c_loads ~   logmerge_mass + logforearm + mean_temp +logcluster + logtouching,  tree = ltrns(0), trait_data=.td27$dat, model = "BM")
  
  predicted <- .plm_BM23$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM23$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_23 <- rbind(PRED_23, preds)
}
View(PRED_23)

write.csv(PRED_23, "PRED_23.csv",row.names = FALSE)
PRED_23<- read.csv("PRED_23.csv")
View(PRED_23)

rmse23 <- sqrt(mean((PRED_23$truevalue - PRED_23$pred)^2)) #1.617258
```


Model 25 with Location
```{r}

ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_25 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM25 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp+ loc,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  
  predicted <- .plm_BM25$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM25$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_25 <- rbind(PRED_25, preds)
}
View(PRED_25)

write.csv(PRED_25, "PRED_25.csv",row.names = FALSE)
PRED_25<- read.csv("PRED_25.csv")
View(PRED_25)

rmse25 <- sqrt(mean((PRED_25$truevalue - PRED_25$pred)^2)) # 1.020179
```


Model 20...All 6 Variables
```{r}

ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_20 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM20 <- phylopars.lm(log_c_loads ~ logmerge_mass + logforearm + mean_temp + logpopsize + logcluster + logtouching ,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  
  predicted <- .plm_BM20$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM20$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_20 <- rbind(PRED_20, preds)
}
View(PRED_20)

write.csv(PRED_20, "PRED_20.csv",row.names = FALSE)
PRED_20<- read.csv("PRED_20.csv")
View(PRED_20)

rmse20 <- sqrt(mean((PRED_20$truevalue - PRED_20$pred)^2))  # 1.345702
```

Model 24- Non phylogenetic model best model

```{r}

ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_24 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM24 <- phylopars.lm(log_c_loads ~  logmerge_mass+logforearm + mean_temp,  tree = ltrns(0), trait_data=.td27$dat, model = "BM")
  
  predicted <- .plm_BM24$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM24$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_24 <- rbind(PRED_24, preds)
}
View(PRED_24)

write.csv(PRED_24, "PRED_24.csv",row.names = FALSE)
PRED_24<- read.csv("PRED_24.csv")
View(PRED_24)

rmse24 <- sqrt(mean((PRED_24$truevalue - PRED_24$pred)^2))  # 1.369235


```
Model 10- Mass and Temp
```{r}

ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_10 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM10 <- phylopars.lm(log_c_loads ~ logmerge_mass+mean_temp,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  
  predicted <- .plm_BM10$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM10$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_10 <- rbind(PRED_10, preds)
}
View(PRED_10)

write.csv(PRED_10, "PRED_10.csv",row.names = FALSE)
PRED_10<- read.csv("PRED_10.csv")
View(PRED_10)

rmse10 <- sqrt(mean((PRED_10$truevalue - PRED_10$pred)^2))  #0.9746389


```
Model 1: Just mass

```{r}

ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_1 <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
.plm_BM1 <- phylopars.lm(log_c_loads ~ logmerge_mass,  tree = .td27$phy, trait_data=.td27$dat, model = "BM")
  
  predicted <- .plm_BM1$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM1$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_1 <- rbind(PRED_1, preds)
}
View(PRED_1)

write.csv(PRED_1, "PRED_1.csv",row.names = FALSE)
PRED_1<- read.csv("PRED_1.csv")
View(PRED_1)

rmse1 <- sqrt(mean((PRED_1$truevalue - PRED_1$pred)^2))  #  1.005352


```

No Traits just tree
```{r}
ltrns <- function(lambda) geiger::lambdaTree(td27$phy, lambda)
PRED_phyloonly <- NULL
for(i in species_with_data){
  .td27 <- td27
  
  refvalue <- .td27$dat$log_c_loads[i]
  .td27$dat$log_c_loads[i]  <- NA
  
  .plm_BM21 <- phylopars.lm(log_c_loads ~ 1,  tree = .td$phy, trait_data=.td$dat, model = "BM")
  summary(.plm_BM21)
  
  predicted <- .plm_BM21$PPE$anc_recon[i,'log_c_loads']
  predicted_var <- .plm_BM21$PPE$anc_var[i,'log_c_loads']
  
  preds <- data.frame(species=.td27$phy$tip.label[i], truevalue=refvalue, pred=predicted, predvar=predicted_var)
  PRED_phyloonly <- rbind(PRED_phyloonly, preds)
}
View(PRED_phyloonly) 

write.csv(PRED_phyloonly, "PRED_phyloonly.csv",row.names = FALSE)
PRED_phyloonly<- read.csv("PRED_phyloonly.csv")
View(PRED_phyloonly)

rmse <- sqrt(mean((PRED_phyloonly$truevalue - PRED_phyloonly$pred)^2)) #1.975135

```

```{r}
mean((PREDgroup$truevalue - PREDgroup$pred)^2) # 1.141533
mean((PREDgroupphy$truevalue - PREDgroupphy$pred)^2) #0.9655226
mean((PRED1$truevalue - PRED1$pred)^2) #0.9064389
mean((PREDfull$truevalue - PREDfull$pred)^2) #  1.53106
mean((PRED_nonphylo$truevalue - PRED_nonphylo$pred)^2) # 2.615667
mean((PRED_phyloonly$truevalue - PRED_phyloonly$pred)^2) # 3.901159

```

using model 12 to make plots
```{r}

library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(tidyr)
library(dplyr)

# Tree
tree <- read.tree("your_tree_file.nwk")
tree <- td$phy

# Trait data with columns: species, male, female
data_long <- PRED1 %>%
  pivot_longer(cols = c(truevalue, pred), names_to = "range", values_to = "LC50")

# Merge to ensure match between tree and data
data_long$species <- factor(data_long$species, levels = tree$tip.label)

# Basic ggtree
p <- ggtree(tree) + geom_tiplab()

# Add LC50 dots per sex using geom_fruit()
p + geom_fruit(
    data = data_long,
    geom = geom_point,
    mapping = aes(y = species, x = LC50, shape = range, fill = range),
    size = 2,
    text.size = 0.8,
    pwidth = 0.5,
    vline=TRUE, 
    show.tiplabels=FALSE
  ) +
  theme_tree2()


```

