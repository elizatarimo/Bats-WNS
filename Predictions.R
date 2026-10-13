```{r}
getwd()
setwd("~/Desktop/Current Research/Batsfile")
library(rotl)
library(treeplyr)
library(ape)
library(dplyr)
library(tidyr)
library(treeplyr)
library(ape)
library(lme4)
library(phytools)
library(Rphylopars)
library(ggplot2)

```

read in data
```{r}
td63 <- readRDS("td63.RDS")
td63$dat <- data.frame("species"=td63$phy$tip.label, td63$dat)
#View(td63$dat)

#td <- readRDS("td27.RDS")
#View(td$dat) 
```

Best Model Predictions..Model 12

```{r}
plm_BM12_63 <- phylopars.lm(log_c_loads~logmerge_mass+logforearm+mean_temp,  tree = td63$phy, trait_data=td63$dat, model = "BM")

plm_BM12_63$PPE
plm_BM12_63$PPE$anc_var
```

fungal load predictions
```{r}
pred_ranges <- read.csv("predictions_ranges.csv",header=TRUE)

predictions <- as.data.frame(plm_BM12_63$PPE$anc_recon[1:63,'log_c_loads'])
predictions$log_c_loads <-predictions$`plm_BM12_63$PPE$anc_recon[1:63, "log_c_loads"]` #All of the log c loads regardless of predicted or observed
predictions$og_log_c_loads <- td63$dat$log_c_loads
predictions$pred_variance <- plm_BM12_63$PPE$anc_var[1:63,'log_c_loads'] #The variance with zero are the ones where you had data...
predictions$lower_95CI <- predictions$log_c_loads - 1.96*sqrt(pred_variance[,1])
predictions$upper_95CI <- predictions$log_c_loads + 1.96*sqrt(pred_variance[,1])
predictions$Location <- pred_ranges$Location 
predictions <- predictions[,2:7]
View(predictions)

write.csv(predictions,"predictions.csv")
predictions <- read.csv("predictions.csv",header=TRUE)
View(predictions)

```

only species with predicted values
```{r}
predicted_only <- predictions %>%
  filter(is.na(og_log_c_loads))

write.csv(predicted_only, "predicted_only.csv")
predicted_only <- read.csv("predicted_only.csv", header=TRUE)
View(predicted_only) #36

predicted_all <- data.frame(species=predicted_only$X,
                            log_loads=predicted_only$log_c_loads,
                            lower_95CI=predicted_only$lower_95CI,
                            upper_95CI=predicted_only$upper_95CI)

predicted_all <- arrange(predicted_all, log_loads)

write.csv(predicted_all, "predicted_all.csv")
predicted_all <- read.csv("predicted_all.csv", header=TRUE)
View(predicted_all)

#A KS TEST? 
#One-sample KS test against standard normal distribution
ks.test(predicted_only[,2], "pnorm", mean = 0, sd = sqrt(predicted_only[,4]))

#reject the null hypothesis.

```

Look into the Confidence Intervals
```{r}
ci_width <- predicted_all$upper_95CI - predicted_all$lower_95CI

quantile(ci_width, c(0.25, 0.5, 0.75))

hist(ci_width, main="CI width distribution", xlab="CI width")

plot(predicted_all$log_loads, ci_width,
     xlab="Predicted value",
     ylab="CI width")

```


US ONLY
```{r}

predictions_only_US <- predictions %>%
  filter(is.na(og_log_c_loads) & Location == "US")

View(predictions_only_US) #21

write.csv(predictions_only_US,"predictions_only_US.csv")
predictions_only_US <- read.csv("predictions_only_US.csv",header=TRUE)
View(predictions_only_US)

predicted_US <- data.frame(species=predictions_only_US$X,
                            log_loads=predictions_only_US$log_c_loads)
View(predicted_US)

predicted_US <- arrange(predicted_US, log_loads)
write.csv(predicted_US,"predicted_US.csv")
```

best model traits
```{r}
predii <- data.frame(log_c_loads=plm_AllBest$PPE$anc_recon[1:63,'log_c_loads'],
                     mass=plm_AllBest$PPE$anc_recon[1:63,'logmerge_mass'],
                     forearm=plm_AllBest$PPE$anc_recon[1:63,'logforearm'],
                      mean_temp=plm_AllBest$PPE$anc_recon[1:63,'mean_temp'])
View(predii)

write.csv(predii,"predii.csv")
predii <- read.csv("predii.csv",header=TRUE)

#line graphs

```


Diagrams
```{r}

#phenogram
phenogram(td63$phy, plm_BM12_63$PPE$anc_recon[,'log_c_loads'], fcols=is.na(td63$dat$log_c_loads), fsize=0.3)
nH <- max(branching.times(td63$phy))
upper <- sqrt(plm_BM12_63$PPE$anc_var[,'log_c_loads']) + plm_BM12_63$PPE$anc_recon[,'log_c_loads']
lower <- -1*sqrt(plm_BM12_63$PPE$anc_var[,'log_c_loads']) + plm_BM12_63$PPE$anc_recon[,'log_c_loads']

COLS = ifelse(td63$dat$disease._state=="epiizootic","blue","orange")
names(COLS) <- td63$phy$tip.label
#plot(td63$phy,tip.color=COLS[td63$phy$tip.label])
object <- contMap(td63$phy,  plm_BM12_63$PPE$anc_recon[1:63,'log_c_loads'], outline = FALSE)

#differentiating predicted and actual
pdf("./preds_cont.pdf", height=10, width=10)
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
       pos=4,cex=0.6,col=cols[status[i]],font=8)

add.simmap.legend(colors=setNames(cols, c("Known", "Model predicted")),prompt=FALSE,
                  x=60,y=-5.4,vertical=FALSE,fsize=0.8)
#tiplabels(td63$phy$tip.label, col=COLS[td63$phy$tip.label])

dev.off()
system("open preds_cont.pdf")
```

gives you a graph to show the correlation between the predicted values and the actual values

```{r}
par(mar=c(5,4,2,2))
plot(PRED$truevalue, PRED$pred, xlim=c(-10,2), ylim=c(-10,2))
abline(0,1)

corval <- cor.test(PRED$pred, PRED$truevalue)
View(corval)
#corval[["estimate"]]
0.5605851 
#corval[["p.value"]]
0.00235435

```

When you remake the plot, you can set 

```{r}
plot(...., xlab = "True Value", ylab="Predicted Value", xlim=c(-10, 0), ylim=c(-10, 0)) 
#and it will give prettier axis labels than PRED$pred  and cut off the extra 2 units that aren't needed.

```

Use the same phylopars lm code, but replace the tree with ltrns(0)
urooting the tree to 0 to 'make' it a normal linear modal

```{r}
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
```

add phylopars command with just variables from PM19 and species name in dat
```{r}

#Known Vs To be predicted
COLS1 = ifelse(Joinedbats1$log_c_loads=="TRUE","blue","orange")
names(COLS1) <- Joinedbats1$species
plot(pruned.tree,tip.color=COLS1[pruned.tree$tip.label],no.margin=TRUE,edge.width=0.5)

#When you remake the plot, you can set 
plot(...., xlab = "True Value", ylab="Predicted Value", xlim=c(-10, 0), ylim=c(-10, 0)) and it will give prettier axis labels than PRED$pred  and cut off the extra 2 units that aren't needed.
 
 
     