
```{r}
getwd()
setwd("~/Desktop/Current Research/Batsfile")

install.packages("rotl")
install.packages("remotes")
install.packages("devtools")
devtools::install_github("uyedaj/treeplyr")
install.packages("ape")
install.packages("lme4")
install.packages("phytools")
install.packages("Rphylopars")
install.packages("xtable")
install.packages("geiger")

```


```{r}
library(ape)
library(remotes)
library(devtools)
library(dplyr)
library(tidyr)
library(readxl)
library(dplyr)
library(rotl)
library(treeplyr)
library(ape)
library(lme4)
library(phytools)
library(Rphylopars)
library(xtable)
library(geiger)
library(ggplot2)
```
load in data
```{r}
td <- readRDS("td27.RDS")
#row.names(td$dat) <- td$phy$tip.label
View(td$dat)


td$dat$loc <- ifelse(td$dat$Location=="US", 1, 0)
td$dat$Group <- ifelse(td$dat$group=="ingroup", 1, 0)

```
How much signal does log_c_loads have?
```{r}
# Create a new data frame
td_loads <- data.frame(td$dat$log_c_loads)
row.names(td_loads) <- td$phy$tip.label
td_loads$log_c_loads <- td_loads$td.dat.log_c_loads

# Calculate phylogenetic signal using phytools
?phylosig
k_result <- phylosig(td$phy, td$dat$log_c_loads, method = "K") # 1.2266 
 
#lambda
lambda_result <- phylosig(td$phy, td$dat$log_c_loads, method = "lambda") # 0.780231 

```

what about other traits?....mass doesn't have the highest phylogenetic signal!
retain all the traits to model their influence on fungal loads
```{r}
# Create a new data frame
td_rownames <- data.frame(td$dat)
row.names(td_rownames) <- td$phy$tip.label
View(td_rownames)


# Calculate phylogenetic signal using phytools
new_dat <- data.frame(species=td$phy$tip.label,
                      log_c_loads=td$dat$log_c_loads,
                      logmerge_mass=td$dat$logmerge_mass,
                      logforearm=td$dat$logforearm,
                      logcluster=td$dat$logcluster, 
                      logtouching=td$dat$logtouching,
                      logpopsize=td$dat$logpopsize,
                      mean_temp=td$dat$mean_temp)


fit_all  <- phylopars(trait_data = new_dat, tree = td$phy, model = "BM")
fit_full <- as.data.frame(fit_all$anc_recon[1:27,]) #to note
fit_full

k_forearm <- phylosig(td$phy, fit_full$logforearm, method = "K")  # 0.954771
k_cluster <- phylosig(td$phy, fit_full$logcluster, method = "K")  # 0.584204 
k_touch <- phylosig(td$phy, fit_full$logtouching, method = "K")  # 0.61642 
k_mass <- phylosig(td$phy, fit_full$logmerge_mass, method = "K") #0.762047  # alone 0.799964 #here 0.798467 
k_pop <- phylosig(td$phy, fit_full$logpopsize, method = "K")  # 0.738781 
k_temp <- phylosig(td$phy, fit_full$mean_temp, method = "K")  # 0.996086 


#lambda
lambda_result <- phylosig(td$phy, td$dat$log_c_loads, method = "lambda") # 0.780231 

```


combine all the models. How is log_c_loads influenced by the predictor variables and what are the evolutionary relationships between them?

```{r}
ltrns <- function(lambda) geiger::lambdaTree(td$phy, lambda)
all_phylopars_models <- list(
  
  #ones
plm_BM1 <- phylopars.lm(log_c_loads ~ logmerge_mass,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM2 <- phylopars.lm(log_c_loads ~ logforearm,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM3 <- phylopars.lm(log_c_loads ~ logcluster,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM4 <- phylopars.lm(log_c_loads ~ logtouching,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM5 <- phylopars.lm(log_c_loads ~ logpopsize,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM6 <- phylopars.lm(log_c_loads ~ mean_temp,  tree = td$phy, trait_data=td$dat, model = "BM"),

#DOUBLE
plm_BM7 <- phylopars.lm(log_c_loads ~ logmerge_mass+logforearm,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM8 <- phylopars.lm(log_c_loads ~ mean_temp+logcluster,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM9 <- phylopars.lm(log_c_loads ~ logmerge_mass+logcluster,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM10 <- phylopars.lm(log_c_loads ~ logmerge_mass+mean_temp,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM11 <- phylopars.lm(log_c_loads ~ logforearm+logcluster,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM12 <- phylopars.lm(log_c_loads ~ logpopsize+mean_temp,  tree = td$phy, trait_data=td$dat, model = "BM"),

#with 3 variables
plm_BM13 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM14_lambda <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp,  tree = td$phy, trait_data=td$dat, model = "lambda"), #13

plm_BM15_EB <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp,  tree = td$phy, trait_data=td$dat, model = "EB"), #14
#plm_BM12sq <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + I(mean_temp^2),  tree = td$phy, trait_data=td$dat, model = "BM"),


plm_BM16 <- phylopars.lm(log_c_loads ~  logmerge_mass  + mean_temp + logcluster,  tree = td$phy, trait_data=td$dat, model = "BM"),

#without mass and forearm_length
plm_BM17 <- phylopars.lm(log_c_loads ~  mean_temp + logcluster + logtouching,  tree = td$phy, trait_data=td$dat, model = "BM"),

#without touching
plm_BM18 <- phylopars.lm(log_c_loads ~ logmerge_mass + logforearm + mean_temp + logcluster,  tree = td$phy, trait_data=td$dat, model = "BM"),

#without forearm
plm_BM19 <- phylopars.lm(log_c_loads ~  logmerge_mass + mean_temp + logcluster + logtouching ,  tree = td$phy, trait_data=td$dat, model = "BM"),

#with all variables
plm_BM20 <- phylopars.lm(log_c_loads ~ logmerge_mass + logforearm + mean_temp + logpopsize + logcluster + logtouching ,  tree = td$phy, trait_data=td$dat, model = "BM"),

plm_BM21 <- phylopars.lm(log_c_loads ~ logmerge_mass + logforearm + mean_temp + logcluster + logtouching,  tree = td$phy, trait_data=td$dat, model = "BM"),

#phylogeny only model
plm_BM22 <- phylopars.lm(log_c_loads ~ 1,  tree = td$phy, trait_data=td$dat, model = "BM"),

#Non Phylogeny model
plm_BM23 <- phylopars.lm(log_c_loads ~   logmerge_mass + logforearm + mean_temp +logcluster + logtouching,  tree = ltrns(0), trait_data=td$dat, model = "BM"),

#Non phylogenetic model best model
plm_BM24 <- phylopars.lm(log_c_loads ~  logmerge_mass+logforearm + mean_temp,  tree = ltrns(0), trait_data=td$dat, model = "BM"),

#variations of best model all phylo
#location best model
plm_BM25_loc <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp+ loc,  tree = td$phy, trait_data=td$dat, model = "BM"),

#endemic vs epizootic
#plm_BM26_ende <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp+ endemism,  tree = td$phy, trait_data=td$dat, model = "BM"),

##outgroup vs ingroup
plm_BM27_group0 <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + Group,  tree = ltrns(0), trait_data=td$dat, model = "BM"),

plm_BM28_group <- phylopars.lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp + Group,  tree = td$phy, trait_data=td$dat, model = "BM")

)


```
get model parameters

```{r}
extract_model_summaries <- function(all_phylopars_models) {
  # Initialize empty lists to store information for each model
  Rsquared <- list()
  PValue <- list()
  Sigma <- list()
  Loglik <- list()
  
  # Loop through each model in the list
  for (i in seq_along(all_phylopars_models)) {
    model <- all_phylopars_models[[i]]
    
    # Extract relevant values from the current model
    Rsquared[[i]] <- model$R2adj
    PValue[[i]] <- model$pval
    Sigma[[i]] <- model$sigma
    Loglik[[i]] <- model$logLik
  }
  
  # Convert each list of values into a data frame
  Rsquared_df <- data.frame(Model = seq_along(Rsquared), Rsquared = unlist(Rsquared))
  PValue_df <- data.frame(Model = seq_along(PValue), PValue = unlist(PValue))
  Sigma_df <- data.frame(Model = seq_along(Sigma), Sigma = unlist(Sigma))
  Loglik_df <- data.frame(Model = seq_along(Loglik), Loglik = unlist(Loglik))
  
  # Combine all data frames into one summary table
  allmodelssummary <- data.frame(
    Model = seq_along(all_phylopars_models),
    Rsquared = unlist(Rsquared),
    PValue = unlist(PValue),
    Sigma = unlist(Sigma),
    Loglik = unlist(Loglik)
  )
  
  # Return the combined summary table
  return(allmodelssummary)
}

```

list of models
```{r}
allmodelssummary <- extract_model_summaries(all_phylopars_models)
print(allmodelssummary)

write.csv(allmodelssummary, "allmodelssummary.csv", row.names = FALSE)
allmodelssummary <- read.csv("allmodelssummary.csv",header = TRUE)
View(allmodelssummary)

```

Save to a tex file
```{r}
summarytable <- xtable(allmodelssummary)
sink("allmodelssummary.tex")
print(summarytable)
sink()
```

arrange models with lowest R squared values to the highest
```{r}
arrange_models <- allmodelssummary%>%arrange(Rsquared)
View(arrange_models)
 
 #top15
arranged_models <- arrange_models[14:28,] 
View(arranged_models)
```


naming the models
```{r}

arranged_models[arranged_models$Model == "1", "Model"] <- "~ mass "
arranged_models[arranged_models$Model == "10", "Model"] <- "~ mass + temperature"
arranged_models[arranged_models$Model == "24", "Model"] <- "~ mass + forearm + temperature (no phylogeny)"
arranged_models[arranged_models$Model == "20", "Model"] <- "~ mass + forearm + temperature + population size + cluster + touch"
arranged_models[arranged_models$Model == "25", "Model"] <- "~ mass + forearm + temperature + continent"
arranged_models[arranged_models$Model == "23", "Model"] <- "~ mass + forearm + temperature + touch + cluster (no phylogeny)"
arranged_models[arranged_models$Model == "26", "Model"] <- "~ mass + forearm + temperature + endemism" 
arranged_models[arranged_models$Model == "7", "Model"] <- "~ mass + forearm"
arranged_models[arranged_models$Model == "21", "Model"] <- "~ mass + forearm + temperature + touch + cluster"
arranged_models[arranged_models$Model == "18", "Model"] <- "~ mass + forearm + temperature + cluster"
arranged_models[arranged_models$Model == "14", "Model"] <- "~ mass + forearm + temperature (lambda)"
arranged_models[arranged_models$Model == "13", "Model"] <- "~ mass + forearm + temperature (BM)"
arranged_models[arranged_models$Model == "15", "Model"] <- "~ mass + forearm + temperature (EB)"
arranged_models[arranged_models$Model == "28", "Model"] <- "~ mass + forearm + temperature + group"
arranged_models[arranged_models$Model == "27", "Model"] <- "~ mass + forearm + temperature + group (no phylogeny)"

View(arranged_models)


```

saving for export
```{r}
write.csv(arranged_models, "arranged_models.csv", row.names = FALSE)
arranged_models <- read.csv("arranged_models.csv",header = TRUE)
View(arranged_models)

pdf("arranged_models.pdf")
grid.table(arranged_models)
dev.off()

```

Get model coefficients
```{r}
extract_all_coefficients <- function(all_phylopars_models) {
  # Extract coefficient vectors
  coef_list <- lapply(all_phylopars_models, function(m) m$coefficients)
  
  # Convert each to a one-row data frame
  coef_df_list <- lapply(coef_list, function(x) as.data.frame(t(x)))
  
  # Bind rows, filling missing columns with NA
  all_coefficients <- dplyr::bind_rows(coef_df_list, .id = "Model")
  
  return(all_coefficients)
}


 all_coefficients <- extract_all_coefficients(all_phylopars_models)
print( all_coefficients)
View(all_coefficients)

write.csv(all_coefficients, "all_coefficients.csv", row.names = FALSE)
all_coefficients <- read.csv("all_coefficients.csv",header = TRUE)
View(all_coefficients)
```

```{r}
setwd("~/Desktop/Current Research/Batsfile/figures")
pdf("./load_temp.pdf", height=10, width=10)
plot(plm_group27$original_data$mean_temp, plm_group27$original_data$log_c_loads, pch=21, bg=plm_group27$original_data$Group, xlab="Mean Hibernation Temperature", ylab="Log c loads", xlim=c(-1, 14), ylim=c(-10, -2))

int_g1 <- mean(plm_group27$original_data$logmerge_mass[plm_group27$original_data$Group==0], na.rm=TRUE)*plm_group27$coefficients[2]  + mean(plm_group27$original_data$logforearm[plm_group27$original_data$Group==0], na.rm=TRUE)*plm_group27$coefficients[3] + plm_group27$coefficients[1]

int_g2 <- mean(plm_group27$original_data$logmerge_mass[plm_group27$original_data$Group==1], na.rm=TRUE)*plm_group27$coefficients[2]  + mean(plm_group27$original_data$logforearm[plm_group27$original_data$Group==1], na.rm=TRUE)*plm_group27$coefficients[3] + plm_group27$coefficients[1] + plm_group27$coefficients[5]

abline(int_g1, plm_group27$coefficients[4], lty=2)
abline(int_g2, plm_group27$coefficients[4], lty=1)

system("open load_temp.pdf")

```


figure to show the relationship between temperature and fungal load for each group. 
```{r}
setwd("~/Desktop/Current Research/Batsfile/figures")
pdf("./load_temp1.pdf", height=10, width=10)
plot(plm_group27$original_data$mean_temp, plm_group27$original_data$log_c_loads, pch=21, bg=plm_group27$original_data$Group, xlab="Mean Hibernation Temperature", ylab="Log c loads", xlim=c(-1, 14), ylim=c(-10, -2))

int_g1 <- mean(plm_group27$original_data$logmerge_mass[plm_group27$original_data$Group==0], na.rm=TRUE)*plm_group27$coefficients[2]  + mean(plm_group27$original_data$logforearm[plm_group27$original_data$Group==0], na.rm=TRUE)*plm_group27$coefficients[3] + plm_group27$coefficients[1]

int_g2 <- mean(plm_group27$original_data$logmerge_mass[plm_group27$original_data$Group==1], na.rm=TRUE)*plm_group27$coefficients[2]  + mean(plm_group27$original_data$logforearm[plm_group27$original_data$Group==1], na.rm=TRUE)*plm_group27$coefficients[3] + plm_group27$coefficients[1] + plm_group27$coefficients[5]

abline(int_g1, plm_group27$coefficients[4], lty=2)
abline(int_g2, plm_group27$coefficients[4], lty=1)

# shade range between x and x
#xrange <- c(5, 9)
#yrange <- par("usr")[3:4]  # y-axis limits of plot

#rect(xrange[1], yrange[1], xrange[2], yrange[2], col=rgb(0,0,1,0.2), border=NA)

dev.off()

system("open load_temp1.pdf")

#Caption: 
# Log c loads plotted against mean hibernation temperature from the original data, with color indicating our two phylogenetic groups (Verspertilionidae vs. every one else). Lines indicate the model prediction from our best-fitting model with other variables (log mass and log forearm length) set at their group means. Over the range of hibernation temperatures in our dataset (0.9 C to 12 C), the model predicts within groups an increase of log c loads of 2.32 units. 

```

Figure with non-phylogenetic model (simply lm with no grouping) to show that without phylogeny this is how the model would look like 
```{r}
ltrns <- function(lambda) geiger::lambdaTree(td$phy, lambda)
model <- lm(log_c_loads ~  logmerge_mass + logforearm + mean_temp, data=td$dat)
plm_BM23 <- phylopars.lm(log_c_loads ~  logmerge_mass+logforearm + mean_temp,  tree = ltrns(0), trait_data=td$dat, model = "BM")

#abline(int_g1, plm_group27$coefficients[4], lty=2)

temp_seq <- seq(
  min(td$dat$mean_temp, na.rm = TRUE),
  max(td$dat$mean_temp, na.rm = TRUE),
  length.out = 100
)

newdata <- data.frame(
  logmerge_mass = mean(td$dat$logmerge_mass, na.rm = TRUE),
  logforearm = mean(td$dat$logforearm, na.rm = TRUE),
  mean_temp = temp_seq
)

coefs <- plm_BM23$coefficients

predicted <- coefs["(Intercept)"] +
             coefs["logmerge_mass"] * newdata$logmerge_mass +
             coefs["logforearm"] * newdata$logforearm +
             coefs["mean_temp"] * newdata$mean_temp


plot(td$dat$mean_temp, td$dat$log_c_loads, pch = 16, col = "gray",
     xlab = "Mean Temperature", ylab = "Log C Loads")
lines(temp_seq, predicted, col = "blue", lwd = 2)

```


Fit a BM model on mass to estimate paramteters.
```{r}
traits_plm_BM12 <- as.data.frame(trait_3 <- as.matrix(plm_BM12$PPE$anc_recon[1:27,-4]))

#get mass
mass <- traits_plm_BM12$logmerge_mass


?fitContinuous
fit_bm <- fitContinuous(
  phy = td$phy, 
  dat = setNames(mass, td$dat[, "species"]), 
  model = "BM"
)

summary(fit_bm)
fit_bm$opt

#sigsq = 0.005737832
#z0 = 2.220846
#lnL = -15.61504
#aicc =  35.73007


#FOREARM LENGTH
forearm <- traits_plm_BM12$logforearm


?fitContinuous
forearm_bm <- fitContinuous(
  phy = td$phy, 
  dat = setNames(forearm, td$dat[, "species"]), 
  model = "BM"
)

forearm_bm$opt

#aicc -31.79524


# TEMP
temp <- traits_plm_BM12$mean_temp


?fitContinuous
temp_bm <- fitContinuous(
  phy = td$phy, 
  dat = setNames(temp, td$dat[, "species"]), 
  model = "BM"
)

temp_bm$opt

# sigsq = 0.16987
#z0 = 7.098403
#lnL = -61.35239
#aicc = 127.2048

```

Simulate the MASS trait evolution 100 times using the estimated sigsq and z0
```{r}

num_simulations <- 100
simulations <- replicate(num_simulations, 
                         fastBM(tree = td$phy, 
                                sig2 = fit_bm$opt$sigsq, 
                                a = fit_bm$opt$z0))

# simulations will be a matrix with 100 columns (simulations) and the number of tips (rows)
# You can convert it to a data frame if needed
simulations_df <- as.data.frame(simulations)

# View the first few simulated datasets
head(simulations_df)
View(simulations_df)
simulations_df$log_c_loads <- td$dat$log_c_loads

simulations_df <- data.frame("species"=td$phy$tip.label, simulations_df)
write.csv(simulations_df, "simulations_df.csv", row.names = FALSE)
simulations_df <- read.csv("simulations_df.csv",header = TRUE)
View(simulations_df)

```


use the simulated values in a phylopars.lm BM model
```{r}

# Initialize result storage
Rsquared <- numeric(100)
PValue   <- numeric(100)
Sigma    <- numeric(100)
Loglik   <- numeric(100)

# Loop through each predictor (V1 to V100)
for (i in 1:100) {
  predictor_col <- paste0("V", i)

  # Prepare trait_data with species, log_c_loads, and current predictor
  sim_dat <- data.frame(
    species       = simulations_df$species,
    log_c_loads   = simulations_df$log_c_loads,
    sim_predictor = simulations_df[[predictor_col]]
  )

  # Fit the model
  model <- phylopars.lm(
    log_c_loads ~ sim_predictor,
    tree = td$phy,
    trait_data = sim_dat,
    model = "BM"
  )

  # Store model stats
  Rsquared[i] <- model$R2adj
  PValue[i]   <- model$pval
  Sigma[i]    <- model$sigma
  Loglik[i]   <- model$logLik
}

# Create final summary table
allmassmodels <- data.frame(
  Predictor = paste0("V", 1:100),
  Rsquared  = Rsquared,
  PValue    = PValue,
  Sigma     = Sigma,
  Loglik    = Loglik
)

# Preview
head(allmassmodels)


write.csv(allmassmodels, "allmassmodels.csv", row.names = FALSE)
allmassmodels <- read.csv("allmassmodels.csv",header = TRUE)

View(allmassmodels)

```

Simulate the FOREARM trait evolution 100 times using the estimated sigsq and z0
```{r}

num_simulations <- 100
simulations_forearm <- replicate(num_simulations, 
                         fastBM(tree = td$phy, 
                                sig2 = forearm_bm$opt$sigsq, 
                                a = forearm_bm$opt$z0))

# simulations will be a matrix with 100 columns (simulations) and the number of tips (rows)
# You can convert it to a data frame if needed
simulations_fa <- as.data.frame(simulations_forearm)

# View the first few simulated datasets
head(simulations_fa)
View(simulations_fa)

#add original log_c_loads column
simulations_fa$log_c_loads <- td$dat$log_c_loads

simulations_fa <- data.frame("species"=td$phy$tip.label, simulations_fa)
write.csv(simulations_fa, "simulations_fa.csv", row.names = FALSE)
simulations_fa <- read.csv("simulations_fa.csv",header = TRUE)
View(simulations_fa)

#NOW RUN PHYLOPARS.LM Models and exctract their R squared values

# Initialize result storage
Rsquared <- numeric(100)
PValue   <- numeric(100)
Sigma    <- numeric(100)
Loglik   <- numeric(100)

# Loop through each predictor (V1 to V100)
for (i in 1:100) {
  predictor_col <- paste0("V", i)

  # Prepare trait_data with species, log_c_loads, and current predictor
  sim_dat <- data.frame(
    species       = simulations_fa$species,
    log_c_loads   = simulations_fa$log_c_loads,
    sim_predictor = simulations_fa[[predictor_col]]
  )

  # Fit the model
  model <- phylopars.lm(
    log_c_loads ~ sim_predictor,
    tree = td$phy,
    trait_data = sim_dat,
    model = "BM"
  )

  # Store model stats
  Rsquared[i] <- model$R2adj
  PValue[i]   <- model$pval
  Sigma[i]    <- model$sigma
  Loglik[i]   <- model$logLik
}

# Create final summary table
allfamodels <- data.frame(
  Predictor = paste0("V", 1:100),
  Rsquared  = Rsquared,
  PValue    = PValue,
  Sigma     = Sigma,
  Loglik    = Loglik
)

# Preview
head(allfamodels)


write.csv(allfamodels, "allfamodels.csv", row.names = FALSE)
allfamodels <- read.csv("allfamodels.csv",header = TRUE)

View(allfamodels)


```

Simulate the TEMP trait evolution 100 times using the estimated sigsq and z0
```{r}

num_simulations <- 100
simulations_temp <- replicate(num_simulations, 
                         fastBM(tree = td$phy, 
                                sig2 = temp_bm$opt$sigsq, 
                                a = temp_bm$opt$z0))

# simulations will be a matrix with 100 columns (simulations) and the number of tips (rows)
# You can convert it to a data frame if needed
simulations_temp <- as.data.frame(simulations_temp)

# View the first few simulated datasets
head(simulations_temp)
View(simulations_temp)

#add original log_c_loads column
simulations_temp$log_c_loads <- td$dat$log_c_loads

simulations_temp <- data.frame("species"=td$phy$tip.label, simulations_temp)
write.csv(simulations_temp, "simulations_temp.csv", row.names = FALSE)
simulations_temp <- read.csv("simulations_temp.csv",header = TRUE)
View(simulations_temp)


#NOW RUN PHYLOPARS.LM Models and exctract their R squared values
# Initialize result storage
Rsquared <- numeric(100)
PValue   <- numeric(100)
Sigma    <- numeric(100)
Loglik   <- numeric(100)

# Loop through each predictor (V1 to V100)
for (i in 1:100) {
  predictor_col <- paste0("V", i)

  # Prepare trait_data with species, log_c_loads, and current predictor
  sim_dat <- data.frame(
    species       = simulations_temp$species,
    log_c_loads   = simulations_temp$log_c_loads,
    sim_predictor = simulations_temp[[predictor_col]]
  )

  # Fit the model
  model <- phylopars.lm(
    log_c_loads ~ sim_predictor,
    tree = td$phy,
    trait_data = sim_dat,
    model = "BM"
  )

  # Store model stats
  Rsquared[i] <- model$R2adj
  PValue[i]   <- model$pval
  Sigma[i]    <- model$sigma
  Loglik[i]   <- model$logLik
}

# Create final summary table
alltempmodels <- data.frame(
  Predictor = paste0("V", 1:100),
  Rsquared  = Rsquared,
  PValue    = PValue,
  Sigma     = Sigma,
  Loglik    = Loglik
)

# Preview
head(alltempmodels)


write.csv(allfamodels, "alltempmodels.csv", row.names = FALSE)
alltempmodels <- read.csv("alltempmodels.csv",header = TRUE)

View(alltempmodels)


```
PLOTS
```{r}
allrsquared <- data.frame(allmassmodels$Rsquared)

pdf("./mass_rsquared.pdf", height=10, width=10)
ggplot(allrsquared, aes(x = allmassmodels.Rsquared)) +
  geom_density() +
  geom_vline(xintercept = 0.06870994, color = "red", linetype = "dashed") +
  labs(title = "Density Plot of Rsquared from 100 mass simulations", x = "R squared", y = "Density")
#the number 0.06870994 is the R squared from the mass phylopars model

dev.off()

system("open mass_rsquared.pdf" )


#FOR FOREARM

all_fa_rsquared <- data.frame(allfamodels$Rsquared)

pdf("./fa_rsquared.pdf", height=10, width=10)
ggplot(all_fa_rsquared, aes(x = allfamodels.Rsquared)) +
  geom_density() +
  geom_vline(xintercept =  -0.029531 , color = "red", linetype = "dashed") +
  labs(title = "Density Plot of Rsquared from 100 forearm-length simulations", x = "R squared", y = "Density")
#the number -0.029531 is the R squared from the mass phylopars model

dev.off()

system("open fa_rsquared.pdf" )


#FOR TEMPERATURE

all_temp_rsquared <- data.frame(alltempmodels$Rsquared)

pdf("./temp_rsquared.pdf", height=10, width=10)
ggplot(all_temp_rsquared, aes(x = alltempmodels.Rsquared)) +
  geom_density() +
  geom_vline(xintercept = 0.02052457 , color = "red", linetype = "dashed") +
  labs(title = "Density Plot of Rsquared from 100 temp simulations", x = "R squared", y = "Density")
#the number 0.02052457 is the R squared from the mass phylopars model

dev.off()

system("open temp_rsquared.pdf" )


```

PHYLO PCA of the best model.. for mass, forearm and temp
PCA gives you which traits have the most variation. 
```{r}
#Perform phylogenetic PCA. Here I am using the anc recon above because it has imputed missing values that phylo pca needs. No NA
trait_3 <- as.matrix(plm_BM12$PPE$anc_recon[1:27,-4])
trait_data <- as.matrix(plm_BM12$PPE$anc_recon[1:27,])
gtrait_data <- as.matrix(plm_group26$PPE$anc_recon[1:27,-5]) #best model with group

pca_3_result <- phyl.pca(td$phy, trait_3, mode = "cov")
pca_result <- phyl.pca(td$phy, trait_data, mode = "cov")
gpca_result <- phyl.pca(td$phy, gtrait_data, mode = "cov")

# Examine and plot results
eigenvalues3 <- pca_3_result$V
loadings3 <- pca_3_result$L
plot(pca_3_result, main = "Phylogenetic PCA")
pca_3_result$Evec
#PC1 Temp, PC2 Mass AND PC3 Forearmlength

#with log_c-loads
eigenvalues <- pca_result$V
loadings <- pca_result$L
plot(pca_result, main = "Phylogenetic PCA")

eigenvalues_g <- gpca_result$V
loadings_g <- gpca_result$L
plot(gpca_result, main = "Phylogenetic PCA")
```

```{r}
par(mfrow=c(2,3))
par(mfrow=c(1,1))
sm <- make.simmap(td$phy, td[["disease._state"]])

pdf("pc_plot.pdf", width = 5, height = 5)  
phytools::phylomorphospace(sm, 
                           pc_BM12$scores[,1:2],
                           xlim=c(-20, 20), 
                           ylim=c(-10,10), 
                           xlab='PC 1',ylab='PC 2',
                           node.by.map=TRUE, 
                           colors=setNames(c("red", "black"), c("endemic", "epiizootic")),
                           ftype='off')
   dev.off()                        

phytools::phylomorphospace(sm, 
                           pc_BM12$scores[,c(2,3)],
                           xlim=c(-15, 10), ylim=c(-10,10), 
                           node.by.map=TRUE, 
                           colors=setNames(c("blue", "red"), c("endemic", "epiizootic")))
                           
phytools::phylomorphospace(sm,
                           pc_BM12$scores[,2:3], 
                           xlim=c(-15, 10), ylim=c(-10,10), 
                           node.by.map=TRUE, colors=setNames(c("blue", "red"), c("endemic", "epiizootic")))

?phytools::phylomorphospace

phytools::fancyTree(td$phy, type="scattergram", X=pc_BM12$scores[,1:2], ftype="off")
phytools::fancyTree(td$phy, type="scattergram", X=pc_BM12$scores[,1:2],colors=setNames(c("red", "black"), c("endemic", "epiizootic")), ftype="off")
```

```{r}
# Save the plot as a TikZ file
tikz("plot.tex", width = 5, height = 5)  # Specify the filename and dimensions
plot(1:10, rnorm(10), main = "Sample Plot", xlab = "X-axis", ylab = "Y-axis")
# Close the TikZ device
dev.off()
```


```{r}
biplot(pc_BM12, choices = c(2,3))

autoplot(data.frame(pc_BM12$scores), colour=td$dat$disease._state, loadings=TRUE)
summary(pc_BM12)
```

Covariances of the variables
```{r}
str(plm_BM12$PPE)
asr <- plm_BM12$PPE$anc_recon
```

how many spp do you have?

```{r}
nrow(plm_BM12$original_data) #27 

```
tiny garphs of the variables against each other to show correlation

```{r}
plot(asr[1:27,'logmerge_mass'], asr[1:27,'log_c_loads'], pch=21, col=1+!is.na(plm_BM12$original_data$log_c_loads))
pairs(asr)

correlations <- pdf(pairs(asr))
sink("correlations.pdf")
print(correlations)
sink()


asr19 <- plm_BM19$PPE$anc_recon
pairs(asr19)

?pch
```

A tree with the bar plots to show the actual fungals-re , and the predicted-black
```{r}
barcols <- 1+!is.na(plm_BM12$original_data$log_c_loads)
plotTree.wBars(td$phy,asr[1:27,"log_c_loads"]-min(asr[1:27,"log_c_loads"]), col=barcols)
```
#args.plotTree=list(ftype="off"),
args.barplot = list(col = barcols, border = barcols, xlab = "", xlim = c(0, 70)) 
args.axis = list(at = seq(0, 1, by = 0.2))

plot(tree, ftype="off")
legend(x="topright",legend=names(colors),pch=22,pt.cex=2,pt.bg=colors,
       box.col="transparent")