### 29th Sept 2021 updated to git repository
#try random split in time

#library's needed
library(survival)
library(KMsurv)
#library(OIsurv)
library(ggplot2)
library(survminer)
library(stringr)
library(ggcorrplot)
library(plotROC)

###################################################

dataFolderMean = "C:\\Users\\alan_\\Desktop\\templateHeart\\run3burnTrueBlurred"
dataFolderMax = "C:\\Users\\alan_\\Desktop\\templateHeart\\run4MaxTrueBlurred"

#location of heart sub-structure doses
heartDataMean <- list.files(dataFolderMean)
heartDataMax <- list.files(dataFolderMax)

# current spreadsheet of clinical data
all <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\Lateral.csv")
# read in dose data for substructures, calculates the median/max of values and appends to the clinical spreadsheet
for(i in 1:length(heartDataMean)) {
  tmpLocation <- paste(dataFolderMean, heartDataMean[i], sep = "\\")
  tmp <- read.csv(tmpLocation, header = FALSE)
  tmp$median <- apply(tmp[,2:6],1, mean)
  #tmp$median_sd <- apply(tmp[,2:6],1, sd)
  #tmp$sem <- tmp$median_sd/sqrt(5)
  
  substructureName1 <- str_replace(heartDataMean[i], ".csv", "_Mean")
  #substructureName2 <- str_replace(heartDataMean[i], ".csv", "_Mean_sd")
  #substructureName3 <- str_replace(heartDataMean[i], ".csv", "_sem")
  colnames(tmp) <- c("ID", "1", "2", "3", "4", "5", substructureName1)
  
  keep <- c("ID", substructureName1)#, substructureName2, substructureName3)
  tmp <- tmp[keep]
  all <- merge(all, tmp, by = "ID")
}


for(i in 1:length(heartDataMax)) {
  tmpLocation <- paste(dataFolderMax, heartDataMax[i], sep = "\\")
  tmp <- read.csv(tmpLocation, header = FALSE)
  #tmp$max <- apply(tmp[,2:6],1, max)
  tmp$max <- apply(tmp[,2:6],1, mean)
  
  #tmp$max_sd <- apply(tmp[,2:6],1, sd)
  #tmp$max_sem <- tmp$max_sd/sqrt(5)
  
  substructureName1 <- str_replace(heartDataMean[i], ".csv", "_Max")
  #substructureName2 <- str_replace(heartDataMean[i], ".csv", "_Max_sd")
  #substructureName3 <- str_replace(heartDataMean[i], ".csv", "_Max_sem")
  
  colnames(tmp) <- c("ID", "1", "2", "3", "4", "5", substructureName1)
  
  keep <- c("ID", substructureName1)#, substructureName2, substructureName3)
  tmp <- tmp[keep]
  all <- merge(all, tmp, by = "ID")
}

#remove the single PS 4 patient
all <- all[all$performance.status < 4, ]
View(all)
write.csv(all, "C:\\Users\\alan_\\Desktop\\testCombineRegion.csv" )

#create a copy of the status and follow up at the end of the DF for ease of creating the required datasets
all$stat <- all$status
all$follow_up <- all$followUp

#split into training and validation
bound <- floor((nrow(all)/3)*2)
all <- all[sample(nrow(all)), ]           #sample rows 
all_test <- all[(bound+1):nrow(all), ]    #get test set
all <- all[1:bound, ]                     #get training set

#check balance between tumour sizes between datasets
###summary(all$tumour_size)
###summary(all_test$tumour_size)
###t.test(all$tumour_size, all_test$tumour_size)

#Heart data and status and follow_up
all2 <- all[ ,46:ncol(all)] # without status (no extra lung stats) 75
#all2_test <- all_test[ ,46:ncol(all)]
View(all2)
#View(all2_test)

#remove lines with NA's for RF
all2$comp <- complete.cases(all2)
all2 <- all2[all2$comp == TRUE,]
all2 <- all2[-ncol(all2)]
#View(all)
View(all2)

#remove lines with NA's
###all2_test$comp <- complete.cases(all2_test)
###all2_test <- all2_test[all2_test$comp == TRUE,]
###all2_test <- all2_test[-ncol(all2_test)]
#View(all2_test)

#######################################################################
write.csv(all2, "C:\\Users\\alan_\\Desktop\\templateHeart\\all_paper_allData.csv")
#write.csv(all_test, "C:\\Users\\alan_\\Desktop\\templateHeart\\all_test_paper.csv")

#read data used for final paper analysis
all <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\all_paper.csv")
all_test <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\all_test_paper.csv")

#######################################################################
# correlation matrix
#need to remove status and follow up
all2_corr <- all2[ ,1:(ncol(all2)-2)]

corr <- round(cor(all2_corr), 1)
head(corr[,1:30])
p.mat <- cor_pmat(all2_corr)
head(p.mat[,1:30])
ggcorrplot(corr)
ggcorrplot(corr, hc.order = TRUE, lab = TRUE, p.mat = p.mat, insig = "blank")

#correlation matrix for median values and max values alone
corr_mean <- all2_corr[ ,1:14]
corr_max <- all2_corr[ ,15:(ncol(all2_corr))]

corr1 <- round(cor(corr_mean), 1)
ggcorrplot(corr1)
corr2 <- round(cor(corr_max), 1)
ggcorrplot(corr2)


#####################################################################

View(all)
# cox proportional hazards model with all heart sub-structures
coxAll <- coxph(Surv(time = followUp, event = status)~rightVentricle_Mean + rightVentricle_Max + rightAtrium_Mean + rightAtrium_Max + leftVentricle_Mean + leftVentricle_Max + leftAtrium_Mean +leftAtrium_Max + aorticValve_Mean +  aorticValve_Max + ascendingAorta_Mean + ascendingAorta_Max + CXA_Mean + CXA_Max + LAD_Mean + LAD_Max + LMC_Mean + LMC_Max + mitralValve_Mean + mitralValve_Max + pulmonaryArtery_Mean + pulmonaryArtery_Max + pulmonicValve_Mean + pulmonicValve_Max + RCA_Mean + RCA_Max + tricuspidValve_Mean + tricuspidValve_Max + log(tumour_size) + Age + gender + factor(performance.status), data = all)
summary(coxAll)

structures <- c("aorticValve_Mean", "ascendingAorta_Mean", "CXA_Mean", "LAD_Mean", "leftAtrium_Mean", "leftVentricle_Mean", "LMC_Mean", "mitralValve_Mean", "pulmonaryArtery_Mean", "pulmonicValve_Mean", "RCA_Mean", "rightAtrium_Mean", "rightVentricle_Mean", "tricuspidValve_Mean", "aorticValve_Max", "ascendingAorta_Max", "CXA_Max", "LAD_Max", "leftAtrium_Max", "leftVentricle_Max", "LMC_Max", "mitralValve_Max", "pulmonaryArtery_Max", "pulmonicValve_Max", "RCA_Max", "rightAtrium_Max", "rightVentricle_Max", "tricuspidValve_Max")
##testing the cox models
for(i in structures){
  print(i)
  coxUni <- coxph(Surv(time = followUp, event = status)~all[[i]] + log(tumour_size), data = all)
  print(summary(coxUni))
}

coxUni <- coxph(Surv(time = followUp, event = status)~all$CXA_Mean, data = all)
coxUni <- coxph(Surv(time = followUp, event = status)~all$CXA_Mean + all$tumour_size, data = all)
print(summary(coxUni))


library(dplyr)
library(finalfit)

count = 1
for(i in structures){
  print(i)
  explanatory = c("tricuspidValve_Mean", "tumour_size")
  dependent = "Surv(followUp, status)"
  
  all %>%
    finalfit(dependent, explanatory) -> tt
  save(tt, file = paste("C:\\Users\\alan_\\Desktop\\templateHeart\\resultsCox\\", count, ".rda", sep = ""))
  count = count + 1
}
##display table
knitr::kable(tt, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

#####################################################################
# LASSO elastic net 
#https://www.r-bloggers.com/variable-selection-with-elastic-net/
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

performElasticNet <- function(dataAll){
  #set coefficents for glm 
  mdlY <- as.factor(as.matrix(dataAll["stat"]))
  mdlX <- as.matrix(dataAll[-ncol(dataAll)])
  
  coefOut <- matrix(NA, nrow = 3, ncol = 29)
  
  #full LASSO
  cv1 <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 1)
  md1 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 1)
  #print(coef(md1))
  tmp_LASSO <- data.frame(coef.name = dimnames(coef(md1))[[1]], coef.value = matrix(coef(md1)))
  coefOut[1,] <- tmp_LASSO[,2]
  
  #rigid
  md2 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 0) 
  #print(coef(md2))
  tmp_rigid <- data.frame(coef.name = dimnames(coef(md2))[[1]], coef.value = matrix(coef(md2)))
  coefOut[2,] <- tmp_rigid[,2]
  
  #elastic net LASSO
  a <- seq(0.1, 0.9, 0.05) #change final number for fine tuning to be faster was 0.01
  search <- foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  md3 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
  #coef(md3)
  tmp_elasticNet <- data.frame(coef.name = dimnames(coef(md3))[[1]], coef.value = matrix(coef(md3)))
  coefOut[3,] <- tmp_elasticNet[,2]
  
  #coefOut <- c(tmp_LASSO[,2], tmp_rigid[,2], tmp_elasticNet[,2])
  
  #print(coefOut)
  return (coefOut)
}

# B number of bootstraps
bootstrap_r <- function(ds, B) {
  ds <- all2_EN
  # Preallocate storage for statistics
  boot_stat_LASSO <- matrix(NA, nrow = B, ncol = 29) #number for ncol - 30 sub-structures collecting coefficients
  boot_stat_rigid <- matrix(NA, nrow = B, ncol = 29) #number for ncol - 30 sub-structures collecting coefficients
  boot_stat_elasticNet <- matrix(NA, nrow = B, ncol = 29) #number for ncol - 30 sub-structures collecting coefficients
  
  # Number of observations
  n <- nrow(ds)
  
  # Perform bootstrap
  for(i in seq_len(B)) {
    print(i)
    # Sample initial data
    gen_data <- ds[ sample(nrow(ds), replace=TRUE), ]
    # Calculate sample data mean and SD
    coefOut2 <- performElasticNet(gen_data)
    #print(coefOut2)
    
    
    boot_stat_LASSO[i,] <- coefOut2[1,]
    boot_stat_rigid[i,] <- coefOut2[2,]
    boot_stat_elasticNet[i,] <- coefOut2[3,]
  }
  
  boot_stat <- rbind(boot_stat_LASSO, boot_stat_rigid, boot_stat_elasticNet)
  
  # Return bootstrap result
  return(boot_stat)
}

#need to remove the last colum, follow up, from the df
all2_EN <- all2[ ,1:(ncol(all2)-1)]
View(all2_EN)

#set a seed for bootstrapping
set.seed(883)
b = 500
resultsAll <- bootstrap_r(all2_EN, b)
#View(resultsAll)

#split into results for each model, add colum names and convert to data frames
boot_LASSO <- resultsAll[1:b,]
boot_rigid <- resultsAll[(b+1):(2*b),]
boot_elastic <- resultsAll[(2*b+1):nrow(resultsAll),]

colnames(boot_rigid) <- c("Intercept", "aorticValve_Mean", "ascendingAorta_Mean", "CXA_Mean", "LAD_Mean", "leftAtrium_Mean", "leftVentricle_Mean", "LMC_Mean", "mitralValve_Mean", "pulmonaryArtery_Mean", "pulmonicValve_Mean", "RCA_Mean", "rightAtrium_Mean", "rightVentricle_Mean", "tricuspidValve_Mean", "aorticValve_Max", "ascendingAorta_Max", "CXA_Max", "LAD_Max", "leftAtrium_Max", "leftVentricle_Max", "LMC_Max", "mitralValve_Max", "pulmonaryArtery_Max", "pulmonicValve_Max", "RCA_Max", "rightAtrium_Max", "rightVentricle_Max", "tricuspidValve_Max")
colnames(boot_LASSO) <- c("Intercept", "aorticValve_Mean", "ascendingAorta_Mean", "CXA_Mean", "LAD_Mean", "leftAtrium_Mean", "leftVentricle_Mean", "LMC_Mean", "mitralValve_Mean", "pulmonaryArtery_Mean", "pulmonicValve_Mean", "RCA_Mean", "rightAtrium_Mean", "rightVentricle_Mean", "tricuspidValve_Mean", "aorticValve_Max", "ascendingAorta_Max", "CXA_Max", "LAD_Max", "leftAtrium_Max", "leftVentricle_Max", "LMC_Max", "mitralValve_Max", "pulmonaryArtery_Max", "pulmonicValve_Max", "RCA_Max", "rightAtrium_Max", "rightVentricle_Max", "tricuspidValve_Max")
colnames(boot_elastic) <- c("Intercept", "aorticValve_Mean", "ascendingAorta_Mean", "CXA_Mean", "LAD_Mean", "leftAtrium_Mean", "leftVentricle_Mean", "LMC_Mean", "mitralValve_Mean", "pulmonaryArtery_Mean", "pulmonicValve_Mean", "RCA_Mean", "rightAtrium_Mean", "rightVentricle_Mean", "tricuspidValve_Mean", "aorticValve_Max", "ascendingAorta_Max", "CXA_Max", "LAD_Max", "leftAtrium_Max", "leftVentricle_Max", "LMC_Max", "mitralValve_Max", "pulmonaryArtery_Max", "pulmonicValve_Max", "RCA_Max", "rightAtrium_Max", "rightVentricle_Max", "tricuspidValve_Max")

boot_rigid <- as.data.frame(boot_rigid)
boot_LASSO <- as.data.frame(boot_LASSO)
boot_elastic <- as.data.frame(boot_elastic)

View(boot_rigid)
View(boot_LASSO)
View(boot_elastic)


generateStats <- function(df){
  for(i in colnames(df)){
    print(i)
    print(mean(df[[i]]))
    print(summary(df[[i]]))
    print(quantile(df[[i]], probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)))
    print(sum(df[[i]] > 0))
    
    #plot histogram
    ##LASSO y 0, 150, seq -0.05, 0.05
    #ggplot(data=df, aes(df[[i]])) + 
      #geom_histogram(breaks=seq(-0.05,0.05, by = 0.001),
      #               col = "skyblue", fill = "lightblue") +
      #labs(title = i, x = "coefficent" ) +
      #ylim(0,500) +
      #theme(panel.background = element_blank())
    
   # ggsave(paste("C:\\Users\\alan_\\Desktop\\templateHeart\\results3\\", i, ".jpg", sep=""))
  }
}


#read data in from bootstarpping
boot_rigid <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\bootstrapped\\bootRigid500.csv")
boot_LASSO <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\bootstrapped\\bootLASSO500.csv")
boot_elastic <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\bootstrapped\\bootElastic500.csv")


#how to select variables from bootstrapping
generateStats(boot_LASSO)




a <- mean(boot_elastic$RCA_Max)#rightAtrium_Max)
s <- sd(boot_elastic$RCA_Max)#rightAtrium_Max)

error <- qnorm(0.95)*s/sqrt(n)  #normal distribution
left <- a-error
right <- a+error
left
right

error <- qt(0.95, df=n-1)*s/sqrt(n)  #t distribution
left <- a-error
right <- a+error
left
right


    


############################################################
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/


library(ranger)
performRandom <- function(df){
  r_fit <- ranger(Surv(follow_up, stat)~rightVentricle_Mean + rightVentricle_Max + rightAtrium_Mean + rightAtrium_Max + leftVentricle_Mean + leftVentricle_Max + leftAtrium_Mean +leftAtrium_Max + aorticValve_Mean +  aorticValve_Max + ascendingAorta_Mean + ascendingAorta_Max + CXA_Mean + CXA_Max + LAD_Mean + LAD_Max + LMC_Mean + LMC_Max + mitralValve_Mean + mitralValve_Max + pulmonaryArtery_Mean + pulmonaryArtery_Max + pulmonicValve_Mean + pulmonicValve_Max + RCA_Mean + RCA_Max + tricuspidValve_Mean + tricuspidValve_Max, data = df, 
                  mtry = 4,
                  importance = "permutation",
                  splitrule = "extratrees",
                  verbose = TRUE
  )
  
  a <- data.frame(round(r_fit$variable.importance, 4))
  
  return(a)
}

# B number of bootstraps
bootstrap_random <- function(ds, B) {
  # Preallocate storage for statistics
  boot_stat_random <- matrix(NA, nrow = B, ncol = 28) #number for ncol - 30 sub-structures collecting coefficients
  
  # Number of observations
  n <- nrow(ds)
  
  # Perform bootstrap
  for(i in seq_len(B)) {
    print(i)
    # Sample initial data
    gen_data <- ds[ sample(nrow(ds), replace=TRUE), ]
    # Calculate sample data mean and SD
    randomOut <- performRandom(gen_data)
    #print(coefOut2)
    
    boot_stat_random[i,] <- randomOut[,1]
    
  }
  
  # Return bootstrap result
  return(boot_stat_random)
}

set.seed(890)
bb = 500
resultsAll_random <- bootstrap_random(all2, bb)

colnames(resultsAll_random) <- c("rightVentricle_Mean", "rightVentricle_Max", "rightAtrium_Mean", "rightAtrium_Max", "leftVentricle_Mean", "leftVentricle_Max", "leftAtrium_Mean", "leftAtrium_Max", "aorticValve_Mean",  "aorticValve_Max", "ascendingAorta_Mean", "ascendingAorta_Max", "CXA_Mean", "CXA_Max", "LAD_Mean", "LAD_Max", "LMC_Mean", "LMC_Max", "mitralValve_Mean", "mitralValve_Max", "pulmonaryArtery_Mean", "pulmonaryArtery_Max", "pulmonicValve_Mean", "pulmonicValve_Max", "RCA_Mean", "RCA_Max", "tricuspidValve_Mean", "tricuspidValve_Max")
View(resultsAll_random)

##write.csv(resultsAll_random, "C:\\alan\\templateHeart\\bootRandom.csv")

bootRandom <- read.csv("C:\\Users\\alan_\\Desktop\\templateHeart\\resultsRandom\\bootRandom.CSV")

summary(bootRandom)

aa <- colnames(bootRandom)[apply(bootRandom,1,which.max)]
summary(factor(aa))

for(i in colnames(bootRandom)){
  ggplot(data=bootRandom, aes(bootRandom[[i]])) + 
    geom_histogram(aes(bootRandom$rightAtrium_Max), breaks=seq(0,0.05, by = 0.001),
                   col = "blue", fill = "blue", alpha=0.5) +
    geom_histogram(aes(bootRandom[[i]]), breaks=seq(0,0.05, by = 0.001),
                   col = "skyblue", fill = "skyblue", alpha=0.5) +
    labs(title = i, x = "importance" ) +
    ylim(0,150) +
    theme(panel.background = element_blank())
  
  ggsave(paste("C:\\Users\\alan_\\Desktop\\templateHeart\\resultsRandom\\", i, ".jpg", sep=""))
  
  
}



######------------------------------------------------------------
library(ranger)
r_fit <- ranger(Surv(follow_up, stat)~rightVentricle_Mean + rightVentricle_Max + rightAtrium_Mean + rightAtrium_Max + leftVentricle_Mean + leftVentricle_Max + leftAtrium_Mean +leftAtrium_Max + aorticValve_Mean +  aorticValve_Max + ascendingAorta_Mean + ascendingAorta_Max + CXA_Mean + CXA_Max + LAD_Mean + LAD_Max + LMC_Mean + LMC_Max + mitralValve_Mean + mitralValve_Max + pulmonaryArtery_Mean + pulmonaryArtery_Max + pulmonicValve_Mean + pulmonicValve_Max + RCA_Mean + RCA_Max + tricuspidValve_Mean + tricuspidValve_Max, data = all2, 
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE
                )

sort(r_fit$variable.importance)

vi <- data.frame(sort(round(r_fit$variable.importance, 4), decreasing = TRUE))
names(vi) <- "importance"
head(vi)


a <- data.frame(round(r_fit$variable.importance, 4))

cat("Prediction Error = 1 - Harrell's c-index = ", r_fit$prediction.error)

importance_pvalues(r_fit, method = c("altmann"), num.permutations = 100, formula = NULL, data = all2)

#validate model for events at one year
#create a variable of patients dead or alive at 12 months for prediction analysis to work
all2_test$reached_event <- ifelse((all2_test$stat == 1 & all2_test$follow_up <= 12),1,0)

#'HeartJK_Max', 'HeartJK_Mean',
survival_predicition <- predict(r_fit, all2_test[, c('rightVentricle_Mean', 'rightVentricle_Max', 'rightAtrium_Mean', 'rightAtrium_Max', 'leftVentricle_Mean', 'leftVentricle_Max', 'leftAtrium_Mean', 'leftAtrium_Max', 'aorticValve_Mean',  'aorticValve_Max', 'ascendingAorta_Mean', 'ascendingAorta_Max', 'CXA_Mean', 'CXA_Max', 'LAD_Mean', 'LAD_Max', 'LMC_Mean', 'LMC_Max', 'mitralValve_Mean', 'mitralValve_Max', 'pulmonaryArtery_Mean', 'pulmonaryArtery_Max', 'pulmonicValve_Mean', 'pulmonicValve_Max', 'RCA_Mean', 'RCA_Max', 'tricuspidValve_Mean', 'tricuspidValve_Max')])
rocRF <- roc(response = all2_test$reached_event, predictor = 1 - survival_predicition$survival[,which(survival_predicition$unique.death.times == 12)])
plot(rocRF, xlim = c(1,0), ylim = c(0,1), col = 'black', lwd = 2, asp = NA, xlab = "1 - Specificity")
legend(1, 1, legend=c("Random Forest"),
       col=c("black", "black"), text.font = 2, box.lty=0, lty=1:2, cex=0.9)
#plot(rocRF)
rocRF

print(r_fit)


################################################################################
#LASSO elastic net shows right atrium max dose  and RCA_max is the only remaining paramter
#RF survival analysis showed right artium max and RCA_max dose had the greatest influence for survival

summary(all$rightAtrium_Max)
summary(all$RCA_Max)

coxUni <- coxph(Surv(time = followUp, event = status)~rightAtrium_Max, data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~RCA_Max, data = all)
summary(coxUni)

#max in the RCA and RA
all$maxComb <- pmax(all$RCA_Max, all$rightAtrium_Max, all$ascendingAorta_Max)
summary(all$maxComb)
all_test$maxComb <- pmax(all_test$RCA_Max, all_test$rightAtrium_Max)
summary(all_test$maxComb)

coxUni <- coxph(Surv(time = followUp, event = status)~maxComb, data = all)
summary(coxUni)

#other clinical univariate analysis
coxUni <- coxph(Surv(time = followUp, event = status)~log(tumour_size), data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~Age, data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~gender, data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~factor(performance.status), data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~factor(N.stage_edit), data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~factor(T.stage_edit), data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~chemo, data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~ mean.lung, data = all)
summary(coxUni)
coxUni <- coxph(Surv(time = followUp, event = status)~ lateral, data = all)
summary(coxUni)

#simple multivariate analysis
coxMulti <- coxph(Surv(time = followUp, event = status)~maxComb + log(tumour_size) + Age + gender + factor(performance.status), data = all)
summary(coxMulti)
#coxMulti <- coxph(Surv(time = followUp, event = status)~rightAtrium_Max + log(tumour_size) + Age + gender + factor(performance.status) + factor(N.stage_edit) + chemo, data = all)
#summary(coxMulti)

#all clinical factors
coxMulti <- coxph(Surv(time = followUp, event = status)~maxComb  + log(tumour_size) + Age + gender + factor(performance.status) + factor(N.stage_edit) + factor(T.stage_edit) + chemo + mean.lung + lateral, data = all)
summary(coxMulti)
#coxMulti <- coxph(Surv(time = followUp, event = status)~rightAtrium_Max  + log(tumour_size) + Age + gender + factor(performance.status) + factor(N.stage_edit) + factor(T.stage_edit) + chemo + mean.lung, data = all)
#summary(coxMulti)
#coxMulti <- coxph(Surv(time = followUp, event = status)~rightAtrium_Max + RCA_Max  + log(tumour_size) + Age + gender + factor(performance.status) + factor(N.stage_edit) + factor(T.stage_edit) + chemo + mean.lung, data = all)
#summary(coxMulti)

#forest plot
ggforest(coxMulti)

#quartiles for right atrium dose
all <- within(all, quartileRA <- as.integer(cut(rightAtrium_Max, quantile(rightAtrium_Max, probs=0:4/4, na.rm = TRUE)), include.lowest=TRUE))
tapply(all$rightAtrium_Max, all$quartile, summary)
SurvRA <- survfit(Surv(time = followUp, event = status)~quartileRA, data = all)
ggsurvplot(SurvRA, risk.table = TRUE, conf.int = TRUE, pval = FALSE, ncensor.plot = FALSE)
#quartiles for right coronary dose
all <- within(all, quartileRC <- as.integer(cut(RCA_Max, quantile(RCA_Max, probs=0:4/4, na.rm = TRUE)), include.lowest=TRUE))
tapply(all$RCA_Max, all$quartile, summary)
SurvRC <- survfit(Surv(time = followUp, event = status)~quartileRC, data = all)
ggsurvplot(SurvRC, risk.table = TRUE, conf.int = TRUE, pval = FALSE, ncensor.plot = FALSE)

#optimal cut for right atrium max 
opt_cut_mean <- survminer::surv_cutpoint()


opt_cut_mean <- survminer::surv_cutpoint(all, time = "followUp", event = "status", "maxComb",  minprop = 0.1, progressbar = TRUE) #smethod="logrank" set within)
summary(opt_cut_mean)
cat_mean <-survminer::surv_categorize(opt_cut_mean)
cat_split_mean <- survfit(Surv(time = followUp, event = status)~maxComb, data = cat_mean)
ggsurvplot(cat_split_mean, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)
#optimal cut for right coronary max 
opt_cut_mean2 <- survminer::surv_cutpoint(all, time = "followUp", event = "status", "RCA_Max",  minprop = 0.4, progressbar = TRUE) #smethod="logrank" set within)
summary(opt_cut_mean2)
cat_mean2 <-survminer::surv_categorize(opt_cut_mean2)
cat_split_mean2 <- survfit(Surv(time = followUp, event = status)~RCA_Max, data = cat_mean2)
ggsurvplot(cat_split_mean2, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)

#optimal cut training data RCA 20.5 and right Atrium 24.5 <- mean 22.5
#split validation dataset based on this value
all_test$RAsplit <- all_test$rightAtrium_Max > 22.5
all_test$RCsplit <- all_test$RCA_Max > 22.5
View(all_test)

#plot kaplan meier curves with test dataframe split at the optimal point
km_RCA <- survfit(Surv(time = followUp, event = status)~RCsplit, data = all_test)
ggsurvplot(km_RCA, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)
km_RA <- survfit(Surv(time = followUp, event = status)~RAsplit, data = all_test)
ggsurvplot(km_RA, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)
km_RCA
km_RA

##############################################################################



opt_cut_mean3 <- survminer::surv_cutpoint(all, time = "followUp", event = "status", "maxComb",  minprop = 0.1, progressbar = TRUE) #smethod="logrank" set within)
summary(opt_cut_mean3)
cat_mean3 <-survminer::surv_categorize(opt_cut_mean3)
cat_split_mean3 <- survfit(Surv(time = followUp, event = status)~maxComb, data = cat_mean3)
ggsurvplot(cat_split_mean3, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)

all_test$CombSplit <- all_test$maxComb > 19.5
km_comb <- survfit(Surv(time = followUp, event = status)~CombSplit, data = all_test)
ggsurvplot(km_comb, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)
km_comb

summary(factor(all$lateral))
tapply(all$maxComb, all$lateral, summary)

##################################################################################

wilcox.test(all$tumour_size, all_test$tumour_size)
wilcox.test(all$Age, all_test$Age)
ks.test(all$Age, all_test$Age)

wilcox.test(all$N.stage_edit, all_test$N.stage_edit)
wilcox.test(all$performance.status, all_test$performance.status)
wilcox.test(all$mean.lung, all_test$mean.lung)

wilcox.test(all$T.stage_edit, all_test$T.stage_edit)

all$gender_flag <- all$gender == 'Male' 
all_test$gender_flag <- all_test$gender == 'Male'
all$gender_flag.num <- as.numeric(all$gender_flag)
all_test$gender_flag.num <- as.numeric(all_test$gender_flag)
wilcox.test(all$gender_flag.num, all_test$gender_flag.num)
typeof(all$gender_flag)

all$chemo_flag <- as.numeric(all$chemo)
all_test$chemo_flag <- as.numeric(all_test$chemo)
wilcox.test(all$chemo_flag, all_test$chemo_flag)


##########################################################################################

ggplot(data= all, aes(all$tumour_size)) +
  geom_histogram(breaks=seq(0,500, by=10), col = "skyblue", fill = "skyblue") +
  labs(x = "tumour volume, cm3") +
  theme(panel.background = element_blank())

ggplot(data= all, aes(log(all$tumour_size))) +
  geom_histogram(breaks=seq(0,7, by=0.1), col = "skyblue", fill = "skyblue") +
  labs(x = "log tumour volume, cm3") +
  theme(panel.background = element_blank())

ggplot(all, aes(x= aorticValve_Mean, y = rightAtrium_Mean)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Aortic Valve mean dose (Gy)", y = "Right Atrium mean dose (Gy)") +
  theme(panel.background = element_blank())

cor.test(all$aorticValve_Mean, all$rightAtrium_Mean, method = c("spearman"))

         