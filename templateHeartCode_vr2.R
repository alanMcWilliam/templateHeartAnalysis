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
  tmp$median_sd <- apply(tmp[,2:6],1, sd)
  tmp$sem <- tmp$median_sd/sqrt(5)
  
  substructureName1 <- str_replace(heartDataMean[i], ".csv", "_Mean")
  substructureName2 <- str_replace(heartDataMean[i], ".csv", "_Mean_sd")
  substructureName3 <- str_replace(heartDataMean[i], ".csv", "_sem")
  colnames(tmp) <- c("ID", "1", "2", "3", "4", "5", substructureName1, substructureName2, substructureName3)
  
  keep <- c("ID", substructureName3)#, substructureName2, substructureName3)
  tmp <- tmp[keep]
  all_tmp <- merge(all_tmp, tmp, by = "ID")
}


for(i in 1:length(heartDataMax)) {
  tmpLocation <- paste(dataFolderMax, heartDataMax[i], sep = "\\")
  tmp <- read.csv(tmpLocation, header = FALSE)
  #tmp$max <- apply(tmp[,2:6],1, max)
  tmp$max <- apply(tmp[,2:6],1, mean)
  
  tmp$max_sd <- apply(tmp[,2:6],1, sd)
  tmp$max_sem <- tmp$max_sd/sqrt(5)
  
  substructureName1 <- str_replace(heartDataMean[i], ".csv", "_Max")
  substructureName2 <- str_replace(heartDataMean[i], ".csv", "_Max_sd")
  substructureName3 <- str_replace(heartDataMean[i], ".csv", "_Max_sem")
  
  colnames(tmp) <- c("ID", "1", "2", "3", "4", "5", substructureName1, substructureName2, substructureName3)
  
  keep <- c("ID", substructureName3)#, substructureName2, substructureName3)
  tmp <- tmp[keep]
  all_tmp <- merge(all_tmp, tmp, by = "ID")
}

#remove the single PS 4 patient
all <- all[all$performance.status < 4, ]
#View(all)
#write.csv(all, "C:\\Users\\alan_\\Desktop\\templateHeart\\test.csv" )

#create a copy of the status and follow up at the end of the DF for ease of creating the required datasets
all$stat <- all$status
all$follow_up <- all$followUp

#split into training and validation
bound <- floor((nrow(all)/3)*2)
all <- all[sample(nrow(all)), ]           #sample rows 
all_test <- all[(bound+1):nrow(all), ]    #get test set
all <- all[1:bound, ]                     #get training set

#check balance between tumour sizes between datasets
summary(all$tumour_size)
summary(all_test$tumour_size)
t.test(all$tumour_size, all_test$tumour_size)

#Heart data and status and follow_up
all2 <- all[ ,46:ncol(all)] # without status (no extra lung stats) 75
all2_test <- all_test[ ,46:ncol(all)]
#View(all2)
#View(all2_test)

#remove lines with NA's for RF
all2$comp <- complete.cases(all2)
all2 <- all2[all2$comp == TRUE,]
all2 <- all2[-ncol(all2)]
#View(all2)
 
#remove lines with NA's
all2_test$comp <- complete.cases(all2_test)
all2_test <- all2_test[all2_test$comp == TRUE,]
all2_test <- all2_test[-ncol(all2_test)]
#View(all2_test)

#######################################################################
write.csv(all, "C:\\Users\\alan_\\Desktop\\templateHeart\\all4.csv")
write.csv(all_test, "C:\\Users\\alan_\\Desktop\\templateHeart\\all_test4.csv")


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

# cox proportional hazards model with all heart sub-structures
coxAll <- coxph(Surv(time = followUp, event = status)~rightVentricle_Median + rightVentricle_Max + rightAtrium_Median + rightAtrium_Max + leftVentricle_Median + leftVentricle_Max + leftAtrium_Median +leftAtrium_Max + aorticValve_Median +  aorticValve_Max + ascendingAorta_Median + ascendingAorta_Max + CXA_Median + CXA_Max + LAD_Median + LAD_Max + LMC_Median + LMC_Max + mitralValve_Median + mitralValve_Max + pulmonaryArtery_Median + pulmonaryArtery_Max + pulmonicValve_Median + pulmonicValve_Max + RCA_Median + RCA_Max + tricuspidValve_Median + tricuspidValve_Max + log(tumour_size) + Age + gender + factor(performance.status), data = all)
summary(coxAll)


#####################################################################
# LASSO elastic net 
#https://www.r-bloggers.com/variable-selection-with-elastic-net/
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

#need to remove the last colum, follow up, from the df
all2_EN <- all2[ ,1:(ncol(all2)-1)]
all2_test_EN <- all2_test[ ,1:(ncol(all2_test)-1)]
#View(all2_EN)
#View(all2_test_EN)

set.seed(2017)
#set coefficents for glm 
mdlY <- as.factor(as.matrix(all2_EN["stat"]))
mdlX <- as.matrix(all2_EN[-ncol(all2_EN)])
newY <- as.factor(as.matrix(all2_test_EN["stat"]))
newX <- as.matrix(all2_test_EN[-ncol(all2_test_EN)])

#full LASSO
cv1 <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 1)
md1 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 1)
coef(md1)
#write.table(md1, "C:\\Users\\alan_\\Desktop\\templateHeart\\results\\md1.txt", sep="\t") 
#tmp <- roc(newY, as.numeric(predict(md1, newX, type = "response")))
#plot(tmp)

#rigid
md2 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 0) 
coef(md2)
tmp2 <- roc(newY, as.numeric(predict(md2, newX, type = "response")))
plot(tmp2)
tmp2

#elastic net LASSO
a <- seq(0.1, 0.9, 0.01)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
coef(md3)
rocEN <- roc(newY, as.numeric(predict(md3, newX, type = "response")))

par(pty="s")
plot(rocEN, xlim = c(1,0), ylim = c(0,1), col = 'black', lwd = 2, asp = NA, xlab = "1 - Specificity")
lines(tmp2, col = 'black', lty = 2, lwd = 2)
legend(1, 1, legend=c("Elastic Net", "Rigid"),
       col=c("black", "black"), text.font = 2, box.lty=0, lty=1:2, cex=0.9)
rocEN


############################################################
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

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
all$maxComb <- pmax(all$RCA_Max, all$rightAtrium_Max)
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


opt_cut_mean <- survminer::surv_cutpoint(all, time = "followUp", event = "status", "rightAtrium_Max",  minprop = 0.1, progressbar = TRUE) #smethod="logrank" set within)
summary(opt_cut_mean)
cat_mean <-survminer::surv_categorize(opt_cut_mean)
cat_split_mean <- survfit(Surv(time = followUp, event = status)~rightAtrium_Max, data = cat_mean)
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

