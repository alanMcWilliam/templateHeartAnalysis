library(OIsurv)
library(ggplot2)
library(survminer)
library(stringr)

dataFolderMean = "C:\\Users\\alan_\\Desktop\\templateHeart\\run3burnTrueBlurred"
dataFolderMax = "C:\\Users\\alan_\\Desktop\\templateHeart\\run4MaxTrueBlurred"

#location of heart sub-structure doses
heartDataMean <- list.files(dataFolderMean)
heartDataMax <- list.files(dataFolderMax)

# current spreadsheet of clinical data
all <- read.csv("C:\\Users\\alan_\\Desktop\\To do\\laterality\\Lateral.csv")

# read in dose data for substructures, calculates the median/max of values and appends to the clinical spreadsheet
for(i in 1:length(heartDataMean)) {
  tmpLocation <- paste(dataFolderMean, heartDataMean[i], sep = "\\")
  tmp <- read.csv(tmpLocation, header = FALSE)
  tmp$median <- apply(tmp[,2:6],1, median)
  
  substructureName <- str_replace(heartDataMean[i], ".csv", "_Median")
  colnames(tmp) <- c("ID", "1", "2", "3", "4", "5", substructureName)
  
  keep <- c("ID", substructureName)
  tmp <- tmp[keep]
  all <- merge(all, tmp, by = "ID")
}

for(i in 1:length(heartDataMax)) {
  tmpLocation <- paste(dataFolderMax, heartDataMax[i], sep = "\\")
  tmp <- read.csv(tmpLocation, header = FALSE)
  tmp$max <- apply(tmp[,2:6],1, max)
  
  substructureName <- str_replace(heartDataMean[i], ".csv", "_Max")
  colnames(tmp) <- c("ID", "1", "2", "3", "4", "5", substructureName)
  
  keep <- c("ID", substructureName)
  tmp <- tmp[keep]
  all <- merge(all, tmp, by = "ID")
}

all <- all[all$performance.status < 4, ]

all$stat <- all$status
all$follow_up <- all$followUp

#Load and merge lung spreadsheet
tmp <- read.csv("C:\\Users\\alan_\\Desktop\\lungV5\\V5results.csv")
all <- merge(all, tmp, by = "ID")
#Remove wrong left lung segmentations
all <- all[all$leftVolume < 9000,]

#split into training and validation
bound <- floor((nrow(all)/4)*3)
all <- all[sample(nrow(all)), ]           #sample rows 
all_test <- all[(bound+1):nrow(all), ]    #get test set
all <- all[1:bound, ]                     #get training set

View(all)

#check balance between sides
summary(all$tumour_size)
summary(all_test$tumour_size)
t.test(all$tumour_size, all_test$tumour_size)

#Create dataframe of only heart values
PCAall <- all[ ,46:76]#ncol(all)] # with status (no extra lung stats)
PCAall <- all[ ,46:75]#ncol(all)] # without status (no extra lung stats)

#Heart data removed
PCAall <- all[ ,46:73]#ncol(all)] # without status (no extra lung stats)
PCAall <- all[ ,46:74]#ncol(all)] # without status (no extra lung stats)

#Heart data and status and follow_up
all2 <- all[ ,46:77]#ncol(all)] # without status (no extra lung stats)
View(all2)
 
# remove lines with NA's
PCAall$comp <- complete.cases(PCAall)
PCAall <- PCAall[PCAall$comp == TRUE,]
PCAall <- PCAall[-ncol(PCAall)]
View(PCAall)

# remove lines with NA's for RF
all2$comp <- complete.cases(all2)
all2 <- all2[all2$comp == TRUE,]
all2 <- all2[-ncol(all2)]
View(all2)


#######################################################################

#install.packages("ggcorrplot")
library(ggcorrplot)

# correlation matrix
corr <- round(cor(PCAall), 1)
head(corr[,1:14])
p.mat <- cor_pmat(PCAall)
head(p.mat[,1:14])
ggcorrplot(corr)
ggcorrplot(corr, hc.order = TRUE, lab = TRUE, p.mat = p.mat, insig = "blank")




##########################################################################

# try and do a PCA....
prin_comp <- prcomp(PCAall, scale. = T)
names(prin_comp)

#outputs the mean of variables
prin_comp$center

#outputs the standard deviation of variables
prin_comp$scale

prin_comp$rotation
prin_comp$rotation[1:15,1:4]

dim(prin_comp$x)

biplot(prin_comp, scale = 1)
biplot(prin_comp,  xlim=c(-5, 10), ylim=c(-5, 5))

#compute standard deviation of each principal component
std_dev <- prin_comp$sdev

#compute variance
pr_var <- std_dev^2

#check variance of first 10 components
pr_var[1:15]
#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

#scree plot
plot(prop_varex, xlab = "Principal Component",
       ylab = "Proportion of Variance Explained",
       type = "b")
#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       type = "b")


######################################################################
# PCA option 2
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")
# need dataframe with only values - PCAall
# remove heart states

# dataframe with mean and max only
PCA_max <- PCAall[, 16:ncol(PCAall)]
View(PCA_max)

PCA_median <- PCAall[ ,1:14]
PCA_median$stat <- PCAall$stat
View(PCA_median)

#change df here...
PCAall <- PCA_max


View(PCAall)
#PCAall[31]
# PCA
res.pca <- PCA(PCAall[,-ncol(PCAall)], graph = FALSE)
res.pca

# get eigenvalues, varience and percentage contribution
eig.val <- get_eigenvalue(res.pca)
eig.val

# Plots a scree plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

#return for variables
var <- get_pca_var(res.pca)
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)


# Coordinates of variables
head(var$coord, 28)
#[plot variables]
fviz_pca_var(res.pca, col.var = "black")

#Quality of variables
head(var$cos2, 31)

install.packages("corrplot")
library(corrplot)
corrplot(var$cos2, is.corr=FALSE)
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:1)

#contribution of variables
head(var$contrib, 31)
corrplot(var$contrib, is.corr=FALSE)    


# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

#The total contribution to PC1 and PC2 is obtained with the following R code:
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)

fviz_pca_var(res.pca, col.var = "contrib",
            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)



# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 4, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             
             legend.title = "Cluster")



#use a grouping variable
# remove the stats variable
heart.pca <- PCA(PCAall[,-31], graph = FALSE)
heart.pca <- PCA(PCAall[,-ncol(PCAall)], graph = FALSE)


fviz_pca_ind(heart.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = factor(PCAall$stat), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             axes = c(2,3)
)

# Variables on dimensions 2 and 3
fviz_pca_var(res.pca, axes = c(2, 3))
# Individuals on dimensions 2 and 3
fviz_pca_ind(res.pca, axes = c(2, 3))

fviz_pca_biplot(heart.pca, 
                col.ind = factor(PCAall$stat), palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species") 


fviz_pca_biplot(heart.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = factor(PCAall$stat),
                col.ind = "black",
               
                
                legend.title = list(fill = "Species", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors


###################################################################

#https://www.statmethods.net/advstats/factor.html


# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors,
# with varimax rotation
fit <- factanal(PCAall, 3, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2
load <- fit$loadings[,1:2]
plot(load,type="n") # set up plot
text(load,labels=names(PCAall),cex=.7) # add variable names 

# Principal Axis Factor Analysis
library(psych)
fit <- fa(PCAall, nfactors = 3)#, rotation = "varimax")
  #factor.pa(PCAall, nfactors=3)#, rotation="varimax")
fit # print results 

# Determine Number of Factors to Extract
install.packages("nFactors")
library(nFactors)
ev <- eigen(cor(PCAall)) # get eigenvalues
ap <- parallel(subject=nrow(PCAall),var=ncol(PCAall),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS) 



# PCA Variable Factor Map
library(FactoMineR)
result <- PCA(PCAall) # graphs generated automatically 




###################################################################################

#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

install.packages("factoextra")
library(factoextra)

res.pca <- prcomp(PCAall, scale = TRUE)
#plot scree plot
fviz_eig(res.pca)

#Graph of individuals. Individuals with a similar profile are grouped together.
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
              )

# graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# biplot of individuals and variables
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


library(factoextra)
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 


PCAall <- all[ ,46:76]
groups <- as.factor(PCAall$stat)
View(PCAall)

#Plot by grouping variable
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

# graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

library(magrittr) # for pipe %>%
library(dplyr)   # everything else
# 1. Individual coordinates
res.ind <- get_pca_ind(res.pca)
# 2. Coordinate of groups
coord.groups <- res.ind$coord %>%
  as_data_frame() %>%
  select(Dim.1, Dim.2) %>%
  mutate(stat = groups) %>%
  group_by(stat) %>%
  summarise(
    Dim.1 = mean(Dim.1),
    Dim.2 = mean(Dim.2)
  )
coord.groups


######################################################################

#http://statweb.stanford.edu/~tibs/superpc/tutorial.html

#install.packages("superpc")
library(superpc)
# PCAall

#split into training and validation
bound <- floor((nrow(all)/4)*3)
PCAall <- PCAall[sample(nrow(PCAall)), ]           #sample rows 
PCAall_test <- PCAall[(bound+1):nrow(PCAall), ]    #get test set
PCAall <- PCAall[1:bound, ]                     #get training set

colnames(PCAall)[colnames(PCAall) == 'stat'] <- 'censoring.status'
View(PCAall2)

PCAall2 <- as.matrix(PCAall)
PCAall.obj<- superpc.train(PCAall, type="survival")





######################################################################
# LASSO elastic net 
#https://www.r-bloggers.com/variable-selection-with-elastic-net/

#install.packages("glmnet")
pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

df1 <- PCAall #need status as final colum 
# remove lines with NA's
df1$comp <- complete.cases(df1)
df1 <- df1[df1$comp == TRUE,]
df1 <- df1[-ncol(df1)]
View(df1)
summary(df1)


# dataframe with mean and max only
df1_max <- df1[, 16:ncol(df1)]
View(df1_max)

df1_median <- df1[ ,1:15]
df1_median$stat <- df1$stat
View(df1_median)

#change df here...
df1 <- df1_median


set.seed(2017)

#split into training and testing dataframes
n <- nrow(df1)
sample <- sample(seq(n), size = n * 0.5, replace = FALSE)
train <- df1[sample, -1]
test <- df1[-sample, -1]

#set coefficents for glm 
mdlY <- as.factor(as.matrix(train["stat"]))
mdlX <- as.matrix(train[-ncol(train)])
newY <- as.factor(as.matrix(test["stat"]))
newX <- as.matrix(test[-ncol(test)])

#full LASSO
cv1 <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 1)
md1 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 1)
coef(md1)
roc(newY, as.numeric(predict(md1, newX, type = "response")))

#rigid
md2 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 0) 
coef(md2)
tmp2 <- roc(newY, as.numeric(predict(md2, newX, type = "response")))
plot(tmp2)

#elastic net LASSO
a <- seq(0.1, 0.9, 0.01)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
coef(md3)
tmp <- roc(newY, as.numeric(predict(md3, newX, type = "response")))
plot (tmp)

################################################################################

#LASSO elastic net shows right atrium max dose is the only remaining paramter
#Right Atrium Max dose
summary(all$rightAtrium_Max)

cox <- coxph(Surv(time = followUp, event = status)~rightAtrium_Max + log(tumour_size) + Age + gender + factor(performance.status), data = all)
summary(cox)

#optimal cut for right atrium max 
opt_cut_mean <- survminer::surv_cutpoint(all, time = "followUp", event = "status", "rightAtrium_Max",  minprop = 0.1, progressbar = TRUE) #smethod="logrank" set within)
summary(opt_cut_mean)
cat_mean <-survminer::surv_categorize(opt_cut_mean)
cat_split_mean <- survfit(Surv(time = followUp, event = status)~rightAtrium_Max, data = cat_mean)
ggsurvplot(cat_split_mean, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)


############################################################
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

install.packages("ranger")
library(ranger)
r_fit <- ranger(Surv(follow_up, stat)~HeartJK_Max + HeartJK_Median + rightVentricle_Median + rightVentricle_Max + rightAtrium_Median + rightAtrium_Max + leftVentricle_Median + leftVentricle_Max + leftAtrium_Median +leftAtrium_Max + aorticValve_Median +  aorticValve_Max + ascendingAorta_Median + ascendingAorta_Max + CXA_Median + CXA_Max + LAD_Median + LAD_Max + LMC_Median + LMC_Max + mitralValve_Median + mitralValve_Max + pulmonaryArtery_Median + pulmonaryArtery_Max + pulmonicValve_Median + pulmonicValve_Max + RCA_Median + RCA_Max + tricuspidValve_Median + tricuspidValve_Max, data = all2, 
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE,
                )


vi <- data.frame(sort(round(r_fit$variable.importance, 4), decreasing = TRUE))
names(vi) <- "importance"
head(vi)

cat("Prediction Error = 1 - Harrell's c-index = ", r_fit$prediction.error)

# do we need training and validation cohorts here - same as the elastic net work
# survival analysis with training - quartiles for right atruim dose?

summary(all$rightAtrium_Max)


#############################################

cox <- coxph(Surv(time = followUp, event = status)~rightAtrium_Max + log(tumour_size) + Age + gender + factor(performance.status) + factor(N.stage_edit) + chemo, data = all)
summary(cox)

cox <- coxph(Surv(time = followUp, event = status)~HeartJK_Median + log(tumour_size) + Age + gender + factor(performance.status) + factor(N.stage_edit) + chemo, data = all)
summary(cox)

