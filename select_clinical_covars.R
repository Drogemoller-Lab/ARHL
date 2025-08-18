library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(olsrr)
library(blorr)
library(BSDA)
library(broom)
library(viridis)
library(stats)
setwd("/home/projects/archive/CLSA/raw_data/clinical_data/")
master <- as.character(unzip("2104035_UManitoba_BDrogemoller_Baseline.zip", list = TRUE)$Name)
PhenData <- read.csv(unz("2104035_UManitoba_BDrogemoller_Baseline.zip",
                         "2104035_UManitoba_BDrogemoller_Baseline/2104035_UManitoba_BDrogemoller_CoP6_Baseline.csv"), header = TRUE,
                     sep = ",")
Dictionary <- read_xlsx("2104035_UManitoba_BDrogemoller_Baseline/2104035_UManitoba_BDrogemoller_CoP6_Baseline-dictionary.xlsx")
metabolic <- read.table("/home/projects/hearing_loss/clsaARHL_SA/phenotypes/metabolic.txt", header = TRUE)
sensory <- read.table("/home/projects/hearing_loss/clsaARHL_SA/phenotypes/sensory.txt", header = TRUE)
#select relevant columns
metabolic.phen <- metabolic %>%
  select(FID, better.rank)
sensory.phen <- sensory %>%
  select(FID, better.rank)
#select clinical covariates: Diabetes, osteoporosis, hypertension, kidney diseases, current home has noise problem, type of smoker,
#served in military
clinical_covar <- select(PhenData, ADM_GWAS3_COM, AGE_NMBR_COM, SEX_ASK_COM, DIA_DIAB_COM, CCC_OSTPO_COM, CCC_HBP_COM,
                         CCC_KIDN_COM, ENV_HMPRB_NOI_MCQ, VET_OCC_COM )
#rename columns
colnames(clinical_covar)[1:2] <- c("FID", "age")
clinical_covar <- subset(clinical_covar, !is.na(FID))
#sex if male Yes, if female No
clinical_covar$sexM <- ifelse(clinical_covar$SEX_ASK_COM == "M", "Yes",
                              ifelse(clinical_covar$SEX_ASK_COM == "F", "No", NA))
#diabetes: 1= yes, 2=no, 8= Don't know or no answer, 9= refused
clinical_covar$diabetes <- ifelse(clinical_covar$DIA_DIAB_COM == 1, "Yes",
                                  ifelse(clinical_covar$DIA_DIAB_COM == 2, "No", NA))
#osteoporosis: 1= yes, 2=no, 8= Don't know or no answer, 9= refused
clinical_covar$osteoporosis <- ifelse(clinical_covar$CCC_OSTPO_COM == 1, "Yes",
                                      ifelse(clinical_covar$CCC_OSTPO_COM == 2, "No", NA))
#hypertension: 1= yes, 2=no, 8= Don't know or no answer, 9= refused
clinical_covar$hypertension <- ifelse(clinical_covar$CCC_HBP_COM == 1, "Yes",
                                      ifelse(clinical_covar$CCC_HBP_COM == 2, "No", NA))
#kidney disease: 1= yes, 2=no, 8= Don't know or no answer, 9= refused
clinical_covar$kidney <- ifelse(clinical_covar$CCC_KIDN_COM == 1, "Yes",
                                ifelse(clinical_covar$CCC_KIDN_COM == 2, "No", NA))
#noise: 1= yes, =no, -8= missing
clinical_covar$noise <- ifelse(clinical_covar$ENV_HMPRB_NOI_MCQ == 1, "Yes",
                               ifelse(clinical_covar$ENV_HMPRB_NOI_MCQ == 0, "No", NA))
#Military: 1= yes in Canada, 2= Yes, outside Canada, 3= yes, inside and outside Canada, 4= No, 8= Don't know or no answer, 9= refused
clinical_covar$military <- ifelse(clinical_covar$VET_OCC_COM == 1 | clinical_covar$VET_OCC_COM == 2 | clinical_covar$VET_OCC_COM == 3, "Yes",
                                  ifelse(clinical_covar$VET_OCC_COM == 4, "No", NA))
#Get rid of unnecessary columns
clinical_covar <- clinical_covar[, - c(3:9)]
#merge the clinical variables with audiogram data
combined.met <- merge(metabolic.phen, clinical_covar, by= "FID")
combined.sen <- merge(sensory.phen, clinical_covar, by = "FID")
#add the sensory column
sensory_to_merge <- select(sensory, FID, better.rank)
colnames(sensory_to_merge)[2] <- "sen.better.rank"
combined.met <- merge(combined.met, sensory_to_merge)
#check the types of data present
sapply(combined.met, class)
# Reclassify incorrectly classified columns
cols <- c(4:10)
combined.met[,cols] <- lapply(combined.met[,cols], factor)
rm(cols)
#separate data into numeric and factor
combined.met_numeric <- combined.met %>% Filter(f = is.numeric)
combined.met_factor <- combined.met %>% Filter(f = is.factor)
#add dependent continuous variable
combined.met_factor <- bind_cols(combined.met_factor, combined.met_numeric['better.rank'])
#Get rid of rows with NA dependent variable column
combined.met_factor <- combined.met_factor[!is.na(combined.met_factor$better.rank), ]
#Identify columns that have only one factor level
colnames(combined.met_factor)[(sapply(combined.met_factor, function(col) length(unique(col))) == 1)]
# Remove them
combined.met_factor <- combined.met_factor[,c(sapply(combined.met_factor, function(col) length(unique(col))) != 1)]
#melt database based on delta max rank for factor
df_met <- reshape2::melt(combined.met_factor, id.vars = "better.rank")
#perform stats and summarize
better.met_factor <- df_met %>% group_by(variable) %>% do(tidy(lm(better.rank ~ value, data=.)))
#melt database based on delta max rank for numeric
df2_met <- reshape2::melt(combined.met_numeric, id.vars = "better.rank")
#perform stats and summarize
better.met_numeric <- df2_met %>% group_by(variable) %>% do(tidy(lm(better.rank ~ value, data=.)))
#forward regression
#metabolic phenotype
#select significant covariants from the previous step
combined.met_sig <- select (combined.met, better.rank, age, sexM, diabetes, osteoporosis, hypertension, sen.better.rank)
#remove missing values
combined.met_sig <- na.omit(combined.met_sig)
#broom::glance(lm(better.rank ~ . , data = combined.met_sig))
model_met <- lm(better.rank ~ ., data = combined.met_sig)
summary(model_met) #age, sex, diabetes
#forward regression
lm<- lm(better.rank ~ 1, data=combined.met_sig)
summary(lm)
#calculate AIC and BIC for this model
AIC(lm)
BIC(lm)
#perform forward regression
step(lm, scope=list(upper=lm(better.rank ~ age + sexM + diabetes + hypertension + osteoporosis +sen.better.rank, data=combined.met_sig)), direction="forward")
#Call:
#  lm(formula = better.rank ~ sen.better.rank + age + diabetes + hypertension, data = combined.met_sig)
############################
################################
#add the sensory column
metabolic_to_merge <- select(metabolic, FID, better.rank)
colnames(metabolic_to_merge)[2] <- "met.better.rank"
combined.sen <- merge(combined.sen, metabolic_to_merge)
#check the types of data present
sapply(combined.sen, class)
# Reclassify incorrectly classified columns
cols <- c(4:10)
combined.sen[,cols] <- lapply(combined.sen[,cols], factor)
rm(cols)
#separate data into numeric and factor
combined.sen_numeric <- combined.sen %>% Filter(f = is.numeric)
combined.sen_factor <- combined.sen %>% Filter(f = is.factor)
#add dependent continuous variable
combined.sen_factor <- bind_cols(combined.sen_factor, combined.sen_numeric['better.rank'])
#Get rid of rows with NA dependent variable column
combined.sen_factor <- combined.sen_factor[!is.na(combined.sen_factor$better.rank), ]
#Identify columns that have only one factor level
colnames(combined.sen_factor)[(sapply(combined.sen_factor, function(col) length(unique(col))) == 1)]
# Remove them
combined.sen_factor <- combined.sen_factor[,c(sapply(combined.sen_factor, function(col) length(unique(col))) != 1)]
#melt database based on delta max rank for factor
df_sen <- reshape2::melt(combined.sen_factor, id.vars = "better.rank")
#perform stats and summarize
better.sen_factor <- df_sen %>% group_by(variable) %>% do(tidy(lm(better.rank ~ value, data=.)))
#melt database based on delta max rank for numeric
df2_sen <- reshape2::melt(combined.sen_numeric, id.vars = "better.rank")
#perform stats and summarize
better.sen_numeric <- df2_sen %>% group_by(variable) %>% do(tidy(lm(better.rank ~ value, data=.)))
#######################################################################
#sensory phenotype
#select significant covariants from the previous step
combined.sen_sig <- select (combined.sen, better.rank, age, sexM, diabetes, hypertension, kidney, military, met.better.rank)
#remove missing values
combined.sen_sig <- na.omit(combined.sen_sig)
#broom::glance(lm(better.rank ~ . , data = combined.sen_sig))
#model_sen <- lm(better.rank ~ . , data = combined.sen_sig)
#summary(model_sen) #age, sex, HTN
#forward regression
lm2<- lm(better.rank ~ 1, data=combined.sen_sig)
summary(lm2)
#calculate AIC and BIC for this model
AIC(lm2)
BIC(lm2)
#perform forward regression
step(lm2, scope=list(upper=lm(better.rank ~ age + sexM + diabetes + hypertension + kidney + military + met.better.rank, data=combined.sen_sig)), direction="forward")
#Call:
#lm(formula = better.rank ~ met.better.rank + age + sexM + hypertension + diabetes, data = combined.sen_sig)

#########################################
#lm(formula = met.better ~ sen.better.rank + age + diabetes + hypertension, data = combined.met_sig)
#lm(formula = sen.better ~ met.sen.better + age + sexM + hypertension + diabetes, data = combined.sen_sig)

#calculate standardized effect sizes
# Load libraries
library(lme4)
library(lmerTest)
library(effectsize) # standard estimates of effect size


# Metabolic analysis (linear model)
MMod <- lm(
  better.rank ~ . - FID, 
  data = combined.met
)
summary(MMod)


# Standard effect sizes 
StMMod = standardize_parameters(MMod)
MSE = StMMod$Std_Coefficient # reorganizing into a vector (opt)
names(MSE) = StMMod$Parameter  # reorganizing into a vector (opt)
round(MSE, 2)

# Sensory analysis (linear model)
SMod <- lm(
  better.rank ~ . - FID, 
  data = combined.sen
)
summary(SMod)


# Standard effect sizes 
StSMod = standardize_parameters(SMod)
SSE = StSMod$Std_Coefficient # reorganizing into a vector (opt)
names(SSE) = StSMod$Parameter  # reorganizing into a vector (opt)
round(SSE, 2)

