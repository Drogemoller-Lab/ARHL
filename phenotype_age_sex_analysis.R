#determine in the difference between males and females in met and sen phenotype is statistically significant
library(rstatix)
met <- read.table("/home/projects/hearing_loss/clsaARHL_SA/phenotypes/metabolic.txt", header = T)
sen <- read.table("/home/projects/hearing_loss/clsaARHL_SA/phenotypes/sensory.txt", header = T)

cor.test(sen$Age_Nm, sen$sen.better)
cor.test(met$Age_Nm, met$met.better)
#try two.sided t-test
t.test(met$met.better ~ met$Sex, alternative = "two.sided") #females have significantly greater metabolic 
t.test(sen$sen.better ~ sen$Sex, alternative = "two.sided") #males have significanlty greater sensory
#calculate the effect size
met |> cohens_d(met.better ~ Sex)
sen |> cohens_d(sen.better ~ Sex)

#perform t test 
t.test(met$met.better ~ met$Sex, alternative = "greater") #females have significantly greater metabolic 
t.test(sen$sen.better ~ sen$Sex, alternative = "less") #males have significantly greater sensory
