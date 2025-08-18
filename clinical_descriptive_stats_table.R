library(dplyr)
#to calculate statistics
#load clinical variables from the covariate file
metabolic <- read.table("/home/projects/hearing_loss/clsa_v2/phenotypes/metabolic.txt", header = TRUE)
sensory <- read.table("/home/projects/hearing_loss/clsa_v2/phenotypes/sensory.txt", header = TRUE)

mean_age<- round(mean(combined.met$age),2)
sd_age <- round(sd(combined.met$age),2)
mean_met <- round(mean(metabolic$met.better),2)
sd_met <- round(sd(metabolic$met.better),2)
mean_sen <- round(mean(sensory$sen.better),2)
sd_sen <- round(sd(sensory$sen.better),2)

sex <- combined.met %>% 
  group_by( sexM ) %>% 
  dplyr::summarise( percent = paste0(round(100 * n() / nrow( combined.met ) ,2), "%"))

diabetes <- combined.met %>% 
  group_by( diabetes ) %>% 
  filter(!is.na(diabetes)) %>%
  dplyr::summarise( percent = paste0(round(100 * n() / nrow( combined.met ) ,2), "%"))

hypertension <- combined.met %>% 
  group_by( hypertension ) %>% 
  filter(!is.na(hypertension)) %>%
  dplyr::summarise( percent = paste0(round(100 * n() / nrow( combined.met ) ,2), "%"))

osteoporosis <- combined.met %>% 
  group_by( osteoporosis ) %>% 
  filter(!is.na(osteoporosis)) %>%
  dplyr::summarise( percent = paste0(round(100 * n() / nrow( combined.met ) ,2), "%"))

noise <- combined.met %>% 
  group_by( noise ) %>% 
  filter(!is.na(noise)) %>%
  dplyr::summarise( percent = paste0(round(100 * n() / nrow( combined.met ) ,2), "%"))

military <-combined.met %>% 
  group_by( military ) %>% 
  filter(!is.na(military)) %>%
  dplyr::summarise( percent = paste0(round(100 * n() / nrow( combined.met ) ,2), "%"))

kidney <-combined.met %>% 
  group_by( kidney ) %>% 
  filter(!is.na(kidney)) %>%
  dplyr::summarise( percent = paste0(round(100 * n() / nrow( combined.met ) ,2), "%"))

summary(combined.met)
