#-----------------------------------------------------------------------------#
# P U L M O N A R Y   R E H A B   C L I N I C A L               s c r i p t   #
#                                                                             #
# Author: Alex                                                                #
# Date created: 29th November 2019                                            #
#-----------------------------------------------------------------------------#




sink() # Just putting this here so that if I run it over again it doesn't create more and more sinks...

filename <- "D:/Alex/PR/PR clinical 2019 interim/Logs/PR_clinical_2019_cleaning_log_"
filedate <- Sys.Date()

sink(file = paste(filename, filedate, ".txt", sep = ""),
     append = FALSE,
     split = TRUE)

cat("\n START \n") # This means that every time I run it it restarts the document instead of getting an
# unuseable document at the end

sink()

sink(file = paste(filename, filedate, ".txt", sep = ""),
     append = TRUE,
     split = TRUE)



library(dplyr)
# library(readstata13)
# library(xlsx)
source("H:/My R functions/MySummary.R")
source("H:/My R functions/lintestOR.R")
source("H:/My R functions/tidyoutput.R")
# library(janitor)
# library(officer)
# library(flextable)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(survsup)
# library(epitools)
library(psych)
library(lme4)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(car)
library(extrafont)
loadfonts()
fonts()

tablex <- function(x, y, z) { x %>% select(!!y, !!z) %>% table(useNA = "ifany") }

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}


meanSumRound <- function(x, variable, roundno) {
  variable <- as.character(variable)
  varcol <- filter(psychic, vars == variable) %>% 
    dplyr::select(vars, N, mean, sd)
  varcol[ ,3:4] <- format(round(varcol[ ,3:4], roundno), nsmall = roundno)
colnames(varcol) <- paste(variable, colnames(varcol), sep = "_")
return(varcol[ , -1])

}

mediSumRound <- function(x, variable, roundno) {
    variable <- as.character(variable)
    varcol <- filter(psychic, vars == variable) %>% 
      dplyr::select(vars, N, median, lo.quart, hi.quart)
    # function updated so that it just gives numbers back rounded according to roundno,
    # without making any exceptions for midway points etc
    varcol[ ,3:5] <- sprintf(paste0("%.", roundno, "f"), 
                             round(varcol[ ,3:5], roundno), nsmall = roundno) # otherwise use 'roundno'
    
    colnames(varcol) <- paste(variable, colnames(varcol), sep = "_")
    return(varcol[ , -1])
  }



FreqSum <- function(x, varname) {
  
  varname <- as.character(varname)
  gen <- x %>% dplyr::select(!!varname) %>% drop_na()
  var_N <- data.frame(nrow(gen))
  colnames(var_N) <- paste0(varname, "_N")
  
# Doesn't matter if nrow = 0... still important to know
# So, I've hashed out the 'if 0' bit, and also the close bracket at the end
  
  
 # if(nrow(gen) == 0) {return(var_N)}
  
 # else {
    
    gen0 <- as.data.frame(table(gen[[1]]))
    gen1 <- as.data.frame(round(prop.table(table(gen[[1]]))*100, 1), nsmall = 1) %>% 
      dplyr::rename(perc = Freq)
    gen2 <- inner_join(gen0, gen1, by = "Var1")
    gen2$perc <- sprintf("%.1f", gen2$perc)
    # gen.E2$England <- paste(gen.E2$Freq, " (", gen.E2$perc, ")", sep = "")
    # gen.E2 <- select(gen.E2, Var1, England)
    for (i in 1:nrow(gen2)) {
      gen3 <- gen2
      gen3$Var1 <- as.character(gen3$Var1)
      gen3 <- gen3[i, ]
      colnames(gen3) <- c("Var1", paste0(varname, "_", gsub(" ", "_", gen3[1,1]), "_n"),
                          paste0(varname, "_", gsub(" ", "_", gen3[1,1]), "_perc")) 
      var_N <- cbind(var_N, gen3[ ,2:3])
    }
    return(var_N)
    
 # }
}



medTable <- function(x, varname) {   
  # x is the dataset, varname is the variable name, val is the value of interest (e.g. males) 
  varname <- as.character(varname)
  
  eng <- x %>% filter(country == "England") %>% dplyr::select(varname)
  EN <- length(eng[!is.na(eng)])
  engIQR <- quantile(eng[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  eng <- paste(engIQR[2], " (", engIQR[1], " to ", engIQR[3], ")", sep = "")
  
  
  wal <- x %>% filter(country == "Wales") %>% dplyr::select(varname)
  WN <- length(wal[!is.na(wal)])
  walIQR <- quantile(wal[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  wal <- paste(walIQR[2], " (", walIQR[1], " to ", walIQR[3], ")", sep = "")
  
  
  scot <- x %>% filter(country == "Scotland") %>% dplyr::select(varname)
  SN <- length(scot[!is.na(scot)])
  scotIQR <- quantile(scot[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  scot <- paste(scotIQR[2], " (", scotIQR[1], " to ", scotIQR[3], ")", sep = "")
  
  
  all <- x %>% dplyr::select(varname)
  AN <- length(all[!is.na(all)])
  allIQR <- quantile(all[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  all <- paste(allIQR[2], " (", allIQR[1], " to ", allIQR[3], ")", sep = "")
  
  ret <- matrix(c(eng, scot, wal, all), nrow = 1, ncol = 4)
  
  colnames(ret) <- c(paste("England (N=", EN, ")", sep = ""),
                     paste("Scotland (N=", SN, ")", sep = ""),
                     paste("Wales (N=", WN, ")", sep = ""),
                     paste("All (N=", AN, ")", sep = ""))
  
  
  return(ret)
}

# And another one that will work for calculatng frequencies:

# Changing this so it's inline with what Sophie wants

myFreqTable <- function(x, varname) {
  
  
  varname <- as.character(varname)
  print(varname)
  gen.E <- x %>% filter(country == "England") %>% dplyr::select(!!varname) %>% drop_na()
  EN <- nrow(gen.E)
  gen.E0 <- as.data.frame(table(gen.E[[1]]))
  gen.E1 <- as.data.frame(round(prop.table(table(gen.E[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.E2 <- inner_join(gen.E0, gen.E1, by = "Var1")
  gen.E2$England <- paste(format(gen.E2$Freq, big.mark=",", trim=TRUE), " (", # N
                          trimws(format(round(gen.E2$perc, 1), nsmall = 1)), "%)", sep = "") # %
  gen.E2 <- select(gen.E2, Var1, England)
  print(gen.E2)
  
  
  gen.W <- x %>% filter(country == "Wales") %>% dplyr::select(!!varname) %>% drop_na()
  WN <- nrow(gen.W)
  gen.W0 <- as.data.frame(table(gen.W[[1]]))
  gen.W1 <- as.data.frame(round(prop.table(table(gen.W[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.W2 <- inner_join(gen.W0, gen.W1, by = "Var1")
  gen.W2$Wales <- paste(format(gen.W2$Freq, big.mark=",", trim=TRUE), " (",
                        trimws(format(round(gen.W2$perc, 1), nsmall = 1)),  "%)", sep = "")
  gen.W2 <- select(gen.W2, Var1, Wales)
  print(gen.W2)
  
  gen.S <- x %>% filter(country == "Scotland") %>% dplyr::select(!!varname) %>% drop_na()
  SN <- nrow(gen.S)
  gen.S0 <- as.data.frame(table(gen.S[[1]]))
  gen.S1 <- as.data.frame(round(prop.table(table(gen.S[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.S2 <- inner_join(gen.S0, gen.S1, by = "Var1")
  gen.S2$Scotland <- paste(format(gen.S2$Freq, big.mark=",", trim=TRUE)," (",
                           trimws(format(round(gen.S2$perc, 1), nsmall = 1)),  "%)", sep = "")
  gen.S2 <- select(gen.S2, Var1, Scotland)
  print(gen.S2)
  
  gen.A <- x %>% dplyr::select(!!varname) %>% drop_na()
  AN <- nrow(gen.A)
  gen.A0 <- as.data.frame(table(gen.A[[1]]))
  gen.A1 <- as.data.frame(round(prop.table(table(gen.A[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.A2 <- inner_join(gen.A0, gen.A1, by = "Var1")
  gen.A2$All <- paste(format(gen.A2$Freq, big.mark=",", trim=TRUE), " (",
                      trimws(format(round(gen.A2$perc, 1), nsmall = 1)),  "%)", sep = "")
  gen.A2 <- select(gen.A2, Var1, All)
  print(gen.A2)
  
  gen.table <- inner_join(gen.E2, gen.S2, by = "Var1") %>% inner_join(gen.W2, by = "Var1") %>%
    inner_join(gen.A2, by = "Var1")
  colnames(gen.table) <- c(varname, paste("England (N=", format(EN, big.mark=",", trim=TRUE), ")", sep = ""),
                           paste("Scotland (N=", format(SN, big.mark=",", trim=TRUE), ")", sep = ""),
                           paste("Wales (N=", format(WN, big.mark=",", trim=TRUE), ")", sep = ""),
                           paste("All (N=", format(AN, big.mark=",", trim=TRUE), ")", sep = ""))
  
  
  
  # row.names(gen.table) <- gen.table$Var1
  
  return(gen.table)
}




histnorm <- function(g) {
  
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g, na.rm = TRUE), max(g, na.rm = TRUE), length = 40) 
  yfit <- dnorm(xfit, mean = mean(g, na.rm = TRUE), sd = sd(g, na.rm = TRUE)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  
  plot(h, ylim = c(0, max(yfit)))
  lines(xfit, yfit, col = "black", lwd = 2)
}


nlc <- function(x) {cat(paste("\n", x, "\n", sep = ""))}
CP <- function(x) {write.table(x, "clipboard", sep = "\t", row.names = FALSE)}
CPwithrn <- function(x) {write.table(x, "clipboard", sep = "\t", row.names = TRUE)}

# When we read this in, if we use stringsasfactors = true then we can say that blank strings count as missing,
# without needing to recode them all. Then we can change all the necessary factors back to characters if need be.

dat <- read.csv("Z:/Group_work/Alex/Encrypted/Alex/PR/PR clinical 2019 interim/Data/rawData/PR-extract-v401-201903-201908/PR-extract-v401.csv", header = TRUE, 
                stringsAsFactors = TRUE, na.strings = c("NA", ""))

# colnames(dat)

# First thing's first - rename column names
# R style guide says variables shouldn't be capitalised.
# NR = not recorded
# init = initial = at assessment
# dis = discharge

dat <- dat %>% rename(study_ID = StudyID,
                      patient_ID = PatientID, 
                      country = Country, 
                      trust_code = TrustCode,
                      org_code = Unit.Org.Code,     
                      nhs_number_valid = NHS.Number.Valid,
                      LSOA = LSOA,
                      age = X.AgeAtAssessment,
                      gender = X1.3.Gender,
                      ethnicity = X1.5.Ethnicity,
                      ref_date = X2.1.Referral.Date,
                      ref_date_NR = X2.1.1.Not.recorded,
                      ref_location = X2.2.Where.was.the.patient.referred.from.,
                      assess_date = X2.3.Initial.PR.Assessment.Appointment,
                      smoke_status = X3.1.Smoking, 
                      FEV1_percpred = X3.2.FEV1...predicted.,
                      FEV1_NR = X3.2.1.Not.recorded,
                      FEV1FVC = X3.3.What.was.the.most.recently.recorded.FEV1.FVC.ratio.,
                      FEV1FVC_NR = X3.3.1.Not.recorded,
                      BMI = X3.4.Patient.s.body.mass.index..BMI.,
                      BMI_NR = X3.4.1.Not.recorded,
                      MRC_score_init = X3.5.MRC.Score,
                      CVD_history_orig = X3.6.Was.a.history.of.cardiovascular.disease.recorded.for.this.patient.,
                      musc_skel_history_orig = X3.7.Was.a.history.of.a.lower.limb.musculoskeletal.disorder.recorded.for.this.patient.,
                      mental_history_orig = X3.8.Was.a.history.of.history.of.mental.illness.recorded.for.this.patient.,
                      mental_history_combined_illness = X3.8a.If..Yes....select.all.mental.health.illnesses.recorded,
                      anxiety_bin = X3.8a.Mental.health...Anxiety,
                      depression_bin = X3.8a.Mental.health...Depression,
                      SMI_bin = X3.8a.Mental.health...Severe.mental.illness,
                      test_type_init = X4.1.Which.test.did.you.record.during.initial.assessment.,
                      test_value_init = X4.1a.Value.in.metres,
                      prac_test_init = X4.1b.Was.a.practice.test.performed.at.the.time.of.the.initial.assessment.,
                      ESWT_at_init = X4.2.Did.you.also.record.the.Endurance.shuttle.walk.test..ESWT..,
                      ESWT_value_init = X4.2a.Value.in.seconds,
                      CRQ_init = X4.3.Chronic.Respiratory.Questionnaire..CRQ.,
                      CRQ_dyspnoea_init = X4.3a.Dyspnoea.score,
                      CRQ_fatigue_init = X4.3b.Fatigue.score,                                                                  
                      CRQ_emotion_init = X4.3c.Emotion.score,                                                                 
                      CRQ_mastery_init = X4.3d.Mastery.score,
                      CAT_init = X4.4.COPD.Assessment.Test..CAT.,
                      CAT_score_init = X4.4a.Total.score,
                      enrolled = X5.1.Post.assessment..was.the.patient.enrolled.onto.a.PR.programme.,
                      start_date = X5.1a.Start.Date,
                      centre_based = X5.2.centre.based.PR.programme,
                      prog_type = X5.2a.Programme.Type,
                      scheduled_sess_centre_no = X5.2b.Centre.based.PR.Sessions.Scheduled,
                      rec_sess_centre_group_no = X5.2c.a..Group.centre.based.sessions.received,
                      rec_sess_centre_indiv_no = X5.2c.b..1.1.centre.based.sessions.received,
                      home_based = X5.3.Home.based.PR.programme,
                      scheduled_sess_home_no = X5.3a.Home.based.PR.Sessions.Scheduled,
                      rec_sess_home_in_person_no = X5.3b.a..Home.based.PR.Sessions...In.Person,
                      rec_sess_home_video_group_no = X5.3b.b..Home.based.PR.Sessions...Video.conferencing...group,
                      rec_sess_home_video_indiv_no = X5.3b.c..Home.based.PR.Sessions...Video.conferencing...1.1,
                      rec_sess_home_phone_no = X5.3b.d..Home.based.PR.Sessions...Phone.calls,
                      rec_sess_home_other_no = X5.3b.e..Home.based.PR.Sessions...Other.Digital.Communications,
                      discharge_assess = X6.1.Was.a.discharge.assessment.performed.,
                      discharge_date = X6.1a.Discharge.assessment.date,
                      exercise_plan = X6.1b.Individualised.discharge.exercise.plan.provided,
                      MRC_score_dis = X7.1.MRC.Score.At.Discharge,
                      test_type_dis = X7.2.Which.walking.test.at.discharge.assessment.,
                      test_value_dis = X7.2a.Value.in.metres,
                      ESWT_at_dis = X7.3.Did.you.also.record.the.ESWT.at.discharge.,
                      ESWT_value_dis = X7.3a.Value.in.seconds,
                      CRQ_dis = X7.4.Chronic.Respiratory.Questionnaire..CRQ.,
                      CRQ_dyspnoea_dis = X7.4a.Dyspnoea.score,
                      CRQ_fatigue_dis = X7.4b.Fatigue.score,
                      CRQ_emotion_dis = X7.4c.Emotion.score,
                      CRQ_mastery_dis = X7.4d.Mastery.score,
                      CAT_dis = X7.5.COPD.Assessment.Test..CAT.,
                      CAT_score_dis = X7.5a.Total.score,
                      dataset_version = Dataset) 

# Phew!

dat <- dat %>% mutate(study_ID = as.character(study_ID),
                      patient_ID = as.character(patient_ID),
                      LSOA = as.character(LSOA))
                      

# It's strange but these do all have different decimal point levels
head(dat$CAT_score_init, 30) # 0 dp
head(dat$CRQ_fatigue_init, 30) # 2 dp
head(dat$CRQ_dyspnoea_init, 30) # 1 dp
head(dat$CRQ_emotion_init, 30) # 2dp
head(dat$CRQ_mastery_init, 30) # 2dp

head(dat$CAT_score_dis, 30) # 0 dp
head(dat$CRQ_fatigue_dis, 30) # 2 dp
head(dat$CRQ_dyspnoea_dis, 30) # 1 dp
head(dat$CRQ_emotion_dis, 30) # 2dp
head(dat$CRQ_mastery_dis, 30) # 2dp


# Everything seems to be alright at this point.

# Now would probably be a good point to link to the LSOAs.

# Read in the IMD quintiles


IMDeng <- read.csv("Z:/Group_work/PS_AA/General UK data/IMD/clean_IMD2019_England.csv",
                   header = TRUE, stringsAsFactors = FALSE)
IMDwales <- read.csv("Z:/Group_work/PS_AA/General UK data/IMD/clean_IMD2019_Wales.csv",
                   header = TRUE, stringsAsFactors = FALSE)
IMDscot <- read.csv("Z:/Group_work/PS_AA/General UK data/IMD/clean_IMD2016_Scotland.csv",
                   header = TRUE, stringsAsFactors = FALSE)


# Create the quintiles for the English IMD data

IMDeng$IMD_quintile <- NA

IMDeng$IMD_quintile[IMDeng$IMD_decile == 1] <- 1
IMDeng$IMD_quintile[IMDeng$IMD_decile == 2] <- 1
IMDeng$IMD_quintile[IMDeng$IMD_decile == 3] <- 2
IMDeng$IMD_quintile[IMDeng$IMD_decile == 4] <- 2
IMDeng$IMD_quintile[IMDeng$IMD_decile == 5] <- 3
IMDeng$IMD_quintile[IMDeng$IMD_decile == 6] <- 3
IMDeng$IMD_quintile[IMDeng$IMD_decile == 7] <- 4
IMDeng$IMD_quintile[IMDeng$IMD_decile == 8] <- 4
IMDeng$IMD_quintile[IMDeng$IMD_decile == 9] <- 5
IMDeng$IMD_quintile[IMDeng$IMD_decile == 10] <- 5


IMDeng <- IMDeng %>% select(LSOA = LSOA_code_2011, IMD_quintile_Eng = IMD_quintile)
IMDwales <- IMDwales %>% select(LSOA = LSOA_Code, IMD_quintile_Wal = WIMD_2019_Overall_Quintile)
IMDscot <- IMDscot %>% select(LSOA = LSOA_Code, IMD_quintile_Scot = IMD_quintile)



# Join them together:

dat <- left_join(dat, IMDeng, by = "LSOA")
dat <- left_join(dat, IMDwales, by = "LSOA")
dat <- left_join(dat, IMDscot, by = "LSOA")


dat <- dat %>% mutate(IMD_quintile_Eng = factor(IMD_quintile_Eng),
                      IMD_quintile_Wal = factor(IMD_quintile_Wal),
                      IMD_quintile_Scot = factor(IMD_quintile_Scot))


# create a variable to represent any IMD present

dat <- dat %>% mutate(anyIMD = factor(ifelse(is.na(IMD_quintile_Eng) & is.na(IMD_quintile_Wal) &
                                        is.na(IMD_quintile_Scot), "No IMD", "IMD present")))

# Issues 

# an 'N/A' option for cvd, musculoskeletal, and mental illness, N = 1-2 for each option
# dat %>% select(CVD_history, musc_skel_history, mental_history) %>% summary()

# inconsistency around home-based and centre-based PR - recode as those marked 'yes' for both as 'mixed modality'?
# N = 28
# dat %>% filter(centre_based == "Yes" & home_based == "Yes") %>% nrow()

# People not marked as yes or no for discharge test type - code as 'no'? N = 87 (for the 3 months worth of data)
# No one who's missing has a value
# dat %>% filter(!is.na(discharge_date) & is.na(test_type_dis)) %>% select(discharge_date, test_type_dis, 
#                                                                          test_value_dis, test_type_init)
# Issue resolved! People who did not perform a test at the initial assessment are not included for this.


# I'd like to rename the country variables:

dat$country <- as.character(dat$country)
dat$country[dat$country == "EN"] <- "England"
dat$country[dat$country == "SC"] <- "Scotland"
dat$country[dat$country == "WA"] <- "Wales"
dat$country <- factor(dat$country, levels = c("England", "Scotland", "Wales"), ordered = FALSE)
levels(dat$country)

# relevel gender

dat$gender <- relevel(dat$gender, ref = "Male")


# dat %>% filter(assess_date < "2019-06-01") %>% filter(home_based == "Yes") %>% select(centre_based:discharge_assess)

# Next thing is to remove the brackets form the study and patient IDs - very unhelpful. 

dat <- dat %>% mutate(study_ID = str_remove_all(study_ID, "[()]"))
dat <- dat %>% mutate(patient_ID = str_remove_all(patient_ID, "[()]"))


#Need to convert dates to dates

dat <- dat %>% mutate(ref_date = as.Date(ref_date, format = "%d/%m/%Y"), 
                      assess_date = as.Date(assess_date, format = "%d/%m/%Y"),
                      start_date = as.Date(start_date, format = "%d/%m/%Y"),
                      discharge_date = as.Date(discharge_date, format = "%d/%m/%Y"))


# Create the 'date to date' variables

dat <- dat %>% mutate(ref_to_start_days = as.numeric(start_date - ref_date),
                      assess_to_start_days = as.numeric(start_date - assess_date),
                      assess_to_discharge_days = as.numeric(discharge_date - assess_date))

# ref to start days needs to be for stable COPD

dat <- dat %>% mutate(ref_to_start_days_stable_COPD = ifelse(
    ref_location == "Primary/Community - stable COPD" |
      ref_location == "Secondary Care - stable COPD", 
    ref_to_start_days, NA))

summary(dat$ref_to_start_days_stable_COPD)


# assess to start days probably also needs to be for stable COPD

dat <- dat %>% mutate(assess_to_start_days_stable_COPD = ifelse(
  ref_location == "Primary/Community - stable COPD" |
    ref_location == "Secondary Care - stable COPD", 
  assess_to_start_days, NA))

summary(dat$assess_to_start_days_stable_COPD)


# Make a variable that says whether or not referal date was recorded
                      
dat <- dat %>% mutate(ref_date_rec = as.character(ref_date_NR))
dat$ref_date_rec[is.na(dat$ref_date_rec)] <- "Known"
dat$ref_date_rec <- factor(dat$ref_date_rec)

# We created a more useful variable so we can drop the other one:

dat$ref_date_NR <- NULL



# Make a variable that says whether or not referal date was recorded

dat <- dat %>% mutate(FEV1_percpred_rec = as.character(FEV1_NR),
                      FEV1FVC_rec = as.character(FEV1FVC_NR),
                      BMI_rec = as.character(BMI_NR))

dat$FEV1_percpred_rec[is.na(dat$FEV1_percpred_rec)] <- "Recorded"
dat$FEV1FVC_rec[is.na(dat$FEV1FVC_rec)] <- "Recorded"
dat$BMI_rec[is.na(dat$BMI_rec)] <- "Recorded"

dat <- dat %>% mutate(FEV1_percpred_rec = factor(FEV1_percpred_rec),
                      FEV1FVC_rec = factor(FEV1FVC_rec),
                      BMI_rec = factor(BMI_rec))

dat$FEV1_NR <- NULL
dat$FEV1FVC_NR <- NULL
dat$BMI_NR <- NULL



# Start date within 90 days for stable COPD

dat <- dat %>% mutate(ninety_day_referral_to_start_for_stable_COPD = ifelse(
                      ref_to_start_days < 90, "<90 days", ">=90 days")) %>%
                        mutate(ninety_day_referral_to_start_for_stable_COPD = factor(ifelse(
                          ref_location == "Primary/Community - stable COPD" |
                            ref_location == "Secondary Care - stable COPD", 
                          ninety_day_referral_to_start_for_stable_COPD, NA)))



dat <- dat %>% mutate(thirty_day_referral_to_start_for_AECOPD = ifelse(
  ref_to_start_days < 30, "<30 days", ">=30 days")) %>%
  mutate(thirty_day_referral_to_start_for_AECOPD = as.factor(ifelse(
    ref_location == "Primary/Community - post treatment for AECOPD" |
      ref_location == "Secondary Care - post admission for AECOPD", 
    thirty_day_referral_to_start_for_AECOPD, NA))) # %>% 

# summary(dat$ninety_day_referral_to_start_for_stable_COPD)
# summary(dat$thirty_day_referral_to_start_for_AECOPD)


# create new health history variables that group 'No' and 'N/A' together

dat <- dat %>% mutate(CVD_history = fct_recode(CVD_history_orig, No = "N/A"),
                      musc_skel_history = fct_recode(musc_skel_history_orig, No = "N/A"),
                      mental_history = fct_recode(mental_history_orig, No = "N/A"),
                      CVD_history_var_for_NA_N_only = fct_recode(CVD_history_orig, `Not N/A` = "No",
                                                                 `Not N/A` = "Yes"),
                      musc_skel_history_var_for_NA_N_only = fct_recode(musc_skel_history_orig, 
                                                                       `Not N/A` = "No", `Not N/A` = "Yes"),
                      mental_history_var_for_NA_N_only = fct_recode(mental_history_orig, 
                                                     `Not N/A` = "No", `Not N/A` = "Yes"))


# recode the mental health variables so they're more easily readable

dat <- dat %>% mutate(anxiety_bin = as.character(anxiety_bin),
                      depression_bin = as.character(depression_bin),
                      SMI_bin = as.character(SMI_bin))

dat$anxiety_bin[dat$anxiety_bin == "0"] <- "No"
dat$anxiety_bin[dat$anxiety_bin == "1"] <- "Yes"
dat$depression_bin[dat$depression_bin == "0"] <- "No"
dat$depression_bin[dat$depression_bin == "1"] <- "Yes"
dat$SMI_bin[dat$SMI_bin == "0"] <- "No"
dat$SMI_bin[dat$SMI_bin == "1"] <- "Yes"

dat <- dat %>% mutate(anxiety_bin = factor(anxiety_bin),
                      depression_bin = factor(depression_bin),
                      SMI_bin = factor(SMI_bin))




# Sort out tests

# Initial tests

dat <- dat %>% mutate(all_3_test_types_init = ifelse(test_type_init == "6MWT" & ESWT_at_init == "No", "6MWT only",
                                              ifelse(test_type_init == "ISWT" & ESWT_at_init == "No", "ISWT only",
                                              ifelse(test_type_init == "None", "None",
                                                     "6MWT or ISWT, and ESWT"))))

dat <- dat %>% mutate(who_only_gets_ESWT_init = ifelse(ESWT_at_init == "No", NA, 
                                                   ifelse(ESWT_at_init == "Yes" & test_type_init == "None", "Only ESWT",
                                                          "ESWT and other walking test")))


dat <- dat %>% mutate(which_other_test_with_ESWT_init = ifelse(ESWT_at_init == "No", NA, 
                                              ifelse(ESWT_at_init == "Yes" & test_type_init == "ISWT", "ISWT",
                                              ifelse(ESWT_at_init == "Yes" & test_type_init == "6MWT", "6MWT", NA))))
 
dat$all_3_test_types_init <- factor(dat$all_3_test_types_init)
dat$who_only_gets_ESWT_init <- factor(dat$who_only_gets_ESWT_init)
dat$which_other_test_with_ESWT_init <- factor(dat$which_other_test_with_ESWT_init)

summary(dat$all_3_test_types_init)
summary(dat$who_only_gets_ESWT_init)
summary(dat$which_other_test_with_ESWT_init)

table(dat$ESWT_at_init, dat$all_3_test_types_init, useNA = "ifany")

# Test value and practice test variable broken down by test type
dat <- dat %>% mutate(test_value_6MWT_init = ifelse(test_type_init == "6MWT", test_value_init, NA),
                      test_value_ISWT_init = ifelse(test_type_init == "ISWT", test_value_init, NA))

dat <- dat %>% mutate(prac_test_6MWT_init = factor(ifelse(test_type_init == "6MWT",
                                                          as.character(prac_test_init), NA)),
                      prac_test_ISWT_init = factor(ifelse(test_type_init == "ISWT",
                                                          as.character(prac_test_init), NA)))
 


# Discharge
dat <- dat %>% mutate(all_3_test_types_dis = ifelse(test_type_dis == "6MWT" & ESWT_at_dis == "No", "6MWT only",
                                                     ifelse(test_type_dis == "ISWT" & ESWT_at_dis == "No", "ISWT only",
                                                            ifelse(test_type_dis == "None", "None",
                                                                   "6MWT or ISWT, and ESWT"))))

dat <- dat %>% mutate(who_is_ESWT_given_to_dis = ifelse(ESWT_at_dis == "Yes" & test_type_dis == "ISWT", "ISWT",
                                                         ifelse(ESWT_at_dis == "Yes" & test_type_dis == "6MWT", "6MWT",
                                                                ifelse(ESWT_at_dis == "Yes" & test_type_dis == "None",
                                                                       "ESWT only", NA))))

# Test value and practice test variable broken down by test type
dat <- dat %>% mutate(test_value_6MWT_dis = ifelse(test_type_dis == "6MWT", test_value_dis, NA),
                      test_value_ISWT_dis = ifelse(test_type_dis == "ISWT", test_value_dis, NA))

# No practise test info for discharge


# Where are patients enrolled?

table(dat$centre_based, dat$home_based, useNA = "ifany")

dat <- dat %>% mutate(PR_location = NA)
dat$PR_location[dat$centre_based == "Yes" & dat$home_based == "No"] <- "Centre-based"
dat$PR_location[dat$centre_based == "No" & dat$home_based == "Yes"] <- "Home-based"
dat$PR_location[dat$centre_based == "Yes" & dat$home_based == "Yes"] <- "Both"

table(dat$PR_location)


summary(dat$prog_type)


# Alternate variable which is a mix of programme type and location

dat <- dat %>% mutate(PR_type_location = NA)
dat$PR_type_location[dat$prog_type == "Rolling" & dat$home_based == "No"] <- "Rolling"
dat$PR_type_location[dat$prog_type == "Cohort" & dat$home_based == "No"] <- "Cohort"
dat$PR_type_location[dat$centre_based == "No" & dat$home_based == "Yes"] <- "Home-based"
dat$PR_type_location[dat$centre_based == "Yes" & dat$home_based == "Yes"] <- "Both"



# table(dat$PR_type_location, useNA = "ifany")
# table(dat$PR_type_location, dat$prog_type, useNA = "ifany")
# table(dat$PR_type_location, dat$centre_based, useNA = "ifany")
# table(dat$PR_type_location, dat$home_based, useNA = "ifany")


# All centre based sessions

dat <- dat %>% mutate(rec_sess_centre_no = rec_sess_centre_group_no + rec_sess_centre_indiv_no)

# All home-based sessions

dat <- dat %>% mutate(rec_sess_home_no = rec_sess_home_in_person_no + rec_sess_home_other_no +
                        rec_sess_home_phone_no + rec_sess_home_video_group_no + 
                        rec_sess_home_video_indiv_no)



# discharge assess bin

summary(dat$discharge_assess)

dat <- dat %>% mutate(discharge_assess_bin = NA)
dat$discharge_assess_bin[dat$discharge_assess == "Yes"] <- "Yes"
dat$discharge_assess_bin[dat$discharge_assess == "No - DNA" | 
                         dat$discharge_assess == "No - drop-out - health reasons" |
                         dat$discharge_assess == "No - drop-out - patient choice" ] <- "No"


table(dat$discharge_assess, dat$discharge_assess_bin)

# Change the discharge assessment variable...

dat <- dat %>% rename(discharge_assess_no_reason = discharge_assess)
dat$discharge_assess_no_reason[dat$discharge_assess_no_reason == "Yes"] <- NA
dat <- dat %>% mutate(discharge_assess_no_reason = fct_drop(discharge_assess_no_reason))

table(dat$discharge_assess_no_reason, dat$discharge_assess_bin, useNA = "ifany")


# Create the variables of discharge assessment by programme type


dat <- dat %>% mutate(discharge_assess_bin_by_rolling = ifelse(PR_type_location == "Rolling",
                                                               discharge_assess_bin, NA))
table(dat$discharge_assess_bin_by_rolling, dat$PR_type_location, useNA = "ifany")


dat <- dat %>% mutate(discharge_assess_bin_by_cohort = ifelse(PR_type_location == "Cohort",
                                                               discharge_assess_bin, NA))
table(dat$discharge_assess_bin_by_cohort, dat$PR_type_location, useNA = "ifany")


dat <- dat %>% mutate(discharge_assess_bin_by_home = ifelse(PR_type_location == "Home-based",
                                                               discharge_assess_bin, NA))
table(dat$discharge_assess_bin_by_home, dat$PR_type_location, useNA = "ifany")


dat <- dat %>% mutate(discharge_assess_bin_by_both = ifelse(PR_type_location == "Both",
                                                               discharge_assess_bin, NA))
table(dat$discharge_assess_bin_by_both, dat$PR_type_location, useNA = "ifany")


# Also create the variable exercise plan by programme type
# To do this, I need to convert exercise plan to a characer variable first

dat <- dat %>% mutate(exercise_plan = as.character(exercise_plan))


dat <- dat %>% mutate(exercise_plan_by_rolling = factor(ifelse(PR_type_location == "Rolling",
                                                               exercise_plan, NA)))
table(dat$exercise_plan_by_rolling, dat$PR_type_location, useNA = "ifany")


dat <- dat %>% mutate(exercise_plan_by_cohort = factor(ifelse(PR_type_location == "Cohort",
                                                              exercise_plan, NA)))
table(dat$exercise_plan_by_cohort, dat$PR_type_location, useNA = "ifany")


dat <- dat %>% mutate(exercise_plan_by_home = factor(ifelse(PR_type_location == "Home-based",
                                                            exercise_plan, NA)))
table(dat$exercise_plan_by_home, dat$PR_type_location, useNA = "ifany")


dat <- dat %>% mutate(exercise_plan_by_both = factor(ifelse(PR_type_location == "Both",
                                                            exercise_plan, NA)))
table(dat$exercise_plan_by_both, dat$PR_type_location, useNA = "ifany")


dat$PR_type_location <- factor(dat$PR_type_location)


sum(table(dat$MRC_score_init, dat$MRC_score_dis)[6, ]) + 
  sum(table(dat$MRC_score_init, dat$MRC_score_dis)[ ,6]) - 
  sum(table(dat$MRC_score_init, dat$MRC_score_dis)[6,6])
  
  

dat <- dat %>% mutate(MRC_score_both = "Both known")
dat$MRC_score_both[dat$MRC_score_init == "Not recorded" | 
                     dat$MRC_score_dis == "Not recorded"] <- "1 or more not known"
dat$MRC_score_both[is.na(dat$MRC_score_init) | is.na(dat$MRC_score_dis)] <- NA


table(dat$MRC_score_both, useNA = "ifany")
table(dat$MRC_score_dis, useNA = "ifany")
table(dat$MRC_score_init, useNA = "ifany")

# Done.

# Now, to create the variable that says whether MRC score has changed.
# to do this safely, need to explicitly give the levels for the MRC score.

levels(dat$MRC_score_init) <- c("Grade 1", "Grade 2", "Grade 3", "Grade 4", "Grade 5", "Not recorded")
levels(dat$MRC_score_dis) <- c("Grade 1", "Grade 2", "Grade 3", "Grade 4", "Grade 5", "Not recorded")

dat <- dat %>% mutate(MRC_score_init_num = as.numeric(MRC_score_init))
dat$MRC_score_init_num[dat$MRC_score_init_num == 6] <- NA

dat <- dat %>% mutate(MRC_score_dis_num = as.numeric(MRC_score_dis))
dat$MRC_score_dis_num[dat$MRC_score_dis_num == 6] <- NA

# This could be useful
dat <- dat %>% mutate(MRC_change_value = MRC_score_dis_num - MRC_score_init_num)

dat <- dat %>% mutate(MRC_change_factor = NA)
dat$MRC_change_factor[dat$MRC_change_value < 0] <- "Improved"
dat$MRC_change_factor[dat$MRC_change_value == 0] <- "Same"
dat$MRC_change_factor[dat$MRC_change_value > 0] <- "Worse"
dat$MRC_change_factor <- factor(dat$MRC_change_factor)



# Create a value for difference in walking tests and helath status questionnaires

dat <- dat %>% mutate(test_value_ISWT_diff = test_value_ISWT_dis - test_value_ISWT_init)
dat <- dat %>% mutate(test_value_6MWT_diff = test_value_6MWT_dis - test_value_6MWT_init)
dat <- dat %>% mutate(test_value_ESWT_diff = ESWT_value_dis - ESWT_value_init)
dat <- dat %>% mutate(CAT_score_diff = CAT_score_dis - CAT_score_init)
dat <- dat %>% mutate(CRQ_dyspnoea_diff = CRQ_dyspnoea_dis - CRQ_dyspnoea_init)
dat <- dat %>% mutate(CRQ_fatigue_diff = CRQ_fatigue_dis - CRQ_fatigue_init)
dat <- dat %>% mutate(CRQ_emotion_diff = CRQ_emotion_dis - CRQ_emotion_init)
dat <- dat %>% mutate(CRQ_mastery_diff = CRQ_mastery_dis - CRQ_mastery_init)


# Create MCID binary values

# MCID is 48 for ISWT
dat <- dat %>% mutate(MCID_ISWT = NA)
dat$MCID_ISWT[dat$test_value_ISWT_diff < 48] <- "MCID not met"
dat$MCID_ISWT[dat$test_value_ISWT_diff >= 48] <- "MCID met"
dat$MCID_ISWT <- factor(dat$MCID_ISWT)
summary(dat$MCID_ISWT)

# MCID is 30 for 6MWT
dat <- dat %>% mutate(MCID_6MWT = NA)
dat$MCID_6MWT[dat$test_value_6MWT_diff < 30] <- "MCID not met"
dat$MCID_6MWT[dat$test_value_6MWT_diff >= 30] <- "MCID met"
dat$MCID_6MWT <- factor(dat$MCID_6MWT)
summary(dat$MCID_6MWT)

dat <- dat %>% mutate(MCID_N_ISWT_6MWT = NA)
dat$MCID_N_ISWT_6MWT[!is.na(dat$MCID_6MWT)] <- 1
dat$MCID_N_ISWT_6MWT[!is.na(dat$MCID_ISWT)] <- 1


# MCID for CAT is -2
dat <- dat %>% mutate(MCID_CAT = NA)
dat$MCID_CAT[dat$CAT_score_diff <= -2] <- "MCID met"
dat$MCID_CAT[dat$CAT_score_diff > -2] <- "MCID not met"
dat$MCID_CAT <- factor(dat$MCID_CAT)
summary(dat$MCID_CAT)

# MCID for each CRQ is 0.5
dat <- dat %>% mutate(MCID_CRQ_dyspnoea = NA)
dat$MCID_CRQ_dyspnoea[dat$CRQ_dyspnoea_diff < 0.5] <- "MCID not met"
dat$MCID_CRQ_dyspnoea[dat$CRQ_dyspnoea_diff >= 0.5] <- "MCID met"
dat$MCID_CRQ_dyspnoea <- factor(dat$MCID_CRQ_dyspnoea)
summary(dat$MCID_CRQ_dyspnoea)

dat <- dat %>% mutate(MCID_CRQ_fatigue = NA)
dat$MCID_CRQ_fatigue[dat$CRQ_fatigue_diff < 0.5] <- "MCID not met"
dat$MCID_CRQ_fatigue[dat$CRQ_fatigue_diff >= 0.5] <- "MCID met"
dat$MCID_CRQ_fatigue <- factor(dat$MCID_CRQ_fatigue)
summary(dat$MCID_CRQ_fatigue)

dat <- dat %>% mutate(MCID_CRQ_emotion = NA)
dat$MCID_CRQ_emotion[dat$CRQ_emotion_diff < 0.5] <- "MCID not met"
dat$MCID_CRQ_emotion[dat$CRQ_emotion_diff >= 0.5] <- "MCID met"
dat$MCID_CRQ_emotion <- factor(dat$MCID_CRQ_emotion)
summary(dat$MCID_CRQ_emotion)

dat <- dat %>% mutate(MCID_CRQ_mastery = NA)
dat$MCID_CRQ_mastery[dat$CRQ_mastery_diff < 0.5] <- "MCID not met"
dat$MCID_CRQ_mastery[dat$CRQ_mastery_diff >= 0.5] <- "MCID met"
dat$MCID_CRQ_mastery <- factor(dat$MCID_CRQ_mastery)
summary(dat$MCID_CRQ_mastery)


dat <- dat %>% mutate(MCID_N_CAT_CRQ = NA)
dat$MCID_N_CAT_CRQ[!is.na(dat$MCID_CAT)] <- 1
dat$MCID_N_CAT_CRQ[!is.na(dat$MCID_CRQ_fatigue)] <- 1



dat <- dat %>% mutate(CAT_CRQ_dis_N = NA)
dat$CAT_CRQ_dis_N[!is.na(dat$CAT_dis)] <- 1
dat$CAT_CRQ_dis_N[!is.na(dat$CRQ_dis)] <- 1


dat <- dat %>% mutate(CAT_CRQ_score_diff_N = NA)
dat$CAT_CRQ_score_diff_N[!is.na(dat$CAT_score_diff)] <- 1
dat$CAT_CRQ_score_diff_N[!is.na(dat$CRQ_score_diff)] <- 1


# recode as factors the variables we left before
dat$discharge_assess_bin <- factor(dat$discharge_assess_bin)

dat$ref_to_start_days_stable_COPD
# Let's create some benchmarking variables

dat <- dat %>% mutate(BM_start_90 = ifelse(ref_to_start_days_stable_COPD < 90, 1,
                                    ifelse(ref_to_start_days_stable_COPD >= 90, 0, NA)))




dat <- dat %>% mutate(BM_prac_test = ifelse(prac_test_init == "No", 0,
                                     ifelse(prac_test_init == "Yes", 1, NA))) 

dat <- dat %>% mutate(BM_discharge_assess = ifelse(discharge_assess_bin == "No", 0,
                                            ifelse(discharge_assess_bin == "Yes", 1, NA)))


dat <- dat %>% mutate(BM_exercise_plan = ifelse(exercise_plan == "No", 0,
                                         ifelse(exercise_plan == "Yes", 1, NA)))


dat <- dat %>% mutate(BM_MCID_exercise = 0)
dat$BM_MCID_exercise[is.na(dat$MCID_6MWT) & is.na(dat$MCID_ISWT)] <- NA
dat$BM_MCID_exercise[dat$MCID_6MWT == "MCID met"] <- 1
dat$BM_MCID_exercise[dat$MCID_ISWT == "MCID met"] <- 1



# table(dat$BM_MCID_exercise, dat$MCID_6MWT, useNA = "ifany")
# table(dat$BM_MCID_exercise, dat$MCID_ISWT, useNA = "ifany")
# table(dat$MCID_ISWT, dat$MCID_6MWT, useNA = "ifany")
# table(dat$BM_MCID_exercise, useNA = "ifany")


dat <- dat %>% mutate(BM_MCID_CAT_CRQ = 0)
dat$BM_MCID_CAT_CRQ[is.na(dat$MCID_CAT) & is.na(dat$MCID_CRQ_dyspnoea)] <- NA

dat$BM_MCID_CAT_CRQ[dat$MCID_CAT == "MCID met"] <- 1
dat$BM_MCID_CAT_CRQ[dat$MCID_CRQ_dyspnoea == "MCID met"] <- 1
dat$BM_MCID_CAT_CRQ[dat$MCID_CRQ_emotion == "MCID met"] <- 1
dat$BM_MCID_CAT_CRQ[dat$MCID_CRQ_fatigue == "MCID met"] <- 1
dat$BM_MCID_CAT_CRQ[dat$MCID_CRQ_mastery == "MCID met"] <- 1


table(dat$BM_MCID_CAT_CRQ, dat$MCID_CAT)
table(dat$BM_MCID_CAT_CRQ, dat$MCID_CRQ_fatigue, useNA = "ifany")
table(dat$BM_MCID_CAT_CRQ, dat$MCID_CRQ_emotion, useNA = "ifany")
table(dat$BM_MCID_CAT_CRQ, dat$MCID_CRQ_dyspnoea, useNA = "ifany")
table(dat$BM_MCID_CAT_CRQ, dat$MCID_CRQ_mastery, useNA = "ifany")
table(dat$BM_MCID_CAT_CRQ, useNA = "ifany")


# We need this for the analysis as well... But I will do it as a factor
dat <- dat %>% mutate(MCID_exercise_cat = BM_MCID_exercise)
dat$MCID_exercise_cat[dat$MCID_exercise_cat == 0] <- "MCID not met"
dat$MCID_exercise_cat[dat$MCID_exercise_cat == 1] <- "MCID met"
dat$MCID_exercise_cat <- factor(dat$MCID_exercise_cat)
dat$MCID_exercise_cat <- relevel(dat$MCID_exercise_cat, ref = "MCID not met")
summary(dat$MCID_exercise_cat)


# We need this for the analysis as well... But I will do it as a factor
dat <- dat %>% mutate(MCID_CAT_CRQ_cat = BM_MCID_CAT_CRQ)
dat$MCID_CAT_CRQ_cat[dat$MCID_CAT_CRQ_cat == 0] <- "MCID not met"
dat$MCID_CAT_CRQ_cat[dat$MCID_CAT_CRQ_cat == 1] <- "MCID met"
dat$MCID_CAT_CRQ_cat <- factor(dat$MCID_CAT_CRQ_cat)
dat$MCID_CAT_CRQ_cat <- relevel(dat$MCID_CAT_CRQ_cat, ref = "MCID not met")
summary(dat$MCID_CAT_CRQ_cat)



# Now let's create some analysis variables



dat <- dat %>% mutate(IMD_quintile_all = IMD_quintile_Eng)
dat$IMD_quintile_all[is.na(dat$IMD_quintile_all)] <- dat$IMD_quintile_Wal[is.na(dat$IMD_quintile_all)]
dat$IMD_quintile_all[is.na(dat$IMD_quintile_all)] <- dat$IMD_quintile_Scot[is.na(dat$IMD_quintile_all)]
dat$IMD_quintile_all <- fct_explicit_na(dat$IMD_quintile_all)
summary(dat$IMD_quintile_all)

dat %>% select(discharge_assess_bin, CAT_init) %>% table()
dat <- dat %>% mutate(agecat = cut(age, breaks = c(35, 45, 55, 65, 75, 85, 105), right = FALSE,
               labels = c("35-44", "45-54", "55-64", "65-74", "75-84", "85+"))) %>%
               mutate(agecat = relevel(agecat, ref = "65-74"))

summary(dat$agecat)





dat$CVD_history <- relevel(dat$CVD_history, ref = "No")
dat$musc_skel_history <- relevel(dat$musc_skel_history, ref = "No")

dat$SMI_bin_no_missing <- as.character(dat$SMI_bin)
dat$SMI_bin_no_missing[is.na(dat$SMI_bin_no_missing)] <- as.character(
                                                         dat$mental_history[is.na(dat$SMI_bin_no_missing)])
dat$SMI_bin_no_missing <- factor(dat$SMI_bin_no_missing)
dat$SMI_bin_no_missing <- relevel(dat$SMI_bin_no_missing, ref = "No")
summary(dat$SMI_bin_no_missing)

dat$anxiety_bin_no_missing <- as.character(dat$anxiety_bin)
dat$anxiety_bin_no_missing[is.na(dat$anxiety_bin_no_missing)] <- as.character(
  dat$mental_history[is.na(dat$anxiety_bin_no_missing)])
dat$anxiety_bin_no_missing <- factor(dat$anxiety_bin_no_missing)
dat$anxiety_bin_no_missing <- relevel(dat$anxiety_bin_no_missing, ref = "No")

dat$depression_bin_no_missing <- as.character(dat$depression_bin)
dat$depression_bin_no_missing[is.na(dat$depression_bin_no_missing)] <- as.character(
  dat$mental_history[is.na(dat$depression_bin_no_missing)])
dat$depression_bin_no_missing <- factor(dat$depression_bin_no_missing)
dat$depression_bin_no_missing <- relevel(dat$depression_bin_no_missing, ref = "No")


summary(dat$depression_bin_no_missing)
summary(dat$SMI_bin_no_missing)
summary(dat$anxiety_bin_no_missing)

dat <- dat %>% mutate(CAT_score_init_cat = cut(CAT_score_init, breaks = c(-1, 10, 20, 30, 40),
                                               labels = c("0-10", "11-20", "21-30", "31-40")))

dat$CAT_score_init_cat <- fct_explicit_na(dat$CAT_score_init_cat, na_level = "CAT not taken")
dat$CAT_score_init_cat <- relevel(dat$CAT_score_init_cat, ref = "21-30")





# plot stuff

dat <- dat %>% mutate(ref_to_start_days_rolling = ifelse(prog_type == "Rolling",
                                                         ref_to_start_days, NA),
                      ref_to_start_days_cohort = ifelse(prog_type == "Cohort",
                                                         ref_to_start_days, NA))
                      

dat <- dat %>% mutate(rec_sess_centre_group_bin = factor(ifelse(is.na(rec_sess_centre_group_no), NA,
                                                  ifelse(rec_sess_centre_group_no == 0, "No", "Yes"))),
                      rec_sess_centre_indiv_bin = factor(ifelse(is.na(rec_sess_centre_indiv_no), NA,
                                                                ifelse(rec_sess_centre_indiv_no == 0, "No", "Yes"))),
                      rec_sess_home_in_person_bin = factor(ifelse(is.na(rec_sess_home_in_person_no), NA,
                                                                ifelse(rec_sess_home_in_person_no == 0, "No", "Yes"))),
                      rec_sess_home_video_group_bin = factor(ifelse(is.na(rec_sess_home_video_group_no), NA,
                                                                ifelse(rec_sess_home_video_group_no == 0, "No", "Yes"))),
                      rec_sess_home_video_indiv_bin = factor(ifelse(is.na(rec_sess_home_video_indiv_no), NA,
                                                                ifelse(rec_sess_home_video_indiv_no == 0, "No", "Yes"))),
                      rec_sess_home_phone_bin = factor(ifelse(is.na(rec_sess_home_phone_no), NA,
                                                                ifelse(rec_sess_home_phone_no == 0, "No", "Yes"))),
                      rec_sess_home_other_bin = factor(ifelse(is.na(rec_sess_home_other_no), NA,
                                                                ifelse(rec_sess_home_other_no == 0, "No", "Yes"))))



# Make separate variables, some with both included, some with both not included, some with only both.
# From now on, the default is to exclude 'both', but I've also created two options for when both
# is included and for when we only want to look at what both is doing.

dat <- dat %>% mutate(prog_type_inc_both = prog_type,
                      scheduled_sess_centre_no_inc_both = scheduled_sess_centre_no,
                      rec_sess_centre_group_no_inc_both = rec_sess_centre_group_no,
                      rec_sess_centre_indiv_no_inc_both = rec_sess_centre_indiv_no,
                      home_based_inc_both = home_based,
                      scheduled_sess_home_no_inc_both = scheduled_sess_home_no,
                      rec_sess_home_in_person_no_inc_both = rec_sess_home_in_person_no,
                      rec_sess_home_video_group_no_inc_both = rec_sess_home_video_group_no,
                      rec_sess_home_video_indiv_no_inc_both = rec_sess_home_video_indiv_no,
                      rec_sess_home_phone_no_inc_both = rec_sess_home_phone_no,
                      rec_sess_home_other_no_inc_both = rec_sess_home_other_no,
                      rec_sess_centre_no_inc_both = rec_sess_centre_no,
                      rec_sess_home_no_inc_both = rec_sess_home_no,
                      
                      rec_sess_centre_group_bin_inc_both = rec_sess_centre_group_bin,
                      rec_sess_centre_indiv_bin_inc_both = rec_sess_centre_indiv_bin,
                      rec_sess_home_in_person_bin_inc_both = rec_sess_home_in_person_bin,
                      rec_sess_home_video_group_bin_inc_both = rec_sess_home_video_group_bin,
                      rec_sess_home_video_indiv_bin_inc_both = rec_sess_home_video_indiv_bin,
                      rec_sess_home_phone_bin_inc_both = rec_sess_home_phone_bin,
                      rec_sess_home_other_bin_inc_both = rec_sess_home_other_bin,
                      
                      prog_type_only_both = prog_type,
                      scheduled_sess_centre_no_only_both = scheduled_sess_centre_no,
                      rec_sess_centre_group_no_only_both = rec_sess_centre_group_no,
                      rec_sess_centre_indiv_no_only_both = rec_sess_centre_indiv_no,
                      home_based_only_both = home_based,
                      scheduled_sess_home_no_only_both = scheduled_sess_home_no,
                      rec_sess_home_in_person_no_only_both = rec_sess_home_in_person_no,
                      rec_sess_home_video_group_no_only_both = rec_sess_home_video_group_no,
                      rec_sess_home_video_indiv_no_only_both = rec_sess_home_video_indiv_no,
                      rec_sess_home_phone_no_only_both = rec_sess_home_phone_no,
                      rec_sess_home_other_no_only_both = rec_sess_home_other_no,
                      rec_sess_centre_no_only_both = rec_sess_centre_no,
                      rec_sess_home_no_only_both = rec_sess_home_no,
                      
                      rec_sess_centre_group_bin_only_both = rec_sess_centre_group_bin,
                      rec_sess_centre_indiv_bin_only_both = rec_sess_centre_indiv_bin,
                      rec_sess_home_in_person_bin_only_both = rec_sess_home_in_person_bin,
                      rec_sess_home_video_group_bin_only_both = rec_sess_home_video_group_bin,
                      rec_sess_home_video_indiv_bin_only_both = rec_sess_home_video_indiv_bin,
                      rec_sess_home_phone_bin_only_both = rec_sess_home_phone_bin,
                      rec_sess_home_other_bin_only_both = rec_sess_home_other_bin)

dat$prog_type[dat$PR_location == "Both"] <- NA
dat$scheduled_sess_centre_no[dat$PR_location == "Both"] <- NA
dat$rec_sess_centre_group_no[dat$PR_location == "Both"] <- NA
dat$rec_sess_centre_indiv_no[dat$PR_location == "Both"] <- NA
dat$home_based[dat$PR_location == "Both"] <- NA
dat$scheduled_sess_home_no[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_in_person_no[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_video_group_no[dat$PR_location == "Both"] <- NA 
dat$rec_sess_home_video_indiv_no[dat$PR_location == "Both"] <- NA 
dat$rec_sess_home_phone_no[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_other_no[dat$PR_location == "Both"] <- NA
dat$rec_sess_centre_no[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_no[dat$PR_location == "Both"] <- NA

dat$rec_sess_centre_group_bin[dat$PR_location == "Both"] <- NA
dat$rec_sess_centre_indiv_bin[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_in_person_bin[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_video_group_bin[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_video_indiv_bin[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_phone_bin[dat$PR_location == "Both"] <- NA
dat$rec_sess_home_other_bin[dat$PR_location == "Both"] <- NA


dat$prog_type_only_both[dat$PR_location != "Both"] <- NA
dat$scheduled_sess_centre_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_centre_group_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_centre_indiv_no_only_both[dat$PR_location != "Both"] <- NA
dat$home_based_only_both[dat$PR_location != "Both"] <- NA 
dat$scheduled_sess_home_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_in_person_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_video_group_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_video_indiv_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_phone_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_other_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_centre_no_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_no_only_both[dat$PR_location != "Both"] <- NA

dat$rec_sess_centre_group_bin_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_centre_indiv_bin_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_in_person_bin_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_video_group_bin_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_video_indiv_bin_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_phone_bin_only_both[dat$PR_location != "Both"] <- NA
dat$rec_sess_home_other_bin_only_both[dat$PR_location != "Both"] <- NA

# Now change all the session count variables so that '0' received sessions is removed

dat$rec_sess_centre_group_no_inc_both[dat$rec_sess_centre_group_no_inc_both == 0] <- NA
dat$rec_sess_centre_indiv_no_inc_both[dat$rec_sess_centre_indiv_no_inc_both == 0] <- NA
dat$rec_sess_home_in_person_no_inc_both[dat$rec_sess_home_in_person_no_inc_both == 0] <- NA
dat$rec_sess_home_video_group_no_inc_both[dat$rec_sess_home_video_group_no_inc_both == 0] <- NA
dat$rec_sess_home_video_indiv_no_inc_both[dat$rec_sess_home_video_indiv_no_inc_both == 0] <- NA
dat$rec_sess_home_phone_no_inc_both[dat$rec_sess_home_phone_no_inc_both == 0] <- NA
dat$rec_sess_home_other_no_inc_both[dat$rec_sess_home_other_no_inc_both == 0] <- NA
# dat$rec_sess_centre_no_inc_both[dat$rec_sess_centre_no_inc_both == 0] <- NA # Not for total variables though
# dat$rec_sess_home_no_inc_both[dat$rec_sess_home_no_inc_both == 0] <- NA

dat$rec_sess_centre_group_no[dat$rec_sess_centre_group_no == 0] <- NA
dat$rec_sess_centre_indiv_no[dat$rec_sess_centre_indiv_no == 0] <- NA
dat$rec_sess_home_in_person_no[dat$rec_sess_home_in_person_no == 0] <- NA
dat$rec_sess_home_video_group_no[dat$rec_sess_home_video_group_no == 0] <- NA
dat$rec_sess_home_video_indiv_no[dat$rec_sess_home_video_indiv_no == 0] <- NA
dat$rec_sess_home_phone_no[dat$rec_sess_home_phone_no == 0] <- NA
dat$rec_sess_home_other_no[dat$rec_sess_home_other_no == 0] <- NA
# dat$rec_sess_centre_no[dat$rec_sess_centre_no == 0] <- NA
# dat$rec_sess_home_no[dat$rec_sess_home_no == 0] <- NA

dat$rec_sess_centre_group_no_only_both[dat$rec_sess_centre_group_no_only_both == 0] <- NA
dat$rec_sess_centre_indiv_no_only_both[dat$rec_sess_centre_indiv_no_only_both == 0] <- NA
dat$rec_sess_home_in_person_no_only_both[dat$rec_sess_home_in_person_no_only_both == 0] <- NA
dat$rec_sess_home_video_group_no_only_both[dat$rec_sess_home_video_group_no_only_both == 0] <- NA
dat$rec_sess_home_video_indiv_no_only_both[dat$rec_sess_home_video_indiv_no_only_both == 0] <- NA
dat$rec_sess_home_phone_no_only_both[dat$rec_sess_home_phone_no_only_both == 0] <- NA
dat$rec_sess_home_other_no_only_both[dat$rec_sess_home_other_no_only_both == 0] <- NA
# dat$rec_sess_centre_no_only_both[dat$rec_sess_centre_no_only_both == 0] <- NA
# dat$rec_sess_home_no_only_both[dat$rec_sess_home_no_only_both == 0] <- NA




# Analyses we're running:
# Likelihood of discharge assessment
# Likelihood of making MCID walking test
# Likelihood of making health status MCID




# Don't need this until further notice

# dat <- dat %>% mutate(CAT_and_CRQ_dis = NA)
# dat$CAT_and_CRQ_dis[dat$CAT_dis == "Yes" & dat$CRQ_dis == "Yes"] <- "Yes"
# dat$CAT_and_CRQ_dis[dat$CAT_dis == "Yes" & dat$CRQ_dis == "No"] <- "No"
# dat$CAT_and_CRQ_dis[dat$CAT_dis == "No" & dat$CRQ_dis == "Yes"] <- "No"
# dat$CAT_and_CRQ_dis[dat$CAT_dis == "No" & dat$CRQ_dis == "No"] <- "No"
# dat$CAT_and_CRQ_dis <- factor(dat$CAT_and_CRQ_dis)
# 
# dat <- dat %>% mutate(CAT_and_CRQ_init = NA)
# dat$CAT_and_CRQ_init[dat$CAT_init == "Yes" & dat$CRQ_init == "Yes"] <- "Yes"
# dat$CAT_and_CRQ_init[dat$CAT_init == "Yes" & dat$CRQ_init == "No"] <- "Just CAT"
# dat$CAT_and_CRQ_init[dat$CAT_init == "No" & dat$CRQ_init == "Yes"] <- "Just CRQ"
# dat$CAT_and_CRQ_init[dat$CAT_init == "No" & dat$CRQ_init == "No"] <- "Neither test"
# dat$CAT_and_CRQ_init <- factor(dat$CAT_and_CRQ_init)






# I've checked that all these variables are okay, and correspond to the number of missing values.

summary(dat$ref_location)

# Do some data cleaning

nlc("Total number of admissions in dataset:")

nrow(dat)

summary(dat$assess_date)

# We are changing the time period to a 3-month time period, from 1st March to 31st May.

nlc("How many people are in the dataset who were assessed before 1st March?")

dat %>% filter(assess_date < "2019-03-01") %>% nrow()

nlc("How many people are assessed after the 31st May?")

dat %>% filter(assess_date > "2019-05-31") %>% nrow()


nlc("We filter these out and are left with this many in the dataset who were assessed between
1st March and 31st May:")

dat <- dat %>% filter(assess_date >= "2019-03-01") %>% filter(assess_date <= "2019-05-31") 

nrow(dat)

summary(dat$assess_date)

nlc("Now we do some data cleaning. Is there anyone who receives their assessment date before their referral date?")

dat %>% filter(assess_date < ref_date) %>% nrow()
dat <- dat %>% filter(assess_date >= ref_date | is.na(assess_date) | is.na(ref_date)) 

nlc("does anyone receive their start date before their assessment date?")

dat %>% filter(start_date < assess_date) %>% nrow()
dat <- dat %>% filter(start_date >= assess_date | is.na(assess_date) | is.na(start_date)) 

nlc("Is anyone discharged before their start date?")
dat %>% filter(discharge_date < start_date) %>% nrow()
dat <- dat %>% filter(discharge_date >= start_date | is.na(discharge_date) | is.na(start_date)) 


nlc("Is anyone marked as missing something who has that same thing encoded?")

# referral date
dat %>% filter(is.na(ref_date) & ref_date_rec == "Known") %>% nrow()
dat %>% filter(is.na(ref_date) & ref_date_rec == "Not known") %>% nrow()
dat %>% filter(!is.na(ref_date) & ref_date_rec == "Known") %>% nrow()
dat %>% filter(!is.na(ref_date) & ref_date_rec == "Not known") %>% nrow()

dat <- dat %>% filter(!(is.na(ref_date) & ref_date_rec == "Known"))
dat <- dat %>% filter(!(!is.na(ref_date) & ref_date_rec == "Not known"))

# This is all fine.

# FEV1 perc pred
dat %>% filter(is.na(FEV1_percpred) & FEV1_percpred_rec == "Recorded") %>% nrow()
dat %>% filter(is.na(FEV1_percpred) & FEV1_percpred_rec == "Not recorded") %>% nrow()
dat %>% filter(!is.na(FEV1_percpred) & FEV1_percpred_rec == "Recorded") %>% nrow()
dat %>% filter(!is.na(FEV1_percpred) & FEV1_percpred_rec == "Not recorded") %>% nrow()

dat <- dat %>% filter(!(is.na(FEV1_percpred) & FEV1_percpred_rec == "Recorded"))
dat <- dat %>% filter(!(!is.na(FEV1_percpred) & FEV1_percpred_rec == "Not recorded"))

# FEV1FVC
dat %>% filter(is.na(FEV1FVC) & FEV1FVC_rec == "Recorded") %>% nrow()
dat %>% filter(is.na(FEV1FVC) & FEV1FVC_rec == "Not recorded") %>% nrow()
dat %>% filter(!is.na(FEV1FVC) & FEV1FVC_rec == "Recorded") %>% nrow()
dat %>% filter(!is.na(FEV1FVC) & FEV1FVC_rec == "Not recorded") %>% nrow()

dat <- dat %>% filter(!(is.na(FEV1FVC) & FEV1FVC_rec == "Recorded"))
dat <- dat %>% filter(!(!is.na(FEV1FVC) & FEV1FVC_rec == "Not recorded"))


# bmi
dat %>% filter(is.na(BMI) & BMI_rec == "Recorded") %>% nrow()
dat %>% filter(is.na(BMI) & BMI_rec == "Not recorded") %>% nrow()
dat %>% filter(!is.na(BMI) & BMI_rec == "Recorded") %>% nrow()
dat %>% filter(!is.na(BMI) & BMI_rec == "Not recorded") %>% nrow()

dat <- dat %>% filter(!(is.na(BMI) & BMI_rec == "Recorded"))
dat <- dat %>% filter(!(!is.na(BMI) & BMI_rec == "Not recorded"))

# No conflicting records.

nlc("How many people are missing IMD quintiles for any country?")
dat %>% filter(is.na(IMD_quintile_Eng) & is.na(IMD_quintile_Scot) & is.na(IMD_quintile_Wal)) %>% nrow()

nlc("This doesn't matter though.")

nlc("Finally, we filter out those with invalid NHS numbers:")
dat %>% filter(nhs_number_valid == 0) %>% nrow()

nlc("Removed to leave this many people in our dataset:")
dat <- dat %>% filter(nhs_number_valid == 1)
nrow(dat)


dat %>% select(country, trust_code, org_code, LSOA, age, gender, ethnicity, ref_date, 
               ref_location, assess_date) %>% nrow()

nlc("assess whether there are duplicate records. Done based on:
country, trust_code, patient_ID, org_code, LSOA, age, gender, ethnicity, ref_date, 
ref_location, assess_date. After duplicates removed, we have this many people:")


dat %>% select(country, trust_code, patient_ID, org_code, LSOA, age, gender, ethnicity, ref_date, 
               ref_location, assess_date) %>% duplicated() %>% sum()


dat <- dat[!duplicated(select(dat, country, trust_code, patient_ID, org_code, LSOA, age, gender, ethnicity, ref_date, 
                               ref_location, assess_date)), ]


nlc("This leaves us with this many people:")
nrow(dat)



# We have a few columns that are character when they should be factor, so we convert them to factor,
# but we leave 'study_ID', 'patient_ID' and 'LSOA' as character.

dat %>% select_if(is.character) %>% colnames()

dat <- dat %>% mutate_if(is.character, factor) %>% 
           mutate_at(c("study_ID", "patient_ID", "LSOA"), as.character)



# For the report...


medTable(dat, "age") %>% CP() # Normal-ish
myFreqTable(dat, "gender") %>% CP()
myFreqTable(dat, "IMD_quintile_Eng") %>% CP()
myFreqTable(dat, "IMD_quintile_Wal") %>% CP()
myFreqTable(dat, "IMD_quintile_Scot") %>% CP()
myFreqTable(dat, "anyIMD") %>% CP()
# myFreqTable(dat, "ethnicity") %>% CP() removed
myFreqTable(dat, "ref_location") %>% CP()
medTable(dat, "ref_to_start_days_stable_COPD") %>% CP() # Not normal

# !!!!!! must be updated in the csv version ^^^^^^^^^

myFreqTable(dat, "ninety_day_referral_to_start_for_stable_COPD") %>% CP()
myFreqTable(dat, "thirty_day_referral_to_start_for_AECOPD") %>% CP()
# medTable(dat, "assess_to_start_days") %>% CP() # Not normal
medTable(dat, "assess_to_start_days_stable_COPD")

summary(dat$ref_date)
summary(dat$assess_date)
summary(dat$start_date)


myFreqTable(dat, "smoke_status") %>% CP()
myFreqTable(dat, "FEV1_percpred_rec") %>% CP()
medTable(dat, "FEV1_percpred") %>% CP() # Not normal
myFreqTable(dat, "FEV1FVC_rec") %>% CP()
medTable(dat, "FEV1FVC") %>% CP() # Not normal
myFreqTable(dat, "BMI_rec") %>% CP()
medTable(dat, "BMI") %>% CP() # Not normal
myFreqTable(dat, "MRC_score_init") %>% CP()
myFreqTable(dat, "CVD_history") %>% CP()
myFreqTable(dat, "CVD_history_var_for_NA_N_only") %>% CP()
myFreqTable(dat, "musc_skel_history") %>% CP()
myFreqTable(dat, "musc_skel_history_var_for_NA_N_only") %>% CP()
myFreqTable(dat, "mental_history") %>% CP()
myFreqTable(dat, "mental_history_var_for_NA_N_only") %>% CP()
myFreqTable(dat, "anxiety_bin") %>% CP()
myFreqTable(dat, "depression_bin") %>% CP()
myFreqTable(dat, "SMI_bin") %>% CP()

dat %>% select(country, all_3_test_types_init) %>% table()
dat %>% select(test_type_init) %>% table()


myFreqTable(dat, "all_3_test_types_init") %>% CP() # Needs extra work

FreqSum(dat, "who_is_ESWT_given_to_init") %>% CP()

str(dat$prac_test_ISWT_init)

medTable(dat, "test_value_ISWT_init") %>% CP() # Not normal
myFreqTable(dat, "prac_test_ISWT_init") %>% CP()
medTable(dat, "test_value_6MWT_init") %>% CP() # Not normal
myFreqTable(dat, "prac_test_6MWT_init") %>% CP()
medTable(dat, "ESWT_value_init") %>% CP() # Not normal
myFreqTable(dat, "who_only_gets_ESWT_init") %>% CP()
myFreqTable(dat, "which_other_test_with_ESWT_init") %>% CP()

myFreqTable(dat, "CAT_init") %>% CP()
myFreqTable(dat, "CRQ_init") %>% CP()
medTable(dat, "CAT_score_init") %>% CP() # Reasonably normal
medTable(dat, "CRQ_dyspnoea_init") %>% CP() # Not normal
medTable(dat, "CRQ_fatigue_init") %>% CP() # 
medTable(dat, "CRQ_emotion_init") %>% CP() # 
medTable(dat, "CRQ_mastery_init") %>% CP()
myFreqTable(dat, "enrolled") %>% CP()   # seems like a lot of issues are getting patients to turn up 
# to assessment maybe?
myFreqTable(dat, "PR_location") %>% CP()
myFreqTable(dat, "prog_type") %>% CP()
medTable(dat, "scheduled_sess_centre_no") %>% CP()
myFreqTable(dat, "rec_sess_centre_group_bin") %>% CP()
myFreqTable(dat, "rec_sess_centre_indiv_bin") %>% CP()
medTable(dat, "rec_sess_centre_group_no") %>% CP()
medTable(dat, "rec_sess_centre_indiv_no") %>% CP()
medTable(dat, "rec_sess_centre_no") %>% CP()



medTable(dat, "scheduled_sess_home_no") %>% CP()
myFreqTable(dat, "rec_sess_home_in_person_bin") %>% CP()
myFreqTable(dat, "rec_sess_home_video_group_bin") %>% CP()
myFreqTable(dat, "rec_sess_home_video_indiv_bin") %>% CP()
myFreqTable(dat, "rec_sess_home_phone_bin") %>% CP()
myFreqTable(dat, "rec_sess_home_other_bin") %>% CP()

medTable(dat, "rec_sess_home_in_person_no") %>% CP()
medTable(dat, "rec_sess_home_video_group_no") %>% CP()
medTable(dat, "rec_sess_home_video_indiv_no") %>% CP()
medTable(dat, "rec_sess_home_phone_no") %>% CP()
medTable(dat, "rec_sess_home_other_no") %>% CP()
medTable(dat, "rec_sess_home_no") %>% CP()
myFreqTable(dat, "discharge_assess_bin") %>% CP()

completion_ratio <- data.frame(completion_ratio = round(nrow(filter(dat, discharge_assess_bin == "Yes"))/
                                                          nrow(filter(dat, discharge_assess_bin == "No")), 1))

completion_ratio


myFreqTable(dat, "discharge_assess_no_reason") %>% CP()
myFreqTable(dat, "PR_type_location") %>% CP()
myFreqTable(dat, "discharge_assess_bin_by_rolling") %>% CP()
myFreqTable(dat, "discharge_assess_bin_by_cohort") %>% CP()
myFreqTable(dat, "discharge_assess_bin_by_home") %>% CP()
myFreqTable(dat, "discharge_assess_bin_by_both") %>% CP()





myFreqTable(dat, "exercise_plan_by_rolling") %>% CP()
myFreqTable(dat, "exercise_plan_by_cohort") %>% CP()
myFreqTable(dat, "exercise_plan_by_home") %>% CP()
myFreqTable(dat, "exercise_plan_by_both") %>% CP()
myFreqTable(dat, "exercise_plan") %>% CP()
medTable(dat, "assess_to_discharge_days") %>% CP()
myFreqTable(dat, "MRC_score_dis") %>% CP()
myFreqTable(dat, "MRC_change_factor") %>% CP()



myFreqTable(dat, "all_3_test_types_dis") %>% CP()
myFreqTable(dat, "who_is_ESWT_given_to_dis") %>% CP()
medTable(dat, "test_value_ISWT_diff") %>% CP() # Not normal
medTable(dat, "test_value_6MWT_diff") %>% CP() # Not normal
medTable(dat, "test_value_ESWT_diff") %>% CP() # Not normal
myFreqTable(dat, "MCID_N_ISWT_6MWT") %>% CP()
myFreqTable(dat, "MCID_ISWT") %>% CP()
myFreqTable(dat, "MCID_6MWT") %>% CP() 
myFreqTable(dat, "CAT_CRQ_dis_N") %>% CP()
myFreqTable(dat, "CAT_dis") %>% CP()
myFreqTable(dat, "CRQ_dis") %>% CP()
myFreqTable(dat, "CAT_CRQ_score_diff_N") %>% CP()
medTable(dat, "CAT_score_diff") %>% CP()
medTable(dat, "CRQ_dyspnoea_diff") %>% CP()
medTable(dat, "CRQ_fatigue_diff") %>% CP()
medTable(dat, "CRQ_emotion_diff") %>% CP()
medTable(dat, "CRQ_mastery_diff") %>% CP()
myFreqTable(dat, "MCID_CAT") %>% CP()
myFreqTable(dat, "MCID_CRQ_dyspnoea") %>% CP() 
myFreqTable(dat, "MCID_CRQ_fatigue") %>% CP() 
myFreqTable(dat, "MCID_CRQ_emotion") %>% CP() 
myFreqTable(dat, "MCID_CRQ_mastery") %>% CP()
myFreqTable(dat, "MCID_N_CAT_CRQ") %>% CP()
myFreqTable(dat, "MCID_N_CAT_CRQ") %>% CP()
myFreqTable(dat, "BM_start_90") %>% CP()
myFreqTable(dat, "BM_prac_test") %>% CP()
myFreqTable(dat, "BM_discharge_assess") %>% CP()
myFreqTable(dat, "BM_exercise_plan") %>% CP()
myFreqTable(dat, "BM_MCID_exercise") %>% CP()
myFreqTable(dat, "BM_MCID_CAT_CRQ") %>% CP()



dat %>% select(MCID_CAT, MCID_CRQ_dyspnoea) %>% table(useNA = "ifany")
dat %>% filter(!is.na(MCID_CAT) & !is.na(MCID_CRQ_dyspnoea)) %>% nrow()



dat %>% select(test_type_init, test_type_dis) %>% table(useNA = "ifany")
dat %>% select(test_type_init, discharge_assess_bin) %>% table(useNA = "ifany")
dat %>% filter(discharge_assess_bin == "Yes") %>% filter(test_type_init == "None") %>%
  select(test_type_dis) %>% table()

dat %>% select(all_3_test_types_init, all_3_test_types_dis) %>% table()

# Cleaning complete! Now let's start...

# But, lets also save the clean dataset from here:
# Actually though we should wait until it's complete, after looking at everything Sally wants.


data.frame(x1 = colnames(dat)) %>% CP()


# . . .   t h e   a n a l y s i s   . . . . . . .


psychic <- psych::describe(dat, skew = FALSE, ranges = FALSE, quant = c(0.25, 0.5, 0.75))
psychic <- as.data.frame(psychic)
psychic$vars <- row.names(psychic)
psychic <- psychic %>% rename(N = n, median = Q0.5, lo.quart = Q0.25, hi.quart = Q0.75)



completion_ratio <- data.frame(completion_ratio = round(nrow(filter(dat, discharge_assess_bin == "Yes"))/
                    nrow(filter(dat, discharge_assess_bin == "No")), 1))


# From here, we start off with all, and then re-run it for the individual countries.



flat <- data.frame(country = "All")

flat <- cbind(flat, mediSumRound(dat, "age", 0), # Normal-ish
              FreqSum(dat, "gender"),
              FreqSum(dat, "IMD_quintile_Eng"),
              FreqSum(dat, "IMD_quintile_Wal"),
              FreqSum(dat, "IMD_quintile_Scot"),
              FreqSum(dat, "anyIMD"),
            #  FreqSum(dat, "ethnicity"), # removed
              FreqSum(dat, "ref_location"),
#             mediSumRound(dat, "ref_to_start_days", 0), # Not normal. Changed to ref2start with stable COPD
              mediSumRound(dat, "ref_to_start_days_stable_COPD", 0),
              FreqSum(dat, "ninety_day_referral_to_start_for_stable_COPD"),
              FreqSum(dat, "thirty_day_referral_to_start_for_AECOPD"),
              mediSumRound(dat, "ref_to_start_days_cohort", 0), # For the bar chart
              mediSumRound(dat, "ref_to_start_days_rolling", 0), # For the bar chart
              mediSumRound(dat, "assess_to_start_days_stable_COPD", 0), # Not normal
              FreqSum(dat, "smoke_status"),
              FreqSum(dat, "FEV1_percpred_rec"),
              mediSumRound(dat, "FEV1_percpred", 0), # Not normal
              FreqSum(dat, "FEV1FVC_rec"),
              mediSumRound(dat, "FEV1FVC", 2), # Not normal
              FreqSum(dat, "BMI_rec"),
              mediSumRound(dat, "BMI", 1), # Not normal
              FreqSum(dat, "MRC_score_init"),
              FreqSum(dat, "CVD_history"),
              FreqSum(dat, "musc_skel_history"),
              FreqSum(dat, "mental_history"),
              FreqSum(dat, "anxiety_bin"),
              FreqSum(dat, "depression_bin"),
              FreqSum(dat, "SMI_bin"),
              FreqSum(dat, "all_3_test_types_init"),
              FreqSum(dat, "who_only_gets_ESWT_init"),
              FreqSum(dat, "which_other_test_with_ESWT_init"),
              mediSumRound(dat, "test_value_ISWT_init", 0), # Not normal
              FreqSum(dat, "prac_test_ISWT_init"),
              mediSumRound(dat, "test_value_6MWT_init", 0), # Not normal
              FreqSum(dat, "prac_test_6MWT_init"),
              mediSumRound(dat, "ESWT_value_init", 0), # Not normal
              FreqSum(dat, "CAT_init"),
              FreqSum(dat, "CRQ_init"),
              mediSumRound(dat, "CAT_score_init", 0), # Reasonably normal
              mediSumRound(dat, "CRQ_dyspnoea_init", 1), # Not normal
              mediSumRound(dat, "CRQ_fatigue_init", 1), # 
              mediSumRound(dat, "CRQ_emotion_init", 1), # 
              mediSumRound(dat, "CRQ_mastery_init", 1),
              FreqSum(dat, "enrolled"),   # seems like a lot of issues are getting patients to turn up 
                                            # to assessment maybe?
              FreqSum(dat, "PR_location"),
              FreqSum(dat, "prog_type"),
              mediSumRound(dat, "scheduled_sess_centre_no", 0),
              FreqSum(dat, "rec_sess_centre_group_bin"),
              FreqSum(dat, "rec_sess_centre_indiv_bin"),
              mediSumRound(dat, "rec_sess_centre_group_no", 0),
              mediSumRound(dat, "rec_sess_centre_indiv_no", 0),
              mediSumRound(dat, "rec_sess_centre_no", 0),
              mediSumRound(dat, "scheduled_sess_home_no", 0),
              FreqSum(dat, "rec_sess_home_in_person_bin"),
              FreqSum(dat, "rec_sess_home_video_group_bin"),
              FreqSum(dat, "rec_sess_home_video_indiv_bin"),
              FreqSum(dat, "rec_sess_home_phone_bin"),
              FreqSum(dat, "rec_sess_home_other_bin"),
              mediSumRound(dat, "rec_sess_home_no", 0),
              mediSumRound(dat, "rec_sess_home_in_person_no", 0),
              mediSumRound(dat, "rec_sess_home_video_group_no", 0),
              mediSumRound(dat, "rec_sess_home_video_indiv_no", 0),
              mediSumRound(dat, "rec_sess_home_phone_no", 0),
              mediSumRound(dat, "rec_sess_home_other_no", 0),
              FreqSum(dat, "discharge_assess_bin"),
              completion_ratio,
              FreqSum(dat, "discharge_assess_no_reason"),
              FreqSum(dat, "PR_type_location"),
              FreqSum(dat, "discharge_assess_bin_by_rolling"),
              FreqSum(dat, "discharge_assess_bin_by_cohort"),
              FreqSum(dat, "discharge_assess_bin_by_home"),
              FreqSum(dat, "discharge_assess_bin_by_both"),
              FreqSum(dat, "exercise_plan"),
              FreqSum(dat, "exercise_plan_by_rolling"),
              FreqSum(dat, "exercise_plan_by_cohort"),
              FreqSum(dat, "exercise_plan_by_home"),
              FreqSum(dat, "exercise_plan_by_both"),
              mediSumRound(dat, "assess_to_discharge_days", 0),
              FreqSum(dat, "MRC_score_dis"),
              FreqSum(dat, "MRC_change_factor"),


data.frame(MRC_init_1_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,1]),
data.frame(MRC_init_1_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,2]),
data.frame(MRC_init_1_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,3]),
data.frame(MRC_init_1_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,4]),
data.frame(MRC_init_1_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,5]),
data.frame(MRC_init_1_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,6]),

data.frame(MRC_init_2_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,1]),
data.frame(MRC_init_2_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,2]),
data.frame(MRC_init_2_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,3]),
data.frame(MRC_init_2_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,4]),
data.frame(MRC_init_2_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,5]),
data.frame(MRC_init_2_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,6]),

data.frame(MRC_init_3_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,1]),
data.frame(MRC_init_3_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,2]),
data.frame(MRC_init_3_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,3]),
data.frame(MRC_init_3_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,4]),
data.frame(MRC_init_3_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,5]),
data.frame(MRC_init_3_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,6]),

data.frame(MRC_init_4_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,1]),
data.frame(MRC_init_4_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,2]),
data.frame(MRC_init_4_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,3]),
data.frame(MRC_init_4_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,4]),
data.frame(MRC_init_4_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,5]),
data.frame(MRC_init_4_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,6]),

data.frame(MRC_init_5_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,1]),
data.frame(MRC_init_5_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,2]),
data.frame(MRC_init_5_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,3]),
data.frame(MRC_init_5_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,4]),
data.frame(MRC_init_5_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,5]),
data.frame(MRC_init_5_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,6]),

data.frame(MRC_init_NR_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,1]),
data.frame(MRC_init_NR_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,2]),
data.frame(MRC_init_NR_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,3]),
data.frame(MRC_init_NR_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,4]),
data.frame(MRC_init_NR_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,5]),
data.frame(MRC_init_NR_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,6]),

FreqSum(dat, "all_3_test_types_dis"),
# FreqSum(dat, "who_is_ESWT_given_to_dis"),
mediSumRound(dat, "test_value_ISWT_diff", 1), # Not normal
mediSumRound(dat, "test_value_6MWT_diff", 0), # Not normal
mediSumRound(dat, "test_value_ESWT_diff", 0), # Not normal
FreqSum(dat, "MCID_ISWT"),
FreqSum(dat, "MCID_6MWT"), 
FreqSum(dat, "CAT_dis"),
FreqSum(dat, "CRQ_dis"),
mediSumRound(dat, "CAT_score_diff", 0),
mediSumRound(dat, "CRQ_dyspnoea_diff", 1),
mediSumRound(dat, "CRQ_fatigue_diff", 1),
mediSumRound(dat, "CRQ_emotion_diff", 1),
mediSumRound(dat, "CRQ_mastery_diff", 1),
FreqSum(dat, "MCID_CAT"),
FreqSum(dat, "MCID_CRQ_dyspnoea"), 
FreqSum(dat, "MCID_CRQ_fatigue"), 
FreqSum(dat, "MCID_CRQ_emotion"), 
FreqSum(dat, "MCID_CRQ_mastery"))



flat.all <- flat

# To keep the row order the same for the org level, we use the national level one
# as a template to follow so we save it as flat.org as well.

flat.org <- flat

# Need to replace the country column with something else and add another for trust code:

flat.org <- flat.org[ ,c(1, 1, 1, 2:ncol(flat.org))]
flat.org <- flat.org %>% rename(org_code = country,
                                trust_code = country.1, 
                                patient_N = country.2) %>% 
                         mutate(org_code = as.character(org_code),
                                trust_code = as.character(trust_code),
                                patient_N = as.character(patient_N))
                                
flat.org$patient_N <- NA
flat.org$patient_N <- as.numeric(flat.org$patient_N)


dat.save <- dat

#### For country...

for (i in unique(dat.save$country)) {
  
  dat <- filter(dat.save, country == i)

psychic <- psych::describe(dat, skew = FALSE, ranges = FALSE, quant = c(0.25, 0.5, 0.75))
psychic <- as.data.frame(psychic)
psychic$vars <- row.names(psychic)
psychic <- psychic %>% rename(N = n, median = Q0.5, lo.quart = Q0.25, hi.quart = Q0.75)



completion_ratio <- data.frame(completion_ratio = round(nrow(filter(dat, discharge_assess_bin == "Yes"))/
                                                          nrow(filter(dat, discharge_assess_bin == "No")), 1))



# From here, we start off with all, and then re-run it for the individual countries.



flat <- data.frame(country = i)

flat <- cbind(flat, mediSumRound(dat, "age", 0), # Normal-ish
              FreqSum(dat, "gender"),
              FreqSum(dat, "IMD_quintile_Eng"),
              FreqSum(dat, "IMD_quintile_Wal"),
              FreqSum(dat, "IMD_quintile_Scot"),
              FreqSum(dat, "anyIMD"),
            #  FreqSum(dat, "ethnicity"), # removed
              FreqSum(dat, "ref_location"),
#             mediSumRound(dat, "ref_to_start_days", 0), # Not normal. Changed to ref2start with stable COPD
              mediSumRound(dat, "ref_to_start_days_stable_COPD", 0),
              FreqSum(dat, "ninety_day_referral_to_start_for_stable_COPD"),
              FreqSum(dat, "thirty_day_referral_to_start_for_AECOPD"),
              mediSumRound(dat, "ref_to_start_days_cohort", 0), # For the bar chart
              mediSumRound(dat, "ref_to_start_days_rolling", 0), # For the bar chart
              mediSumRound(dat, "assess_to_start_days_stable_COPD", 0), # Not normal
              FreqSum(dat, "smoke_status"),
              FreqSum(dat, "FEV1_percpred_rec"),
              mediSumRound(dat, "FEV1_percpred", 0), # Not normal
              FreqSum(dat, "FEV1FVC_rec"),
              mediSumRound(dat, "FEV1FVC", 2), # Not normal
              FreqSum(dat, "BMI_rec"),
              mediSumRound(dat, "BMI", 1), # Not normal
              FreqSum(dat, "MRC_score_init"),
              FreqSum(dat, "CVD_history"),
              FreqSum(dat, "musc_skel_history"),
              FreqSum(dat, "mental_history"),
              FreqSum(dat, "anxiety_bin"),
              FreqSum(dat, "depression_bin"),
              FreqSum(dat, "SMI_bin"),
              FreqSum(dat, "all_3_test_types_init"),
              FreqSum(dat, "who_only_gets_ESWT_init"),
      #        FreqSum(dat, "who_is_ESWT_given_to_init"),
              mediSumRound(dat, "test_value_ISWT_init", 0), # Not normal
              FreqSum(dat, "prac_test_ISWT_init"),
              mediSumRound(dat, "test_value_6MWT_init", 0), # Not normal
              FreqSum(dat, "prac_test_6MWT_init"),
              mediSumRound(dat, "ESWT_value_init", 0), # Not normal
              FreqSum(dat, "CAT_init"),
              FreqSum(dat, "CRQ_init"),
              mediSumRound(dat, "CAT_score_init", 0), # Reasonably normal
              mediSumRound(dat, "CRQ_dyspnoea_init", 1), # Not normal
              mediSumRound(dat, "CRQ_fatigue_init", 1), # 
              mediSumRound(dat, "CRQ_emotion_init", 1), # 
              mediSumRound(dat, "CRQ_mastery_init", 1),
              FreqSum(dat, "enrolled"),   # seems like a lot of issues are getting patients to turn up 
              # to assessment maybe?
              FreqSum(dat, "PR_location"),
              FreqSum(dat, "prog_type"),
              mediSumRound(dat, "scheduled_sess_centre_no", 0),
              FreqSum(dat, "rec_sess_centre_group_bin"),
              FreqSum(dat, "rec_sess_centre_indiv_bin"),
              mediSumRound(dat, "rec_sess_centre_group_no", 0),
              mediSumRound(dat, "rec_sess_centre_indiv_no", 0),
              mediSumRound(dat, "rec_sess_centre_no", 0),
              mediSumRound(dat, "scheduled_sess_home_no", 0),
              FreqSum(dat, "rec_sess_home_in_person_bin"),
              FreqSum(dat, "rec_sess_home_video_group_bin"),
              FreqSum(dat, "rec_sess_home_video_indiv_bin"),
              FreqSum(dat, "rec_sess_home_phone_bin"),
              FreqSum(dat, "rec_sess_home_other_bin"),
              mediSumRound(dat, "rec_sess_home_no", 0),
              mediSumRound(dat, "rec_sess_home_in_person_no", 0),
              mediSumRound(dat, "rec_sess_home_video_group_no", 0),
              mediSumRound(dat, "rec_sess_home_video_indiv_no", 0),
              mediSumRound(dat, "rec_sess_home_phone_no", 0),
              mediSumRound(dat, "rec_sess_home_other_no", 0),
              FreqSum(dat, "discharge_assess_bin"),
              completion_ratio,
              FreqSum(dat, "discharge_assess_no_reason"),
              FreqSum(dat, "PR_type_location"),
              FreqSum(dat, "discharge_assess_bin_by_rolling"),
              FreqSum(dat, "discharge_assess_bin_by_cohort"),
              FreqSum(dat, "discharge_assess_bin_by_home"),
              FreqSum(dat, "discharge_assess_bin_by_both"),
              FreqSum(dat, "exercise_plan"),
              FreqSum(dat, "exercise_plan_by_rolling"),
              FreqSum(dat, "exercise_plan_by_cohort"),
              FreqSum(dat, "exercise_plan_by_home"),
              FreqSum(dat, "exercise_plan_by_both"),
              mediSumRound(dat, "assess_to_discharge_days", 0),
              FreqSum(dat, "MRC_score_dis"),
              FreqSum(dat, "MRC_change_factor"),
              
              
              data.frame(MRC_init_1_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,1]),
              data.frame(MRC_init_1_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,2]),
              data.frame(MRC_init_1_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,3]),
              data.frame(MRC_init_1_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,4]),
              data.frame(MRC_init_1_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,5]),
              data.frame(MRC_init_1_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,6]),
              
              data.frame(MRC_init_2_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,1]),
              data.frame(MRC_init_2_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,2]),
              data.frame(MRC_init_2_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,3]),
              data.frame(MRC_init_2_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,4]),
              data.frame(MRC_init_2_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,5]),
              data.frame(MRC_init_2_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,6]),
              
              data.frame(MRC_init_3_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,1]),
              data.frame(MRC_init_3_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,2]),
              data.frame(MRC_init_3_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,3]),
              data.frame(MRC_init_3_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,4]),
              data.frame(MRC_init_3_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,5]),
              data.frame(MRC_init_3_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,6]),
              
              data.frame(MRC_init_4_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,1]),
              data.frame(MRC_init_4_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,2]),
              data.frame(MRC_init_4_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,3]),
              data.frame(MRC_init_4_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,4]),
              data.frame(MRC_init_4_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,5]),
              data.frame(MRC_init_4_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,6]),
              
              data.frame(MRC_init_5_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,1]),
              data.frame(MRC_init_5_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,2]),
              data.frame(MRC_init_5_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,3]),
              data.frame(MRC_init_5_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,4]),
              data.frame(MRC_init_5_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,5]),
              data.frame(MRC_init_5_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,6]),
              
              data.frame(MRC_init_NR_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,1]),
              data.frame(MRC_init_NR_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,2]),
              data.frame(MRC_init_NR_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,3]),
              data.frame(MRC_init_NR_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,4]),
              data.frame(MRC_init_NR_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,5]),
              data.frame(MRC_init_NR_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,6]),
              
              FreqSum(dat, "all_3_test_types_dis"),
       #       FreqSum(dat, "who_is_ESWT_given_to_dis"),
              mediSumRound(dat, "test_value_ISWT_diff", 1), # Not normal
              mediSumRound(dat, "test_value_6MWT_diff", 0), # Not normal
              mediSumRound(dat, "test_value_ESWT_diff", 0), # Not normal
              FreqSum(dat, "MCID_ISWT"),
              FreqSum(dat, "MCID_6MWT"), 
              FreqSum(dat, "CAT_dis"),
              FreqSum(dat, "CRQ_dis"),
              mediSumRound(dat, "CAT_score_diff", 0),
              mediSumRound(dat, "CRQ_dyspnoea_diff", 1),
              mediSumRound(dat, "CRQ_fatigue_diff", 1),
              mediSumRound(dat, "CRQ_emotion_diff", 1),
              mediSumRound(dat, "CRQ_mastery_diff", 1),
              FreqSum(dat, "MCID_CAT"),
              FreqSum(dat, "MCID_CRQ_dyspnoea"), 
              FreqSum(dat, "MCID_CRQ_fatigue"), 
              FreqSum(dat, "MCID_CRQ_emotion"), 
              FreqSum(dat, "MCID_CRQ_mastery"))

flat.all <- bind_rows(flat.all, flat)

}


# National level.

glimpse(flat.all)
# write.csv(flat.all,
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/PR_clinical_national_report_data_2020-01-14.csv",
#           row.names = FALSE)


# write.csv(flat.all,
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/PR_clinical_national_report_data_2020-04-01.csv",
#           row.names = FALSE)


# write.csv(flat.all,
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/PR_clinical_national_report_data_2020-04-03.csv",
#           row.names = FALSE)


# And then we get our old dat back:

dat <- dat.save

# Hospital level

#### For organisational...

for (i in unique(dat.save$org_code)) {
  
  dat <- filter(dat.save, org_code == i)
  
  psychic <- psych::describe(dat, skew = FALSE, ranges = FALSE, quant = c(0.25, 0.5, 0.75))
  psychic <- as.data.frame(psychic)
  psychic$vars <- row.names(psychic)
  psychic <- psychic %>% rename(N = n, median = Q0.5, lo.quart = Q0.25, hi.quart = Q0.75)
  
  
  
  completion_ratio <- data.frame(completion_ratio = round(nrow(filter(dat, discharge_assess_bin == "Yes"))/
                                                            nrow(filter(dat, discharge_assess_bin == "No")), 1))
  
 
  
  # From here, we start off with all, and then re-run it for the individual countries.
  
  
  
  flat <- data.frame(org_code = i)
  flat$trust_code <- as.character(dat$trust_code[1])
  flat$patient_N <- nrow(dat)
  
  flat <- cbind(flat, mediSumRound(dat, "age", 0), # Normal-ish
                FreqSum(dat, "gender"),
                FreqSum(dat, "IMD_quintile_Eng"),
                FreqSum(dat, "IMD_quintile_Wal"),
                FreqSum(dat, "IMD_quintile_Scot"),
                FreqSum(dat, "anyIMD"),
                #  FreqSum(dat, "ethnicity"), # removed
                FreqSum(dat, "ref_location"),
                #             mediSumRound(dat, "ref_to_start_days", 0), # Not normal. Changed to ref2start with stable COPD
                mediSumRound(dat, "ref_to_start_days_stable_COPD", 0),
                FreqSum(dat, "ninety_day_referral_to_start_for_stable_COPD"),
                FreqSum(dat, "thirty_day_referral_to_start_for_AECOPD"),
                mediSumRound(dat, "ref_to_start_days_cohort", 0), # For the bar chart
                mediSumRound(dat, "ref_to_start_days_rolling", 0), # For the bar chart
                mediSumRound(dat, "assess_to_start_days_stable_COPD", 0), # Not normal
                FreqSum(dat, "smoke_status"),
                FreqSum(dat, "FEV1_percpred_rec"),
                mediSumRound(dat, "FEV1_percpred", 0), # Not normal
                FreqSum(dat, "FEV1FVC_rec"),
                mediSumRound(dat, "FEV1FVC", 2), # Not normal
                FreqSum(dat, "BMI_rec"),
                mediSumRound(dat, "BMI", 1), # Not normal
                FreqSum(dat, "MRC_score_init"),
                FreqSum(dat, "CVD_history"),
                FreqSum(dat, "musc_skel_history"),
                FreqSum(dat, "mental_history"),
                FreqSum(dat, "anxiety_bin"),
                FreqSum(dat, "depression_bin"),
                FreqSum(dat, "SMI_bin"),
                FreqSum(dat, "all_3_test_types_init"),
                FreqSum(dat, "who_only_gets_ESWT_init"),
                #        FreqSum(dat, "who_is_ESWT_given_to_init"),
                mediSumRound(dat, "test_value_ISWT_init", 0), # Not normal
                FreqSum(dat, "prac_test_ISWT_init"),
                mediSumRound(dat, "test_value_6MWT_init", 0), # Not normal
                FreqSum(dat, "prac_test_6MWT_init"),
                mediSumRound(dat, "ESWT_value_init", 0), # Not normal
                FreqSum(dat, "CAT_init"),
                FreqSum(dat, "CRQ_init"),
                mediSumRound(dat, "CAT_score_init", 0), # Reasonably normal
                mediSumRound(dat, "CRQ_dyspnoea_init", 1), # Not normal
                mediSumRound(dat, "CRQ_fatigue_init", 1), # 
                mediSumRound(dat, "CRQ_emotion_init", 1), # 
                mediSumRound(dat, "CRQ_mastery_init", 1),
                FreqSum(dat, "enrolled"),   # seems like a lot of issues are getting patients to turn up 
                # to assessment maybe?
                FreqSum(dat, "PR_location"),
                FreqSum(dat, "prog_type"),
                mediSumRound(dat, "scheduled_sess_centre_no", 0),
                FreqSum(dat, "rec_sess_centre_group_bin"),
                FreqSum(dat, "rec_sess_centre_indiv_bin"),
                mediSumRound(dat, "rec_sess_centre_group_no", 0),
                mediSumRound(dat, "rec_sess_centre_indiv_no", 0),
                mediSumRound(dat, "rec_sess_centre_no", 0),
                mediSumRound(dat, "scheduled_sess_home_no", 0),
                FreqSum(dat, "rec_sess_home_in_person_bin"),
                FreqSum(dat, "rec_sess_home_video_group_bin"),
                FreqSum(dat, "rec_sess_home_video_indiv_bin"),
                FreqSum(dat, "rec_sess_home_phone_bin"),
                FreqSum(dat, "rec_sess_home_other_bin"),
                mediSumRound(dat, "rec_sess_home_no", 0),
                mediSumRound(dat, "rec_sess_home_in_person_no", 0),
                mediSumRound(dat, "rec_sess_home_video_group_no", 0),
                mediSumRound(dat, "rec_sess_home_video_indiv_no", 0),
                mediSumRound(dat, "rec_sess_home_phone_no", 0),
                mediSumRound(dat, "rec_sess_home_other_no", 0),
                FreqSum(dat, "discharge_assess_bin"),
                completion_ratio,
                FreqSum(dat, "discharge_assess_no_reason"),
                FreqSum(dat, "PR_type_location"),
                FreqSum(dat, "discharge_assess_bin_by_rolling"),
                FreqSum(dat, "discharge_assess_bin_by_cohort"),
                FreqSum(dat, "discharge_assess_bin_by_home"),
                FreqSum(dat, "discharge_assess_bin_by_both"),
                FreqSum(dat, "exercise_plan"),
                FreqSum(dat, "exercise_plan_by_rolling"),
                FreqSum(dat, "exercise_plan_by_cohort"),
                FreqSum(dat, "exercise_plan_by_home"),
                FreqSum(dat, "exercise_plan_by_both"),
                mediSumRound(dat, "assess_to_discharge_days", 0),
                FreqSum(dat, "MRC_score_dis"),
                FreqSum(dat, "MRC_change_factor"),
                
                
                data.frame(MRC_init_1_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,1]),
                data.frame(MRC_init_1_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,2]),
                data.frame(MRC_init_1_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,3]),
                data.frame(MRC_init_1_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,4]),
                data.frame(MRC_init_1_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,5]),
                data.frame(MRC_init_1_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[1,6]),
                
                data.frame(MRC_init_2_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,1]),
                data.frame(MRC_init_2_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,2]),
                data.frame(MRC_init_2_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,3]),
                data.frame(MRC_init_2_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,4]),
                data.frame(MRC_init_2_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,5]),
                data.frame(MRC_init_2_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[2,6]),
                
                data.frame(MRC_init_3_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,1]),
                data.frame(MRC_init_3_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,2]),
                data.frame(MRC_init_3_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,3]),
                data.frame(MRC_init_3_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,4]),
                data.frame(MRC_init_3_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,5]),
                data.frame(MRC_init_3_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[3,6]),
                
                data.frame(MRC_init_4_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,1]),
                data.frame(MRC_init_4_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,2]),
                data.frame(MRC_init_4_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,3]),
                data.frame(MRC_init_4_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,4]),
                data.frame(MRC_init_4_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,5]),
                data.frame(MRC_init_4_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[4,6]),
                
                data.frame(MRC_init_5_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,1]),
                data.frame(MRC_init_5_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,2]),
                data.frame(MRC_init_5_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,3]),
                data.frame(MRC_init_5_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,4]),
                data.frame(MRC_init_5_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,5]),
                data.frame(MRC_init_5_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[5,6]),
                
                data.frame(MRC_init_NR_dis_1_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,1]),
                data.frame(MRC_init_NR_dis_2_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,2]),
                data.frame(MRC_init_NR_dis_3_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,3]),
                data.frame(MRC_init_NR_dis_4_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,4]),
                data.frame(MRC_init_NR_dis_5_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,5]),
                data.frame(MRC_init_NR_dis_NR_N = table(dat$MRC_score_init, dat$MRC_score_dis)[6,6]),
                
                FreqSum(dat, "all_3_test_types_dis"),
                #       FreqSum(dat, "who_is_ESWT_given_to_dis"),
                mediSumRound(dat, "test_value_ISWT_diff", 1), # Not normal
                mediSumRound(dat, "test_value_6MWT_diff", 0), # Not normal
                mediSumRound(dat, "test_value_ESWT_diff", 0), # Not normal
                FreqSum(dat, "MCID_ISWT"),
                FreqSum(dat, "MCID_6MWT"), 
                FreqSum(dat, "CAT_dis"),
                FreqSum(dat, "CRQ_dis"),
                mediSumRound(dat, "CAT_score_diff", 0),
                mediSumRound(dat, "CRQ_dyspnoea_diff", 1),
                mediSumRound(dat, "CRQ_fatigue_diff", 1),
                mediSumRound(dat, "CRQ_emotion_diff", 1),
                mediSumRound(dat, "CRQ_mastery_diff", 1),
                FreqSum(dat, "MCID_CAT"),
                FreqSum(dat, "MCID_CRQ_dyspnoea"), 
                FreqSum(dat, "MCID_CRQ_fatigue"), 
                FreqSum(dat, "MCID_CRQ_emotion"), 
                FreqSum(dat, "MCID_CRQ_mastery"))
  
  flat.org <- bind_rows(flat.org, flat)
  
}

# Remove these confusing unnecessary columns:

flat.org <- flat.org %>% select(
  -(who_only_gets_ESWT_init_N:which_other_test_with_ESWT_init_ISWT_perc))

# and remove the first line that is used as a template:

flat.org <- flat.org[-1, ]
glimpse(flat.org)

# write.csv(flat.org,
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/PR_clinical_org_level_data_2020-04-17.csv",
#           row.names = FALSE)


# let's recreate 'dat'.

dat <- dat.save


dat %>% select(ref_location) %>% table()




# Use summarise function to get necessary columns
bmk <- dat %>% dplyr::group_by(org_code) %>%
      summarise(trust_code = first(trust_code), 
                cases.audited = n(),
                start_90_denom = sum(!is.na(BM_start_90)),
                start_90_nume = sum(BM_start_90, na.rm = TRUE),
                start_90_perc = (start_90_nume/start_90_denom)*100,

BM_prac_test_denom = sum(!is.na(BM_prac_test)),
BM_prac_test_nume = sum(BM_prac_test, na.rm = TRUE),
BM_prac_test_perc = (BM_prac_test_nume/BM_prac_test_denom)*100,

BM_discharge_assess_denom = sum(!is.na(BM_discharge_assess)),
BM_discharge_assess_nume = sum(BM_discharge_assess, na.rm = TRUE),
BM_discharge_assess_perc = (BM_discharge_assess_nume/BM_discharge_assess_denom)*100,

BM_exercise_plan_denom = sum(!is.na(BM_exercise_plan)),
BM_exercise_plan_nume = sum(BM_exercise_plan, na.rm = TRUE),
BM_exercise_plan_perc = (BM_exercise_plan_nume/BM_exercise_plan_denom)*100,

BM_MCID_exercise_denom = sum(!is.na(BM_MCID_exercise)),
BM_MCID_exercise_nume = sum(BM_MCID_exercise, na.rm = TRUE),
BM_MCID_exercise_perc = (BM_MCID_exercise_nume/BM_MCID_exercise_denom)*100,

BM_MCID_CAT_CRQ_denom = sum(!is.na(BM_MCID_CAT_CRQ)),
BM_MCID_CAT_CRQ_nume = sum(BM_MCID_CAT_CRQ, na.rm = TRUE),
BM_MCID_CAT_CRQ_perc = (BM_MCID_CAT_CRQ_nume/BM_MCID_CAT_CRQ_denom)*100)

bmk
# quartz1 is for calculating stuff, quartz_fmt is the well-formatted one



quartz1 <- matrix(data = NA, nrow = 3, ncol = 7)
quartz1[1:3, 1] <- c("lower.quartile", "median", "upper.quartile")

quartz1[1:3, 2] <- quantile(bmk$start_90_perc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4)
quartz1[1:3, 3] <- quantile(bmk$BM_prac_test_perc,probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4)
quartz1[1:3, 4] <- quantile(bmk$BM_discharge_assess_perc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4)
quartz1[1:3, 5] <- quantile(bmk$BM_exercise_plan_perc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4)
quartz1[1:3, 6] <- quantile(bmk$BM_MCID_exercise_perc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4)
quartz1[1:3, 7] <- quantile(bmk$BM_MCID_CAT_CRQ_perc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4)


colnames(quartz1) <- c("statistic", "start_90_perc", "BM_prac_test_perc", "BM_discharge_assess_perc", 
                      "BM_exercise_plan_perc", "BM_MCID_exercise_perc", "BM_MCID_CAT_CRQ_perc") 

quartz1 <- as.data.frame(quartz1)
# quartz1 %>% mutate_if(is.factor, as.character(.)) %>% mutate_at(~vars(-statistic), ~as.numeric)

quartz1 <- quartz1 %>% mutate_at(.vars = vars(-statistic), .funs = ~as.numeric(as.character(.)))

quartz1 <- quartz1 %>% mutate_at(.vars = vars(-statistic), .funs = ~round(., 0))

# Now that we're rounding the medians anyway, this is a very long-winded way to do it and I could have 
# just used quartz1 to make quartz_fmt

quartz_fmt <- matrix(data = NA, nrow = 3, ncol = 7)
quartz_fmt[1:3, 1] <- c("lower.quartile", "median", "upper.quartile")

quartz_fmt[1:3, 2] <- sprintf("%.0f", round(quantile(bmk$start_90_perc,
                                                 probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4),0))
quartz_fmt[1:3, 3] <- sprintf("%.0f", round(quantile(bmk$BM_prac_test_perc,
                                                 probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4),0))
quartz_fmt[1:3, 4] <- sprintf("%.0f", round(quantile(bmk$BM_discharge_assess_perc,
                                                 probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4),0))
quartz_fmt[1:3, 5] <- sprintf("%.0f", round(quantile(bmk$BM_exercise_plan_perc,
                                                 probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4),0))
quartz_fmt[1:3, 6] <- sprintf("%.0f", round(quantile(bmk$BM_MCID_exercise_perc,
                                                 probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4),0))
quartz_fmt[1:3, 7] <- sprintf("%.0f", round(quantile(bmk$BM_MCID_CAT_CRQ_perc,
                                                 probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 4),0))


colnames(quartz_fmt) <- c("statistic", "start_90_perc", "BM_prac_test_perc", "BM_discharge_assess_perc", 
                      "BM_exercise_plan_perc", "BM_MCID_exercise_perc", "BM_MCID_CAT_CRQ_perc") 

quartz_fmt <- as.data.frame(quartz_fmt)

# write.csv(quartz_fmt, file =
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/benchmarking_quartiles_2020-04-07.csv",
#           row.names = FALSE)


# It's at this point that we round BMK so it can be compared to the medians. 

colnames(bmk)

bmk <- bmk %>% mutate_at(.vars = vars(contains("perc")), .funs = ~round(., 0))



# Now, using quartz1, we add in the BMK colour code.

bmk <- bmk %>% mutate(start_90_colour_end = ifelse(start_90_denom < 5 | is.na(start_90_denom) == TRUE, "Grey",
                                               ifelse(start_90_perc < quartz1$start_90_perc[1], "Red",
                                                      ifelse(start_90_perc >= quartz1$start_90_perc[3], "Green", 
                                                             "Yellow"))),
                      BM_prac_test_colour_end = ifelse(BM_prac_test_denom < 5 | is.na(BM_prac_test_denom) == TRUE, "Grey",
                                               ifelse(BM_prac_test_perc < quartz1$BM_prac_test_perc[1], "Red",
                                                      ifelse(BM_prac_test_perc >= quartz1$BM_prac_test_perc[3], "Green", 
                                                             "Yellow"))),
                      BM_discharge_assess_colour_end = ifelse(BM_discharge_assess_denom < 5 | is.na(BM_discharge_assess_denom) == TRUE, "Grey",
                                               ifelse(BM_discharge_assess_perc < quartz1$BM_discharge_assess_perc[1], "Red",
                                                      ifelse(BM_discharge_assess_perc >= quartz1$BM_discharge_assess_perc[3], "Green", 
                                                             "Yellow"))),
                      BM_exercise_plan_colour_end = ifelse(BM_exercise_plan_denom < 5 | is.na(BM_exercise_plan_denom) == TRUE, "Grey",
                                               ifelse(BM_exercise_plan_perc < quartz1$BM_exercise_plan_perc[1], "Red",
                                                      ifelse(BM_exercise_plan_perc >= quartz1$BM_exercise_plan_perc[3], "Green", 
                                                             "Yellow"))),
                      BM_MCID_exercise_colour_end = ifelse(BM_MCID_exercise_denom < 5 | is.na(BM_MCID_exercise_denom) == TRUE, "Grey",
                                               ifelse(BM_MCID_exercise_perc < quartz1$BM_MCID_exercise_perc[1], "Red",
                                                      ifelse(BM_MCID_exercise_perc >= quartz1$BM_MCID_exercise_perc[3], "Green", 
                                                             "Yellow"))),
                      BM_MCID_CAT_CRQ_colour_end = ifelse(BM_MCID_CAT_CRQ_denom < 5 | is.na(BM_MCID_CAT_CRQ_denom) == TRUE, "Grey",
                                               ifelse(BM_MCID_CAT_CRQ_perc < quartz1$BM_MCID_CAT_CRQ_perc[1], "Red",
                                                      ifelse(BM_MCID_CAT_CRQ_perc >= quartz1$BM_MCID_CAT_CRQ_perc[3], "Green", 
                                                             "Yellow"))))






bmk <- bmk %>% add_column(start_90_colour = bmk$start_90_colour_end, .after = "start_90_perc") %>% 
  add_column(BM_prac_test_colour = bmk$BM_prac_test_colour_end, .after = "BM_prac_test_perc") %>%
  add_column(BM_discharge_assess_colour = bmk$BM_discharge_assess_colour_end, .after = "BM_discharge_assess_perc") %>% 
  add_column(BM_exercise_plan_colour = bmk$BM_exercise_plan_colour_end, .after = "BM_exercise_plan_perc") %>% 
  add_column(BM_MCID_exercise_colour = bmk$BM_MCID_exercise_colour_end, .after = "BM_MCID_exercise_perc") %>%
  add_column(BM_MCID_CAT_CRQ_colour = bmk$BM_MCID_CAT_CRQ_colour_end, .after = "BM_MCID_CAT_CRQ_perc") %>%
  select(-start_90_colour_end, -BM_prac_test_colour_end, -BM_discharge_assess_colour_end, -BM_exercise_plan_colour_end,
         -BM_MCID_exercise_colour_end, -BM_MCID_CAT_CRQ_colour_end)





bmk_all <- dat %>%
  summarise(org_code = "National",
    trust_code = "National", 
            cases.audited = n(),
            start_90_denom = sum(!is.na(BM_start_90)),
            start_90_nume = sum(BM_start_90, na.rm = TRUE),
            start_90_perc = round((start_90_nume/start_90_denom)*100, 0),
            
            BM_prac_test_denom = sum(!is.na(BM_prac_test)),
            BM_prac_test_nume = sum(BM_prac_test, na.rm = TRUE),
            BM_prac_test_perc = round((BM_prac_test_nume/BM_prac_test_denom)*100, 0),
            
            BM_discharge_assess_denom = sum(!is.na(BM_discharge_assess)),
            BM_discharge_assess_nume = sum(BM_discharge_assess, na.rm = TRUE),
            BM_discharge_assess_perc = round((BM_discharge_assess_nume/BM_discharge_assess_denom)*100, 0),
            
            BM_exercise_plan_denom = sum(!is.na(BM_exercise_plan)),
            BM_exercise_plan_nume = sum(BM_exercise_plan, na.rm = TRUE),
            BM_exercise_plan_perc = round((BM_exercise_plan_nume/BM_exercise_plan_denom)*100, 0),
            
            BM_MCID_exercise_denom = sum(!is.na(BM_MCID_exercise)),
            BM_MCID_exercise_nume = sum(BM_MCID_exercise, na.rm = TRUE),
            BM_MCID_exercise_perc = round((BM_MCID_exercise_nume/BM_MCID_exercise_denom)*100, 0),
            
            BM_MCID_CAT_CRQ_denom = sum(!is.na(BM_MCID_CAT_CRQ)),
            BM_MCID_CAT_CRQ_nume = sum(BM_MCID_CAT_CRQ, na.rm = TRUE),
            BM_MCID_CAT_CRQ_perc = round((BM_MCID_CAT_CRQ_nume/BM_MCID_CAT_CRQ_denom)*100, 0))

# We want to keep the column order of the site-level table
# We then need to change the row order so that the national analysis is at the top.
# We therefore put the last row at the top using the indexing below

bmk <- bind_rows(bmk, bmk_all)
bmk <- bmk[c(nrow(bmk), 1:(nrow(bmk)-1)), ]

bmk <- bmk %>% mutate_at(.vars = vars(matches("perc")), .funs = ~sprintf("%.0f", round(., 0)))



bmk

str(bmk)

# write.csv(bmk, file =
# "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/benchmarking_2020-04-07.csv",
#   row.names = FALSE)


# That's it for everything apart from the analyses!




#---------------------------------------------#
# A n a l y s i s   s e c t i o n
#---------------------------------------------#

# Analyses we're running:
# Likelihood of discharge assessment
# Likelihood of making MCID walking test
# Likelihood of making health status MCID

# Variables in model:
# Gender
# IMD quintile
# Age
# Physical comorbidity
# mental comorbidity
# CAT score


# Analysis 1: Likelihood of discharge assessment

# It says it's not converging, but that's just due to the small number of alternate genders, everything else
# is fine. I compared it to when those genders were dropped and the results were exactly the same.

summary(dat$gender)

# Normal one - fine but does come with warnings.

m1.1 <- glmer(discharge_assess_bin ~ 1 + gender + IMD_quintile_all + agecat + CVD_history +
                musc_skel_history + SMI_bin_no_missing + anxiety_bin_no_missing + depression_bin_no_missing + 
                CAT_score_init_cat + (1|org_code),
              data = filter(dat, gender == "Male" | gender == "Female") , family = "binomial")


# You can use this optimizer with no convergence issues.

m1.3 <- glmer(discharge_assess_bin ~ 1 + gender + IMD_quintile_all + agecat + CVD_history +
          musc_skel_history + SMI_bin_no_missing + anxiety_bin_no_missing + depression_bin_no_missing + 
            CAT_score_init_cat + (1|org_code), control=glmerControl(optimizer="bobyqa",
                                                                    optCtrl=list(maxfun=2e5)),
          data = filter(dat, gender == "Male" | gender == "Female") , family = "binomial")

# Filters out the non-males/females etc.



# This is just comparing it with if you don't take PR course into account.

m1.2 <- glm(discharge_assess_bin ~ 1 + gender + IMD_quintile_all + agecat + CVD_history +
                musc_skel_history + SMI_bin_no_missing + anxiety_bin_no_missing + depression_bin_no_missing + 
                CAT_score_init_cat,
              data = filter(dat, gender == "Male" | gender == "Female") , family = "binomial")


m1.4 <- glmer(discharge_assess_bin ~ 1 + gender + IMD_quintile_all + agecat + CVD_history +
                musc_skel_history + SMI_bin_no_missing +  
                CAT_score_init_cat + (1|org_code), control=glmerControl(optimizer="bobyqa",
                                                                        optCtrl=list(maxfun=2e5)),
              data = filter(dat, gender == "Male" | gender == "Female") , family = "binomial")


summary(m1.4)
test1.4 <- tidyoutput(m1.4, meth = "Wald", MEM = TRUE)
test1.4


summary(m1.1)
test1.1 <- tidyoutput(m1.1, meth = "Wald", MEM = TRUE)
test1.1

summary(m1.3)
test1.3 <- tidyoutput(m1.3, meth = "Wald", MEM = TRUE)
test1.3

cbind(test1.1, test1.3)

identical(test1.1, test1.3)


summary(m1.2)
test1.2 <- tidyoutput(m1.2, meth = "Wald", MEM = FALSE)
test1.2

cbind(test1.1, test1.2)


tt <- getME(m1.1,"theta")
ll <- getME(m1.1,"lower")
min(tt[ll==0])



m1 <- glmer(discharge_assess_bin ~ 1 + gender + IMD_quintile_all + agecat + CVD_history + 
            musc_skel_history + SMI_bin_no_missing + anxiety_bin_no_missing + depression_bin_no_missing +
            CAT_score_init_cat + (1 | org_code),
            data = dat, family = "binomial")

tidym1 <- tidyoutput(m1, meth = "Wald", MEM = TRUE)

vif(m1)

tidym1 <- tidyoutput(m1, meth = "Wald", MEM = TRUE)[-1, ] %>% unite(CI, -est, sep = " to ")
colnames(tidym1) <- c("Adjusted_estimate", "Adjusted_95_CI")

m1.1 <- glmer(discharge_assess_bin ~ 1 + gender + (1 | org_code),
            data = dat, family = "binomial")

m1.2 <- glmer(discharge_assess_bin ~ 1 + IMD_quintile_all + (1 | org_code),
              data = dat, family = "binomial")

m1.3 <- glmer(discharge_assess_bin ~ 1  + agecat + (1 | org_code),
              data = dat, family = "binomial")

m1.4 <- glmer(discharge_assess_bin ~ 1 +  CVD_history + (1 | org_code),
              data = dat, family = "binomial")

m1.5 <- glmer(discharge_assess_bin ~ 1 + musc_skel_history + (1 | org_code),
              data = dat, family = "binomial")

m1.6 <- glmer(discharge_assess_bin ~ 1 + SMI_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m1.7 <- glmer(discharge_assess_bin ~ 1 + anxiety_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m1.8 <- glmer(discharge_assess_bin ~ 1 + depression_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m1.9 <- glmer(discharge_assess_bin ~ 1 + CAT_score_init_cat + (1 | org_code),
              data = dat, family = "binomial")


tidym1un <- rbind(tidyoutput(m1.1, meth = "Wald", MEM = TRUE)[-1, ],    
                  tidyoutput(m1.2, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m1.3, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m1.4, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m1.5, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m1.6, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m1.7, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m1.8, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m1.9, meth = "Wald", MEM = TRUE)[-1, ])

tidym1un <- tidym1un %>% unite(CI, -est, sep = " to ") 
colnames(tidym1un) <- c("Unadjusted_estimate", "Unadjusted_95_CI")

discharge_assess_OR_table <- cbind(tidym1, tidym1un)
discharge_assess_OR_table <- discharge_assess_OR_table %>% mutate(Variable = as.character(row.names(.))) %>%
                             filter(Variable %!in% c("genderNot recorded/Preferred not to say",
                                                     "genderOther", "genderTransgender")) %>%
                             select(Variable, Unadjusted_estimate, Unadjusted_95_CI,
                                    Adjusted_estimate, Adjusted_95_CI)


discharge_assess_OR_table <- insertRow(discharge_assess_OR_table, 
                                       newrow = c("IMD_quintile_all1", 1, "-", 1, "-"),
                                       2)

discharge_assess_OR_table <- insertRow(discharge_assess_OR_table, 
                                       newrow = c("agecat65-74", 1, "-", 1, "-"),
                                       11)

discharge_assess_OR_table <- insertRow(discharge_assess_OR_table, 
                                       newrow = c("CAT_score_init_cat21-30", 1, "-", 1, "-"),
                                       21)

discharge_assess_OR_table


# write.csv(discharge_assess_OR_table,
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/discharge_assess_OR_table.csv",
#           row.names = FALSE)


# summary(dat$discharge_assess_no_reason)

summary(dat$test_value_init)

dat$baseline_test_value_cat <- cut(dat$test_value_init, breaks = 10)

lintestOR(dat, depvar = "MCID_exercise_cat", indvar = "baseline_test_value_cat")

# Analysis 2: Likelihood of making MCID walking test

m2 <- glmer(MCID_exercise_cat ~ 1 + gender + IMD_quintile_all + agecat + CVD_history + 
              musc_skel_history + SMI_bin_no_missing + anxiety_bin_no_missing + depression_bin_no_missing +
              CAT_score_init_cat + (1 | org_code),
            data = dat, family = "binomial")

summary(m2)

vif(m2)







tidym2 <- tidyoutput(m2, meth = "Wald", MEM = TRUE)[-1, ] %>% unite(CI, -est, sep = " to ")
colnames(tidym2) <- c("Adjusted_estimate", "Adjusted_95_CI")

m2.1 <- glmer(MCID_exercise_cat ~ 1 + gender + (1 | org_code),
              data = dat, family = "binomial")

m2.2 <- glmer(MCID_exercise_cat ~ 1 + IMD_quintile_all + (1 | org_code),
              data = dat, family = "binomial")

m2.3 <- glmer(MCID_exercise_cat ~ 1  + agecat + (1 | org_code),
              data = dat, family = "binomial")

m2.4 <- glmer(MCID_exercise_cat ~ 1 +  CVD_history + (1 | org_code),
              data = dat, family = "binomial")

m2.5 <- glmer(MCID_exercise_cat ~ 1 + musc_skel_history + (1 | org_code),
              data = dat, family = "binomial")

m2.6 <- glmer(MCID_exercise_cat ~ 1 + SMI_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m2.7 <- glmer(MCID_exercise_cat ~ 1 + anxiety_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m2.8 <- glmer(MCID_exercise_cat ~ 1 + depression_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m2.9 <- glmer(MCID_exercise_cat ~ 1 + CAT_score_init_cat + (1 | org_code),
              data = dat, family = "binomial")


tidym2un <- rbind(tidyoutput(m2.1, meth = "Wald", MEM = TRUE)[-1, ],    
                  tidyoutput(m2.2, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m2.3, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m2.4, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m2.5, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m2.6, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m2.7, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m2.8, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m2.9, meth = "Wald", MEM = TRUE)[-1, ])

tidym2un <- tidym2un %>% unite(CI, -est, sep = " to ") 
colnames(tidym2un) <- c("Unadjusted_estimate", "Unadjusted_95_CI")

MCID_exercise_OR_table <- cbind(tidym2, tidym2un)
MCID_exercise_OR_table <- MCID_exercise_OR_table %>% mutate(Variable = as.character(row.names(.))) %>%
  filter(Variable %!in% c("genderNot recorded/Preferred not to say",
                          "genderOther", "genderTransgender")) %>%
  select(Variable, Unadjusted_estimate, Unadjusted_95_CI,
         Adjusted_estimate, Adjusted_95_CI)

MCID_exercise_OR_table <- insertRow(MCID_exercise_OR_table, 
                                       newrow = c("IMD_quintile_all1", 1, "-", 1, "-"),
                                       2)

MCID_exercise_OR_table <- insertRow(MCID_exercise_OR_table, 
                                       newrow = c("agecat65-74", 1, "-", 1, "-"),
                                       11)

MCID_exercise_OR_table <- insertRow(MCID_exercise_OR_table, 
                                       newrow = c("CAT_score_init_cat21-30", 1, "-", 1, "-"),
                                       21)

MCID_exercise_OR_table

# write.csv(MCID_exercise_OR_table,
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/MCID_exercise_OR_table.csv",
#           row.names = FALSE)




# Analysis 3: Likelihood of making MCID health status

m3 <- glmer(MCID_CAT_CRQ_cat ~ 1 + gender + IMD_quintile_all + agecat + CVD_history + 
              musc_skel_history + SMI_bin_no_missing + anxiety_bin_no_missing + depression_bin_no_missing +
              CAT_score_init_cat + (1 | org_code),
            data = dat, family = "binomial")

summary(m3)

vif(m3)


tidym3 <- tidyoutput(m3, meth = "Wald", MEM = TRUE)[-1, ] %>% unite(CI, -est, sep = " to ")
colnames(tidym3) <- c("Adjusted_estimate", "Adjusted_95_CI")

m3.1 <- glmer(MCID_CAT_CRQ_cat ~ 1 + gender + (1 | org_code),
              data = dat, family = "binomial")

m3.2 <- glmer(MCID_CAT_CRQ_cat ~ 1 + IMD_quintile_all + (1 | org_code),
              data = dat, family = "binomial")

m3.3 <- glmer(MCID_CAT_CRQ_cat ~ 1  + agecat + (1 | org_code),
              data = dat, family = "binomial")

m3.4 <- glmer(MCID_CAT_CRQ_cat ~ 1 +  CVD_history + (1 | org_code),
              data = dat, family = "binomial")

m3.5 <- glmer(MCID_CAT_CRQ_cat ~ 1 + musc_skel_history + (1 | org_code),
              data = dat, family = "binomial")

m3.6 <- glmer(MCID_CAT_CRQ_cat ~ 1 + SMI_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m3.7 <- glmer(MCID_CAT_CRQ_cat ~ 1 + anxiety_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m3.8 <- glmer(MCID_CAT_CRQ_cat ~ 1 + depression_bin_no_missing + (1 | org_code),
              data = dat, family = "binomial")

m3.9 <- glmer(MCID_CAT_CRQ_cat ~ 1 + CAT_score_init_cat + (1 | org_code),
              data = dat, family = "binomial")


tidym3un <- rbind(tidyoutput(m3.1, meth = "Wald", MEM = TRUE)[-1, ],    
                  tidyoutput(m3.2, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m3.3, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m3.4, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m3.5, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m3.6, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m3.7, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m3.8, meth = "Wald", MEM = TRUE)[-1, ],
                  tidyoutput(m3.9, meth = "Wald", MEM = TRUE)[-1, ])

tidym3un <- tidym3un %>% unite(CI, -est, sep = " to ") 
colnames(tidym3un) <- c("Unadjusted_estimate", "Unadjusted_95_CI")

MCID_health_status_OR_table <- cbind(tidym3, tidym3un)
MCID_health_status_OR_table <- MCID_health_status_OR_table %>% mutate(Variable = as.character(row.names(.))) %>%
  filter(Variable %!in% c("genderNot recorded/Preferred not to say",
                          "genderOther", "genderTransgender")) %>%
  select(Variable, Unadjusted_estimate, Unadjusted_95_CI,
         Adjusted_estimate, Adjusted_95_CI)

MCID_health_status_OR_table <- insertRow(MCID_health_status_OR_table, 
                                    newrow = c("IMD_quintile_all1", 1, "-", 1, "-"),
                                    2)

MCID_health_status_OR_table <- insertRow(MCID_health_status_OR_table, 
                                    newrow = c("agecat65-74", 1, "-", 1, "-"),
                                    11)

MCID_health_status_OR_table <- insertRow(MCID_health_status_OR_table, 
                                    newrow = c("CAT_score_init_cat21-30", 1, "-", 1, "-"),
                                    21)

MCID_health_status_OR_table

# write.csv(MCID_health_status_OR_table,
#           "Z:/Group_work/PS_AA/PR/PR clinical 2019/Data/tidyData/MCID_health_status_OR_table.csv",
#           row.names = FALSE)




# Plots

summary(dat$ref_location)


# Do this to get the fonts in.

font_import() # This I don't know if you have to do it every time if you want to import calibri.
loadfonts(device = "win")




# non-AECOPD referall to start

datKMref_to_start_non_AECOPD <- dat %>% filter(ref_location %in% c("Primary/Community - stable COPD",
                                                                   "Secondary Care - stable COPD", 
                                                                   "Self-referral"))
datKMref_to_start_non_AECOPD$seen <- 1

summary(datKMref_to_start_non_AECOPD$ref_to_start_days)

datKMref_to_start_non_AECOPD %>% select(ref_to_start_days) %>% filter(!is.na(ref_to_start_days)) %>% 
  arrange(ref_to_start_days) %>% tail()

dat %>% select(ref_date, start_date, assess_date, ref_to_start_days, assess_to_start_days) %>% 
  filter(!is.na(ref_to_start_days)) %>% 
  arrange(ref_to_start_days) %>% tail(20)
  
survfit(Surv(datKMref_to_start_non_AECOPD$ref_to_start_days, datKMref_to_start_non_AECOPD$seen) ~ 
          datKMref_to_start_non_AECOPD$country, data = datKMref_to_start_non_AECOPD) %>% 
  plot_survfit(ci = TRUE, legend.title = "Country", xmax = 365) + #, xbreaks = seq(0, 48, 6)) + 
  labs(x = "Time (days)", y = "Percentage of patients referred with stable COPD who have started PR (%)") +
  theme(text=element_text(size=13, family="Calibri"))
                          #       family="CM Roman"))
                          #       family="TT Times New Roman"))
                          #       family="Sans"))
                          # family="Serif"))

library(ggplot2)





# Initial assessment to discharge

datKM_assess_to_discharge <- dat %>% filter(!is.na(discharge_assess_bin))


datKM_assess_to_discharge$discharge_assess_bin_num <- NA
datKM_assess_to_discharge$discharge_assess_bin_num[datKM_assess_to_discharge$discharge_assess_bin == "No"] <- 0
datKM_assess_to_discharge$discharge_assess_bin_num[datKM_assess_to_discharge$discharge_assess_bin == "Yes"] <- 1

# If we do this, those that didn't have a discharge assessment are included:
# 
datKM_assess_to_discharge$assess_to_discharge_days[is.na(datKM_assess_to_discharge$assess_to_discharge) == TRUE] <- 1000

# Otherwise, hash this out if we don't want them included and unhash the bit below

datKM_assess_to_discharge <- datKM_assess_to_discharge %>% filter(discharge_assess_bin == "Yes")



# # # # datKM_assess_to_discharge$seen <- 1

dat %>% select(assess_to_discharge_days) %>% arrange(assess_to_discharge_days) %>% 
  filter(!is.na(assess_to_discharge_days)) %>%tail(5)


datKM_assess_to_discharge1 <- filter(datKM_assess_to_discharge, assess_date < "2019-04-01")

survfit(Surv(datKM_assess_to_discharge1$assess_to_discharge_days, datKM_assess_to_discharge1$discharge_assess_bin_num) ~ 
          datKM_assess_to_discharge1$country, data =  datKM_assess_to_discharge1) %>% 
  plot_survfit(ci = TRUE, legend.title = "Country") + # , xmax = 90) + #, xbreaks = seq(0, 48, 6)) + 
  labs(x = "Time (days)", y = "Percentage of patients who have been discharged from PR following assessment (%)") +
  theme(text=element_text(size=13, family="Calibri"), panel.grid.major.x = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5), 
        panel.grid.minor.x = element_line(colour="grey", size=0.5),
        panel.grid.minor.y = element_line(colour="grey", size=0.5))




survfit(Surv(datKM_assess_to_discharge$assess_to_discharge_days, datKM_assess_to_discharge$discharge_assess_bin_num) ~ 
          datKM_assess_to_discharge$country, data = datKM_assess_to_discharge) %>% 
  plot_survfit(ci = TRUE, legend.title = "Country") + # , xmax = 90) + #, xbreaks = seq(0, 48, 6)) + 
  labs(x = "Time (days)", y = "Percentage of patients who have been discharged from PR following assessment (%)") +
  theme(text=element_text(size=13, family="Calibri"), panel.grid.major.x = element_line(colour="grey", size=0.5),
        panel.grid.major.y = element_line(colour="grey", size=0.5), 
        panel.grid.minor.x = element_line(colour="grey", size=0.5),
        panel.grid.minor.y = element_line(colour="grey", size=0.5))


summary(datKM_assess_to_discharge$assess_to_discharge_days)


              
# Is there an issue including those marked as 'centre-based' along with those marked as 'both'in the same
# group
dat %>% filter(PR_location == "Centre-based") %>% select(scheduled_sess_centre_no) %>% summary()
dat %>% filter(PR_location == "Both") %>% select(scheduled_sess_centre_no) %>% summary()


dat %>% filter(PR_location == "Centre-based") %>% select(scheduled_sess_centre_no) %>% pull() %>% hist()
dat %>% filter(PR_location == "Both") %>% select(scheduled_sess_centre_no) %>% pull() %>% hist()
# No doesn't look like it

# Is there an issue including those marked as 'home-based' along with those marked as 'both' in the same
# group
dat %>% filter(PR_location == "Home-based") %>% select(scheduled_sess_home_no) %>% summary()
dat %>% filter(PR_location == "Both") %>% select(scheduled_sess_home_no) %>% summary()


dat %>% filter(PR_location == "Home-based") %>% select(scheduled_sess_home_no) %>% pull() %>% hist()
dat %>% filter(PR_location == "Both") %>% select(scheduled_sess_home_no) %>% pull() %>% hist()


# wilcox.test((dat %>% filter(PR_location == "Home-based") %>% select(scheduled_sess_home_no) %>% pull()),
#             (dat %>% filter(PR_location == "Both") %>% select(scheduled_sess_home_no) %>% pull()),
#              alternative = "two.sided")
# Yes - difference in number of sessions scheduled for those who are just home-based and those who are both

histnorm(dat$CRQ_dyspnoea_init) 
hist(dat$ESWT_value_init)
curve(dnorm(dat, mean=mean(dat$age, na.rm = TRUE), sd=sd(dat$age, na.rm = TRUE)))


hist(dat$CAT_score_init) # Normal
histnorm(dat$CAT_score_init)
hist(dat$CRQ_dyspnoea_init) # Skewed
hist(dat$CRQ_fatigue_init) # Skewed
histnorm(dat$CRQ_emotion_init) # Skewed
hist(dat$CRQ_mastery_init) # Skewed



