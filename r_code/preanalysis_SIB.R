library(readxl)
library(dplyr)
library(tableone)
library(survminer)
library(survival)
library(MASS)
library(DescTools)
library(survMisc)
#data from mail on 30 August:
dat <- read_excel("data/210830_210721_210526 SinonasaltabelleSIB_for stats_short_2_VMEI_sent to Sandar.xlsx",sheet = "Daten")

# Tumor code:
dat$TumorCode <- dat$TumorCode_simplyfied...9
dat$TumorCode <- ifelse(is.na(dat$TumorCode), "other", dat$TumorCode)
dat$TumorCode <- ifelse(dat$TumorCode==0, "epithelia", dat$TumorCode)
dat$TumorCode <- ifelse(dat$TumorCode==1, "mesenchymal", dat$TumorCode)

# calculate follow-up time
dat$followup <- as.numeric(dat$DateLast_FU-dat$FirstRT)

# configure data types
dat$Stage4_y_n <- factor(dat$Stage4_y_n)
dat$Epistaxis <- factor(dat$Epistaxis)
dat$Tumor_undifferentiated <- factor(dat$Tumor_undifferentiated)

names(dat)

#Events: Progresion and Death
table(dat$Protocol)

dat %>% group_by(Protocol) %>%
  count(Progression_or_Death_y_n)
dat %>% group_by(Protocol) %>%
  count(Progression_only_y_n)
dat %>% group_by(Protocol) %>%
  count(Dead_y_n)

# Calculating PFS OS stratified by Protocol
# ------------------
dat_regular <- dat %>% filter(Protocol==0)
dat_boost <- dat %>% filter(Protocol==1)

#PFS: regular, boost, all
MedianCI(dat_regular$PFS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat_boost$PFS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat$PFS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)

#OS: regular, boost, all
MedianCI(dat_regular$OS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat_boost$OS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat$OS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)

# Calculating PFS OS stratified by Stage I-III vs Stage IV
# ------------------
dat_s1_3 <- dat %>% filter(Stage4_y_n==0)
dat_s4 <- dat %>% filter(Stage4_y_n==1)
# PFS
MedianCI(dat_s1_3$PFS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat_s4$PFS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
# OS
MedianCI(dat_s1_3$OS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat_s4$OS,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)


# Calculating proportions of dogs alive at 1 year and 2 year (OS)
# ------------------
dat %>%
  group_by(Protocol) %>%
  count(Dead_y_n)
##OS 1 year:
dat %>%
  group_by(Protocol) %>%
  count(OS>365.25)  %>%
  mutate(prop = n / sum(n), 
         lower = lapply(n, prop.test, n = sum(n)), 
         upper = sapply(lower, function(x) x$conf.int[2]), 
         lower = sapply(lower, function(x) x$conf.int[1]))
#26 dogs alive at the end of 1 year. 
# 8 survive, 18 died eventually
##OS 2 year:
dat %>%
  group_by(Protocol) %>%
  count(OS>365.25*2)  %>%
  mutate(prop = n / sum(n), 
         lower = lapply(n, prop.test, n = sum(n)), 
         upper = sapply(lower, function(x) x$conf.int[2]), 
         lower = sapply(lower, function(x) x$conf.int[1]))



# Calculating proportions of dogs only progression-free at 1 year and 2 year (PFS)
#I would also like to know the proportion of dogs free of progression at 1 and at 2 years. This is what some studies looked at in the past. Since for PFS the endpoint was either progression or death, you cannot do the analysis from the table I gave you (=from PFS and the variable progression/death yes/no). I therefore added another column in the excel file in the attachment (highlighted in yellow). <-- Progression_only_y_n

dat %>%
  group_by(Protocol) %>%
  count(Progression_only_y_n)

dat <- dat %>%
  mutate(PFSonly = DateProgr_clin - FirstRT,
         PFSonly = if_else(is.na(PFSonly),DateLast_FU - FirstRT , PFSonly)) # PFSonly <- either only PFS or max. follow up time
##only PFS 1 year:
dat %>%
  group_by(Protocol) %>%
  count(PFSonly>365.25)  %>%
  mutate(prop = n / sum(n), 
         lower = lapply(n, prop.test, n = sum(n)), 
         upper = sapply(lower, function(x) x$conf.int[2]), 
         lower = sapply(lower, function(x) x$conf.int[1]))
##only PFS 2 year:
dat %>%
  group_by(Protocol) %>%
  count(PFSonly>365.25*2)  %>%
  mutate(prop = n / sum(n), 
         lower = lapply(n, prop.test, n = sum(n)), 
         upper = sapply(lower, function(x) x$conf.int[2]), 
         lower = sapply(lower, function(x) x$conf.int[1]))



#relevant parameters for table 1
#-------------------
myVars <- c("Age", 
            "Weight",
            "Sex",
            "Stage", 
            "OS",
            "PFS", 
            "followup", 
            "TumorCode",
            "Tumor_undifferentiated")
catVars <- c("Sex", "TumorCode","Tumor_undifferentiated", "Stage")

tab2 <- CreateTableOne(strata = "Protocol", vars = myVars, factorVars = catVars, data = dat,test = TRUE )
print(tab2)
CreateTableOne(strata = "Protocol", vars = myVars, factorVars = catVars, data = dat,test = TRUE )
CreateTableOne(vars = myVars, data = dat,  factorVars = catVars,test = TRUE )


#Follow-up times
MedianCI(dat$followup,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat_regular$followup,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)
MedianCI(dat_boost$followup,
         conf.level = 0.95,
         na.rm = FALSE,
         method = "exact",
         R = 10000)


