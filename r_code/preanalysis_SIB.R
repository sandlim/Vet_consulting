library(readxl)
library(dplyr)
library(tableone)
library(survminer)
library(survival)
library(MASS)
library(DescTools)

dat_og <- read_excel("data/210526 SinonasaltabelleSIB_for stats_short_2_VMEI.xlsx", sheet = "Daten")
additional_dat <- read_excel("data/210721_210526 SinonasaltabelleSIB_for stats_short_2_VMEI_sent to Sandar.xlsx",sheet = "Daten")
names(dat_og)
names(additional_dat)
dat <- full_join(dat_og, additional_dat)
str(dat)
table(dat$TumorCode_simplyfied...9)
dat$TumorCode <- dat$TumorCode_simplyfied...9
dat$TumorCode <- ifelse(is.na(dat$TumorCode), "other", dat$TumorCode)
dat$TumorCode <- ifelse(dat$TumorCode==0, "epithelia", dat$TumorCode)
dat$TumorCode <- ifelse(dat$TumorCode==1, "mesenchymal", dat$TumorCode)

table(is.na(dat$DateLast_FU))
table(is.na(dat$FirstRT))
dat$followup <- as.numeric(dat$DateLast_FU-dat$FirstRT)

dat$Stage4_y_n <- factor(dat$Stage4_y_n)
dat$Epistaxis <- factor(dat$Epistaxis)
dat$Tumor_undifferentiated <- factor(dat$Tumor_undifferentiated)

#Events
dat %>% group_by(Protocol) %>%
  count(Progression_y_n)
dat %>% group_by(Protocol) %>%
  count(Dead_y_n)
table(dat$Protocol)

##PFS OS stratified by Protocol
dat_regular <- dat %>% filter(Protocol==0)
dat_boost <- dat %>% filter(Protocol==1)
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


##OS 1 year:
dat %>%
  group_by(Protocol) %>%
  count(OS>365)  %>%
  mutate(prop = n / sum(n), 
         lower = lapply(n, prop.test, n = sum(n)), 
         upper = sapply(lower, function(x) x$conf.int[2]), 
         lower = sapply(lower, function(x) x$conf.int[1]))



##PFS OS stratified by Stage
dat_s1_3 <- dat %>% filter(Stage4_y_n==0)
dat_s4 <- dat %>% filter(Stage4_y_n==1)
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


#relevant parameters for table 1


myVars <- c("Age", 
            "Weight",
            "Sex",
            "Stage", 
            "OS","PFS", "followup", 
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


