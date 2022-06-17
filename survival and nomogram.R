library(dplyr)
library(plyr)
library(lubridate)
library(survival)
library(finalfit)
library(survminer)
library(rms)
library(pec)
library(randomForestSRC)
library(rpart)
library(ranger)

# date preprocessing --------------------------------------------------------

df<-read.csv('cancer survival analysis.csv', na.strings = c('NULL','#VALUE!', ' '))
df$Status[df$Status==2]<-0

# Cause of Death simplification -------------------------------------------

df<-df%>%mutate(CauseofDeathSimple= ifelse(CauseofDeath=='NULL', 'alive',
                                           ifelse(CauseofDeath=='NCD', 'non_cancer', 'cancer')))

df<-df %>% mutate(StatusCancer = ifelse(Status==1&CauseofDeathSimple=='cancer',1, 0))


#remove patients younger than 18
df<-df[which(df$AgeAtDiagnosis>=18),]

#Year of Diagnosis and Death
df$DateofDiagnosis<-dmy(df$DateofDiagnosis)
df$DateOfDeath<-dmy(df$DateOfDeath)
df<-df %>% mutate(YearofDiagnosis = format(DateofDiagnosis, format = '%Y'))
df<-df %>% mutate(YearofDeath = format(DateOfDeath, format = '%Y'))


#change factors into number 
df$Differentiation<-factor(df$Differentiation,
                           levels = c('Well Differentiated',
                                             'Moderately Differentiated',
                                             'Poorly Differentiated',
                                             'Undifferentiated',
                                             'Not Stated/Unknown'),
                           labels = c(1:5))

df$PrimaryTumourSite<-factor(df$PrimaryTumourSite, 
                             levels = c('Commissure of lip',
                                        'Floor of mouth',
                                        'Cheek and Vestibule',
                                        'Tongue and Lingual tonsil',
                                        'Retromolar area',
                                        'Mouth',
                                        'Gum','Tonsil',
                                        'Hard and Soft palate',
                                        'Overlapping lesion of lip, oral cavity and pharynx',
                                        'Oropharynx'),
                             labels = c(1:11))

df$Gender<-factor(df$Gender, levels = c('Female', 'Male'), labels = c(1, 2))


# descriptive analysis ----------------------------------------------------

mean(df$AgeAtDiagnosis)
sd(df$AgeAtDiagnosis)
table<-lapply(df[,2:4], table)
lapply(table, prop.table)


# univariable cox regression ----------------------------------------------

agecox<-coxph(Surv(SurvivalTimeMonth, Status)~AgeAtDiagnosis, data = df)
summary(agecox)



gendercox<-coxph(Surv(SurvivalTimeMonth,Status)~Gender, data = df)
summary(gendercox)

diffcox<-coxph(Surv(SurvivalTimeMonth, Status)~Differentiation, data = df)
summary(diffcox)

sitecox<-coxph(Surv(SurvivalTimeMonth, Status)~PrimaryTumourSite, data = df)


#multivariable cox regression
dependent_os<-'Surv(SurvivalTimeMonth, Status==1)'
dependent_cfs<-'Surv(SurvivalTimeMonth, StatusCancer==1)'

explanatory<-c('AgeAtDiagnosis', 'Gender', 'Differentiation',
               'PrimaryTumourSite')

surtest<-df %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE, digits=c(3,3,4)) %>% 
  rename("Overall survival" = label) %>% 
  rename(" " = levels) %>% 
  rename("  " = all)

df %>%
  hr_plot(dependent_os, explanatory, dependent_label = "Survival")

# survival analysis -------------------------------------------------------


coxmodel<-coxph(Surv(SurvivalTimeMonth, StatusCancer==1)~AgeAtDiagnosis+Gender+
                  Differentiation+PrimaryTumourSite,data = df)
summary(coxmodel)

# Cox proportional hazard regression --------------------------------------

#cox regression and nomogram

ss <- c(0.05,0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

ddist<-datadist(df)
options(datadist='ddist')

timeinc = 12
overallcox<-cph(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation, 
                x=T, y=T, surv=T, data = df, time.inc = timeinc)
overallcox

overallnom<-nomogram(overallcox, fun = list(surv1y,function(x) surv(36, lp=x),
                                            function(x) surv(60, lp=x), function(x) surv(120, lp=x)),
                     funlabel=c("1 year survival", "3 year survival", "5 year survival", '10 year survival'),
                     lp=F, fun.at = ss)

plot(overallnom)

validate(overallnom, B=150, dxy=T)
concordance(overallcox)

#cancer specific survival analysis
cancercox<-cph(Surv(SurvivalTimeMonth, StatusCancer==1)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation, 
               x=T, y=T, surv=T, data = df, time.inc = timeinc)


cancernom<-nomogram(cancercox, fun = list(surv1y,function(x) surv(36, lp=x),
                                          function(x) surv(60, lp=x), function(x) surv(120, lp=x)),
                    funlabel=c("1 year cancer-free survival", "3 year cancer-free survival", 
                               "5 year cancer-free survival", '10 year cancer-free survival'),
                    lp=F, fun.at = ss)


validate(cancernom, B=150, dxy=T)
concordance(cancercox)


png('nomogram overall survival.png', width=700, height=1000, unit="px")
plot(overallnom)
dev.off()

png('nomogram cancer-free survival.png', width=700, height=1000, unit="px")
plot(cancernom)
dev.off()


# nomogram and its calibration----------------------------------------------------------------

timeinc = 120
overallcox<-cph(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation, 
                x=T, y=T, surv=T, data = df, time.inc = timeinc)

overallcal1<-calibrate(overallcox, cmethod="KM", method='boot', u=timeinc, m=500, B=200, pr=F, subtitles=F)

par(mar=c(10,6,5,5), cex=1.0)

plot(overallcal1, xlim = c(0, 1), ylim = c(0, 1),
     xlab = '',
     ylab='actual ten year overall survival')
title(xlab = 'nomogram predicted probability', line=2)


timeinc = 120
cancercox<-cph(Surv(SurvivalTimeMonth, StatusCancer==1)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation, 
               x=T, y=T, surv=T, data = df, time.inc = timeinc)


cancercal1<-calibrate(cancercox, cmethod="KM", method='boot', u=timeinc, m=500, B=200, pr=F, subtitles=F)

par(mar=c(10,6,5,5), cex=1.0)

plot(cancercal1, xlim = c(0, 1), ylim = c(0, 1),
     xlab = '',
     ylab='actual ten year cancer-free survival')
title(xlab = 'nomogram predicted probability', line=2)


# Random survival forest --------------------------------------------------
#fill the censors
df2<-df

for (i in 1:length(df2$DateOfDeath)){
  if (is.na(df2$SurvivalTimeMonth[i])) {
    df2$SurvivalTimeMonth[i]<-
      interval(df2$DateofDiagnosis[i], df2$DateofDiagnosis[length(df2$AgeAtDiagnosis)])%/%months(1)
  }
}

r_fit<-ranger(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation,
              data = df2)


r_fit$prediction.error

# recursive partition analysis ------------------------------------------
df2<-df2[which(df2$SurvivalTimeMonth!=0),]
overallrpa <-rpart(Surv(SurvivalTimeMonth, Status)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation,
                   data = df2,x=T,y=T)
summary(overallrpa)
overallrpa

