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

#Fill the Censor

for (i in 1:length(df$DateOfDeath)){
  if (is.na(df$SurvivalTimeMonth[i])) {
    df$SurvivalTimeMonth[i]<-
      interval(df$DateofDiagnosis[i], df$DateofDiagnosis[length(df$AgeAtDiagnosis)])%/%months(1)
  }
}

# descriptive analysis ----------------------------------------------------

mean(df$AgeAtDiagnosis)
sd(df$AgeAtDiagnosis)
table<-lapply(df[,2:4], table)
lapply(table, prop.table)


# univariable  and cox regression ----------------------------------------------

dependent_os<-'Surv(SurvivalTimeMonth, Status==1)'
dependent_cfs<-'Surv(SurvivalTimeMonth, StatusCancer==1)'

explanatory<-c('AgeAtDiagnosis', 'Gender', 'Differentiation',
               'PrimaryTumourSite')

suros<-df %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE, digits=c(3,3,4)) 

surcfs<-df %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE, digits=c(3,3,4)) 


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


# Calibration of the nomogram----------------------------------------------------------------

#the calibration of overall survival
timeinc = 12
overallcox<-cph(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
               x=T, y=T, surv=T, data = df, time.inc = timeinc)
overallcal1<-calibrate(overallcox, cmethod="KM", method='boot', u=12, m=500, B=1000, pr=F, subtitles=F)

timeinc = 36
overallcox<-cph(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
                x=T, y=T, surv=T, data = df, time.inc = timeinc)
overallcal3<-calibrate(overallcox, cmethod="KM", method='boot', u=36, m=500, B=1000, pr=F, subtitles=F)

timeinc = 60
overallcox<-cph(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
                x=T, y=T, surv=T, data = df, time.inc = timeinc)
overallcal5<-calibrate(overallcox, cmethod="KM", method='boot', u=60, m=500, B=1000, pr=F, subtitles=F)

timeinc = 120
overallcox<-cph(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
                x=T, y=T, surv=T, data = df, time.inc = timeinc)
overallcal10<-calibrate(overallcox, cmethod="KM", method='boot', u=120, m=500, B=1000, pr=F, subtitles=F)


png('overall survival calibration.png', width = 16,
    height    = 4,
    units     = "in",
    res       = 1200,
    pointsize = 10)
par(mfcol=c(1,4))

plot(overallcal1, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual one-year overall survival', mgp=c(2,0,0), cex.lab=2)

plot(overallcal3, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual three-year overall survival', mgp=c(2,0,0), cex.lab=2)

plot(overallcal5, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual five-year overall survival', mgp=c(2,0,0), cex.lab=2)

plot(overallcal10, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual ten-year overall survival', mgp=c(2,0,0), cex.lab=2)

dev.off()
                    
#calibration of the cancer-free survival

timeinc = 12
cancercox<-cph(Surv(SurvivalTimeMonth, StatusCancer==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
               x=T, y=T, surv=T, data = df, time.inc = timeinc)
cancercal1<-calibrate(cancercox, cmethod="KM", method='boot', u=12, m=500, B=200, pr=F, subtitles=F)

timeinc = 36
cancercox<-cph(Surv(SurvivalTimeMonth, StatusCancer==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
               x=T, y=T, surv=T, data = df, time.inc = timeinc)
cancercal3<-calibrate(cancercox, cmethod="KM", method='boot', u=36, m=500, B=200, pr=F, subtitles=F)

timeinc = 60
cancercox<-cph(Surv(SurvivalTimeMonth, StatusCancer==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
               x=T, y=T, surv=T, data = df, time.inc = timeinc)
cancercal5<-calibrate(cancercox, cmethod="KM", method='boot', u=60, m=500, B=200, pr=F, subtitles=F)

timeinc = 120
cancercox<-cph(Surv(SurvivalTimeMonth, StatusCancer==1)~AgeAtDiagnosis+Sex+PrimaryTumourSite+Differentiation, 
               x=T, y=T, surv=T, data = df, time.inc = timeinc)

cancercal10<-calibrate(cancercox, cmethod="KM", method='boot', u=timeinc, m=500, B=200, pr=F, subtitles=F)

png('cancer-specific survival calibration.png', width = 16,
    height    = 4,
    units     = "in",
    res       = 1200,
    pointsize = 10)
par(mfcol=c(1,4))
plot(cancercal1, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual one-year cancer-specific survival', mgp=c(2,0,0), cex.lab=2)


plot(cancercal3, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual three-year cancer-specific survival', mgp=c(2,0,0), cex.lab=2)


plot(cancercal5, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual five-year cancer-specific survival', mgp=c(2,0,0), cex.lab=2)

plot(cancercal10, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'nomogram predicted probability',
     ylab='actual ten-year cancer-specific survival', mgp=c(2,0,0), cex.lab=2)
dev.off()




# Random survival forest --------------------------------------------------
#fill the censors

r_fit<-ranger(Surv(SurvivalTimeMonth, Status==1)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation,
              data = df)


r_fit$prediction.error

# recursive partition analysis ------------------------------------------
df<-df[which(df$SurvivalTimeMonth!=0),]
overallrpa <-rpart(Surv(SurvivalTimeMonth, Status)~AgeAtDiagnosis+Gender+PrimaryTumourSite+Differentiation,
                   data = df,x=T,y=T)
summary(overallrpa)
overallrpa

