install.packages("vars")
install.packages("TSstudio")
install.packages("DescTools")
library(DescTools)
library(vars)
library(mFilter)
library(tseries)
library(TSstudio)
library(forecast)
library(tidyverse)
library(fpp2)
library(ggplot2)
library(urca)
library(astsa)
library(TSA)
library(tsoutliers)
library(seasonal)
library(tseries)
library(readr)
data <- read_csv("updated.csv")
almor <- ts(data$almond, start = c(1999,11), frequency = 12)
peanr <- ts(data$peanut, start = c(1999,11), frequency = 12)
pecar <- ts(data$pecan, start = c(1999,11), frequency = 12)
walnr <- ts(data$walnut, start = c(1999,11), frequency = 12)

fit_mstl_almo<-mstl(almor)
autoplot(fit_mstl_almo)+ggtitle("mstl decomposition of almond")
#

#
fit_mstl_pean<-mstl(peanr)
autoplot(fit_mstl_pean)+ggtitle("mstl decomposition of peanut")
#
fit_mstl_peca<-mstl(pecar)
autoplot(fit_mstl_peca)+ggtitle("mstl decomposition of pecan")
#
fit_mstl_waln<-mstl(walnr)
autoplot(fit_mstl_waln)+ggtitle("mstl decomposition of walnut")
#

almo<-diff(almor)
pean<-diff(peanr)
peca<-diff(pecar)
waln<-diff(walnr)
overall <-cbind(almo, pean, peca, waln)


adf.test(almo)
adf.test(pean)
adf.test(peca)
adf.test(waln)

train.almo <-subset(almo, end=length(almo)-24)
train.pean <-subset(pean, end=length(pean)-24)
train.peca <-subset(peca, end=length(peca)-24)
train.waln <-subset(waln, end=length(waln)-24)

test.almo <-subset(almo, start=length(almo)-23)
test.pean <-subset(pean, start=length(almo)-23)
test.peca <-subset(peca, start=length(almo)-23)
test.waln <-subset(waln, start=length(almo)-23)




ts_plot(almo)
ts_plot(pean)
ts_plot(peca)
ts_plot(waln)


autoplot(overall)

###########

pp.test(almo)
pp.test(pean)
pp.test(peca)
pp.test(waln)

ndiffs(train.almo, alpha = 0.05, test = c("adf"))
ndiffs(train.pean, alpha = 0.05, test = c("adf"))
ndiffs(train.peca, alpha = 0.05, test = c("adf"))
ndiffs(train.waln, alpha = 0.05, test = c("adf"))

ndiffs(test.almo, alpha = 0.05, test = c("adf"))
ndiffs(test.pean, alpha = 0.05, test = c("adf"))
ndiffs(test.peca, alpha = 0.05, test = c("adf"))
ndiffs(test.waln, alpha = 0.05, test = c("adf"))

adf.test(train.almo)
adf.test(train.pean)
adf.test(train.peca)
adf.test(train.waln)
pp.test(test.almo)
pp.test(test.pean)
pp.test(test.peca)
pp.test(test.waln)


#dif_almo = diff(train.almo, differences = 1)
#dif_pean = diff(train.pean, differences = 1)
#dif_peca = diff(train.peca, differences = 1)
#dif_waln = diff(train.waln, differences = 1)

#dif_almo.test = diff(test.almo, differences = 1)
#dif_pean.test = diff(test.pean, differences = 1)
#dif_peca.test = diff(test.peca, differences = 1)
#dif_waln.test = diff(test.waln, differences = 1)

################################################
#everything is stationary
################################################


#bind############################################

cb <-cbind(train.almo, train.pean, train.peca, train.waln)
colnames(cb) <- cbind("almo","pean", "peca", "waln")
#
tt <-cbind(test.almo, test.pean, test.peca, test.waln)
colnames(tt) <- cbind("almo.t","pean.t", "peca.t", "waln.t")
##################################################
#var model
##################################################
lagselect <- VARselect(cb, lag.max = 15, type = "const")
lagselect$selection
####################################################


#playing with seasonality
####################################################
Var_model1 <- VAR(cb, p = 1, type = "const", season = 12, exog = NULL)
Var_model1
summary(Var_model1)

####Forecast#######################################

forecast <- predict(Var_model1, n.ahead = 24)
forecast

###################################################
#calculating MAPE##################################
##################################################
almo.f<-c(-2.68138554, 069702040, 0.70937161, 3.82980195, 1.49344521, -0.71877891, 3.12615790, 9.69997155, 1.66317335, -1.99109628, -1.94645719, 0.62403281, 0.07874091, -0.80053142, 1.14049749, 3.73649408, 1.50966985, -0.72093827, 3.12630775, 9.69999860, 1.66315866, -1.99109223, -1.94645804, 0.62403296)
pean.f<-c(10.2422642, -0.3122977, 0.1429483, -4.1044193, -2.3706067, 5.7732192, -2.1054015, 1.4474322, -0.9368589, 2.7011173, -4.0977534, -0.1739087, -1.1742972, 1.9344739, -0.2628326, -4.0402309, -2.3785276, 5.7736604, -2.1052700, 1.4473716, -0.9368429, 2.7011140, -4.0977529, -0.1739088)
peca.f<-c(-2.51614378, -5.48370728, 0.89655101, 1.08242890, 6.59968789, -0.08154445, 20.02427039, 0.54759183, 9.43157651, -13.75216592, -1.28770007, 2.28876918, -0.35474032, -7.47633850, 1.81516688, 0.80408307, 6.66498525, -0.09407613, 20.02621311, 0.54737969, 9.43157881, -13.75215903, -1.28770251, 2.28876976)
waln.f<-c(-17.3139984, 8.3538926, 2.1149443, 0.4284578, 1.4464255, 1.3184587, 3.9156583, 1.1788967, -12.0736622, 0.1194217, -2.2711042, 2.6327144, -1.6270768, 4.7522419, 2.8357730, 0.2964303, 1.4673178, 1.3159238, 3.9157827, 1.1789449, -12.0736831, 0.1194271, -2.2711053, 2.6327146)


almo.a<-c(431.7, 431.7, 431.7, 431.7, 431.7, 431.7, 405.7, 405.7, 405.7, 405.7, 405.7, 405.7, 405.7, 405.7, 405.7, 405.7, 405.7, 405.7, 305.2, 305.2, 305.2, 305.2, 305.2)
pean.a<-c(99.1, 105.5, 104.4, 97.9, 104.7, 95.1, 100.5, 90.4, 95.8, 97.8, 104.2, 100.5, 97.9, 99.1, 100.7, 103.7, 100.2, 107.1, 105.9, 103.7, 101.4, 102.9, 100.6, 103.1)
peca.a<-c(712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 712.7, 922.6, 624.8, 624.8, 624.8, 624.8, 624.8)
waln.a<-c(174.9, 191.2, 191.2, 191.2, 191.2, 191.2, 191.2, 254.2, 254.2, 254.2, 254.2, 254.2, 254.2, 254.2, 219.2, 219.2, 219.2, 219.2, 219.2, 191.2, 191.2, 218.4, 184.2, 184.2)


accuracy(almo.f, almo.a)
accuracy(pean.f, pean.a)
accuracy(peca.f, peca.a)
accuracy(waln.f, waln.a)


##############################
#Error decomposition
###############################
fore_nut <- fevd(Var_model1, n.ahead = 24)
fore_nut
plot(fore_nut)

#######################################
#extracting data for DAG model########
#######################################
residuals(Var_model1)
#######################################

#Var_model2 <- VAR(tt, p = 1, type = "const", season = 12, exog = NULL)
#summary(Var_model2)

forecast <- predict(Var_model1, n.ahead = 24)

plot(forecast)
print(forecast)

##################
forecast <- predict(Var_model1, n.ahead = 24)
plot(forecast)
print(forecast)

fanchart(forecast, names = "almo", main = "Fanchart for almo", xlab = "Horizon", ylab = "almo")
fanchart(forecast, names = "pean", main = "Fanchart for pean", xlab = "Horizon", ylab = "pean")
fanchart(forecast, names = "peca", main = "Fanchart for peca", xlab = "Horizon", ylab = "peca")
fanchart(forecast, names = "waln", main = "Fanchart for waln", xlab = "Horizon", ylab = "waln")
forecast
###############################
###############################
###############################
# forecast issue
###########


Serial1 <- serial.test(Var_model1, lags.pt = 12, type = "PT.asymptotic")
Serial1
#no serial correlation from test


Norm1 <- normality.test(Var_model1, multivariate.only = TRUE)
Norm1

####
#didnt pass normaility
#

Stability1 <- stability(Var_model1, type = "OLS-CUSUM")
plot(Stability1)
####
#no structural breaks


#save residuals and estimate cyclic graph

####
#granger, > contemparaneous 
#####

Granger_almo<- causality(Var_model1, cause = "almo")
Granger_almo
#1
Granger_pean <- causality(Var_model1, cause = "pean")
Granger_pean
#2
Granger_peca <- causality(Var_model1, cause = "peca")
Granger_peca
#3
Granger_waln <- causality(Var_model1, cause = "waln")
Granger_waln
#4


######
#impulse
######

almo_irf <- irf(Var_model1, impulse = "almo", response = "almo", n.ahead = 12, boot = TRUE)
plot(almo_irf, ylab = "almo", main = "almo's shock to almo")
pean_irf <- irf(Var_model1, impulse = "almo", response = "pean", n.ahead = 12, boot = TRUE)
plot(pean_irf, ylab = "pean", main = "almo's shock to pean")
peca_irf <- irf(Var_model1, impulse = "almo", response = "peca", n.ahead = 12, boot = TRUE)
plot(peca_irf, ylab = "peca", main = "almo's shock to peca")
waln_irf <- irf(Var_model1, impulse = "almo", response = "waln", n.ahead = 12, boot = TRUE)
plot(waln_irf, ylab = "waln", main = "almo's shock to waln")
print(almo_irf)
print(pean_irf)
print(peca_irf)
print(waln_irf)
#######
#repeating
########
almo_irf <- irf(Var_model1, impulse = "pean", response = "almo", n.ahead = 12, boot = TRUE)
plot(almo_irf, ylab = "almo", main = "pean's shock to almo")
pean_irf <- irf(Var_model1, impulse = "pean", response = "pean", n.ahead = 12, boot = TRUE)
plot(pean_irf, ylab = "pean", main = "pean's shock to pean")
peca_irf <- irf(Var_model1, impulse = "pean", response = "peca", n.ahead = 12, boot = TRUE)
plot(peca_irf, ylab = "peca", main = "pean's shock to peca")
waln_irf <- irf(Var_model1, impulse = "pean", response = "waln", n.ahead = 12, boot = TRUE)
plot(waln_irf, ylab = "waln", main = "pean's shock to waln")
print(almo_irf)
print(pean_irf)
print(peca_irf)
print(waln_irf)
########

almo_irf <- irf(Var_model1, impulse = "peca", response = "almo", n.ahead = 12, boot = TRUE)
plot(almo_irf, ylab = "almo", main = "peca's shock to almo")
pean_irf <- irf(Var_model1, impulse = "peca", response = "pean", n.ahead = 12, boot = TRUE)
plot(pean_irf, ylab = "pean", main = "peca's shock to pean")
peca_irf <- irf(Var_model1, impulse = "peca", response = "peca", n.ahead = 12, boot = TRUE)
plot(peca_irf, ylab = "peca", main = "peca's shock to peca")
waln_irf <- irf(Var_model1, impulse = "peca", response = "waln", n.ahead = 12, boot = TRUE)
plot(waln_irf, ylab = "waln", main = "peca's shock to waln")
print(almo_irf)
print(pean_irf)
print(peca_irf)
print(waln_irf)
#########
almo_irf <- irf(Var_model1, impulse = "waln", response = "almo", n.ahead = 12, boot = TRUE)
plot(almo_irf, ylab = "almo", main = "waln's shock to almo")
pean_irf <- irf(Var_model1, impulse = "waln", response = "pean", n.ahead = 12, boot = TRUE)
plot(pean_irf, ylab = "pean", main = "waln's shock to pean")
peca_irf <- irf(Var_model1, impulse = "waln", response = "peca", n.ahead = 12, boot = TRUE)
plot(peca_irf, ylab = "peca", main = "waln's shock to peca")
waln_irf <- irf(Var_model1, impulse = "waln", response = "waln", n.ahead = 12, boot = TRUE)
plot(waln_irf, ylab = "waln", main = "waln's shock to waln")
print(almo_irf)
print(pean_irf)
print(peca_irf)
print(waln_irf)
#####
#analysis of supply
#####
###########################################################
###########################################################
###########################################################

fore_nut <- fevd(Var_model1, n.ahead = 24)
fore_nut
plot(fore_nut)
#forecast error decomposition 
###
#estimate forecast error decomposition in r
####
#
####
#
####
forecast <- predict(Var_model1, n.ahead = 24)
plot(forecast)
print(forecast)

fanchart(forecast, names = "almo", main = "Fanchart for almo", xlab = "Horizon", ylab = "almo")
fanchart(forecast, names = "pean", main = "Fanchart for pean", xlab = "Horizon", ylab = "pean")
fanchart(forecast, names = "peca", main = "Fanchart for peca", xlab = "Horizon", ylab = "peca")
fanchart(forecast, names = "waln", main = "Fanchart for waln", xlab = "Horizon", ylab = "waln")
forecast

model = sm.tsa.VARMAX(y_train, order=(5), trend='c')
model_result = model.fit(maxiter=1000, disp=False)
model_result.summary()


# not appropriate here, must have same length as p
####
#
# 238 and use 12 model for forecast accuracy
# point forcast and use MAPE
####
#playing with accuracy
#


#type = c("response", "terms"),
#terms = NULL, na.action = na.pass,
#pred.var = res.var/weights, weights = 1, .)



plot(result.mean, type = "l")

MAPE <- mean(result.mean)
MAPE

####################
#adf stat Model ????
# estimate var model diff # done
# save residuals and then do cyclic graph from resid #done?
# it will tell contemporaneous series
# take var and covert impulse response function
# estimate forecast error variance decomposition for each # done
# end
# 2004 bessler dharmasena tea market


