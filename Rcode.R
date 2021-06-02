library(fpp2)
library(ggplot2)
library(urca)
library(astsa)
library(TSA)
library(tsoutliers)
library(seasonal)
library(tseries)
library(readr)
train <- read_csv("Almond.csv")
train<-na.omit(train)
na.omit(train)
as.numeric(train$ppi)
#as.Date(train$date, "%m/%d/%Y")

#myts=ts(train[,-1], start=c(2003,1), frequency = 12)
myts=ts(train[,-1], start=c(1991,12), frequency = 12)

#create ts object
tsdisplay(myts)
#plot acf and pacf

autoplot(myts)
adf.test(myts)
#data is non-stationary
myts %>% ur.kpss() %>% summary()
# we need diffrerence (unit root test)
ndiffs(myts)
# one diff is required
nsdiffs(myts)
#non seasonal data
dmyts<-diff(myts)

#repeating this process
autoplot(dmyts)
adf.test(dmyts)
#data is stationary
dmyts %>% ur.kpss() %>% summary()
# we need diffrerence (unit root test)
ndiffs(dmyts)
# diff is not required
nsdiffs(dmyts)
#non seasonal data


#############################
# creating training and test data with windows()
#############################

myts.train<-window(dmyts,end=c(2018,10))
myts.test<-window(dmyts,start=c(2018,11))

#####plot
autoplot(myts) +
  autolayer(myts.train, series="Training") +
  autolayer(myts.test, series="Test")


###################
###################
#transformation
###################
###################
lambda <- BoxCox.lambda(myts)
print(lambda)
autoplot(BoxCox(myts,lambda))
#A value of lambda=0 is multiplicative and lambda=1 is additive

plot(decompose(myts))
plot(decompose(myts, type = "multiplicative"))

#SEATS error
seas_ts<-seas(myts)
autoplot(seas_ts)+ggtitle("SEATS decomposition")
#x11

fitdecomp<-seas(myts,x11="")
myts %>% seas(x11="") -> fitdecomp
autoplot(fitdecomp) +
  ggtitle("X11 decomposition")
#stl error
myts %>%
  stl(t.window=13, s.window="periodic", robust=TRUE) %>%
  autoplot()+ggtitle("stl decomposition")
#mstl
fit_mstl<-mstl(myts, lambda = -0.293802)
autoplot(fit_mstl)+ggtitle("mstl decomposition")
#Should I use lambda?
fit_mstl<-mstl(myts)
autoplot(fit_mstl)+ggtitle("mstl decomposition of pean")
#
fit_mstl <- mstl(myts.train, t.window=13, s.window="periodic",
           robust=TRUE)

fit_mstl %>% seasadj() %>% naive() %>%
  autoplot() + ylab("New orders index") +
  ggtitle("Naive forecasts of seasonally adjusted data")

naive<-forecast(fit_mstl, method="naive", h=24)
autoplot(naive) + ylab("New orders index")

accuracy(naive,myts.test)

############################
######exponential smoothing
############################
fcast <- stlf(myts.train, method='ets', h=24)
fcast
autoplot(fcast)+ggtitle("Exponential Smoothing (Figure 6)")+ylab("PPI")
accuracy(fcast,myts.test)
#ETS(M,N,N): simple exponential smoothing 
# with multiplicative errors
##############################
#chain rule of forecast check that!!!!!
## check in R
########################

########################
#benchmark forecast/
########################
########################
mean_mod <- rwf(fit_mstl, h=24)
naive_mod <- naive(fit_mstl, h=24)
snaive_mod <- snaive(fit_mstl, h=24)
autoplot(subset(myts, end = 530)) +
  autolayer(mean_mod, PI=FALSE, series="mean") +
  autolayer(naive_mod, PI=FALSE, series="snaive") +
  autolayer(snaive_mod, PI=FALSE, series="naive")+xlab("Time") + ylab("PPI") +
  ggtitle("Benchmark forecast") +
  guides(colour=guide_legend(title="Forecast"))


accuracy(mean_mod,myts.test)
accuracy(naive_mod,myts.test)
accuracy(snaive_mod,myts.test)

########################
#forecast arima here
########################
#myts.train<-window(myts,start=c(2003,1),end=c(2018,10))
#myts.test<-window(myts,start=c(2018,11))

fit_arima=auto.arima(myts.train, stepwise = F, approximation = F, trace = TRUE )
#my_arima<-Arima(myts.train, order = c(1,1,2))
arimafore <- forecast(fit_arima, h = 24)
autoplot(subset(myts, end = 530))+
  autolayer(arimafore, PI=FALSE)+ylab("PPI")+ggtitle("ARIMA(1,1,2) (Figure 7) ")
plot(arimafore)
checkresiduals(arimafore)
accuracy(arimafore,myts.test)



######
tso_res <- tso(myts)
tso_res
plot(tso_res)

########################
#length of our time series
########################
n <- length(myts)
tso_res$outliers
outliers_idx <- tso_res$outliers$ind
tso_res$outliers$time
mo_tc <- outliers("TC",outliers_idx)
mo_ls <- outliers("LS",outliers_idx)
mo_ao <- outliers("AO",outliers_idx)

tc <- outliers.effects(mo_tc, n)
ls <- outliers.effects(mo_ls, n)
ao <- outliers.effects(mo_ao, n)


#############
data.train_adj <- tso_res$yadj
plot(data.train_adj)

########
Y<-auto.arima(myts, d=1, D=1, stepwise = FALSE, approximation = FALSE,trace=TRUE)
Y
checkresiduals(Y)
fff<-forecast(Y,h=50)
autoplot(fff)
autoplot(subset(myts, end = 530))+ggtitle("ARIMA(0,1,2)(2,1,0)[12] (Figure 9)")+autolayer(fff, PI=FALSE)+ ylab("PPI")
accuracy(myts)

### aic=3534.148

model_auto2 <- auto.arima(data.train_adj, trace=TRUE,  ic="aic", stepwise = FALSE,seasonal = FALSE, xreg= fourier(data.train_adj, K = 6))
model_auto2
autoplot(forecast(model_auto2, xreg=fourier(data.train_adj, K=6, h=50)))+ggtitle("ARIMA 0,1,1 (Figure 6)")
accuracy(model_auto2)
# just fourier without outliers
harmonics <- fourier(myts.train, K = 6)
fit <- auto.arima(myts.train, xreg = harmonics, seasonal = FALSE)
newharmonics <- fourier(myts, K = 6, h = 50)
fc <- forecast(fit, xreg = newharmonics)
fit
accuracy(fc, myts)
autoplot(fc)
autoplot(subset(myts, end = 530))+ autolayer(fc, PI=FALSE)+ggtitle("ARIMA 1,1,1 (Figure 8)")+ylab("PPI")

#other forecast

# now check= accuracy(fc,myts.test)
accuracy(fit_arima)

###############################
#tsCV
###############################
fun_arima <- function(x, h) {
  forecast(auto.arima(x), h=h)
}

# Compute CV errors for ARIMA as e2
e <- tsCV(myts, fun_arima, h=24)
# Find MSE of each model class
mean(e^2, na.rm=TRUE)



myts %>% ets() %>% forecast() %>% autoplot()
##############################
# TSCV
##############################
gg <- tsCV(myts, fun_arima, h=1)
gg_ets <-tsCV(myts,fun_ets, h=1)


mean(gg^2, na.rm=TRUE)
mean(gg_ets^2, na.rm=TRUE)

autoplot(gg)
# how to check for residuals?
###beer <- window(ausbeer, start=1992)
#fc <- snaive(beer)
#autoplot(fc)
# <- residuals(fc)
#autoplot(res)
###############
#fit <- tslm(Demand ~ Temperature, data=daily20)
#checkresiduals(fit)
#forecast(fit, newdata=data.frame(Temperature=c(15,35)))
###
#market integration tetrad, var or ecm
### rats reg analysis for ts, cats coin anl for ts
### use PPI 6 column
### excel to text file in tetrad
### var model 2 var, vecm
### 

almond <- read_csv("Almond.csv")
na.omit(almond)
as.numeric(almond$PPI)
almond1=ts(almond[,-1], start=c(1991,12), frequency = 12)
autoplot(almond2)+ ylab("PPI")+ggtitle("Almond (Figure 1)")

adf.test(almond1)
#data is non-stationary
almond1 %>% ur.kpss() %>% summary()
# we need diffrerence (unit root test)
ndiffs(almond1)
# one diff is required
nsdiffs(almond1)
#non seasonal data

#mstl
fit_mstl<-mstl(almond1)
autoplot(fit_mstl)+ggtitle("mstl decomposition")


##peanut

Peanut <- read_csv("Peanut.csv")
Peanut<-na.omit(Peanut)
na.omit(Peanut)
as.numeric(Peanut$ppi)
#as.Date(train$date, "%m/%d/%Y")

#myts=ts(train[,-1], start=c(2003,1), frequency = 12)
Peanut2=ts(Peanut[,-1], start=c(1947,1), frequency = 12)
autoplot(Peanut2)+ ylab("PPI")+ggtitle("Peanut (Figure 2)")

adf.test(Peanut2)
#data is non-stationary
Peanut2 %>% ur.kpss() %>% summary()
# we need diffrerence (unit root test)
ndiffs(Peanut2)
# one diff is required
nsdiffs(Peanut2)
#non seasonal data
dmyts<-diff(Peanut2)

fit_mstl<-mstl(Peanut2)
autoplot(fit_mstl)+ggtitle("mstl decomposition")

########
########
########
Pecans <- read_csv("Pecans.csv")
Pecans<-na.omit(Pecans)
na.omit(Pecans)
as.numeric(Pecans$ppi)
#as.Date(train$date, "%m/%d/%Y")

#myts=ts(train[,-1], start=c(2003,1), frequency = 12)
Pecans2=ts(Pecans[,-1], start=c(1999,11), frequency = 12)
autoplot(Pecans2)+ ylab("PPI")+ggtitle("Pecans (Figure 3)")

fit_mstl1<-mstl(Pecans2)
autoplot(fit_mstl1)+ggtitle("mstl decomposition")

adf.test(Pecans2)
#data is non-stationary
Pecans2 %>% ur.kpss() %>% summary()
# we need diffrerence (unit root test)
ndiffs(Pecans2)
# one diff is required
nsdiffs(Pecans2)
#non seasonal data
dmyts<-diff(Pecans2)

#######
######
WALNUT <- read_csv("WALNUT.csv")
WALNUT<-na.omit(WALNUT)
na.omit(WALNUT)
as.numeric(WALNUT$ppi)
#as.Date(train$date, "%m/%d/%Y")

#myts=ts(train[,-1], start=c(2003,1), frequency = 12)
WALNUT2=ts(WALNUT[,-1], start=c(1998,11), frequency = 12)
autoplot(WALNUT2)+ ylab("PPI")+ggtitle("Walnut (Figure 4)")

adf.test(WALNUT2)
#data is non-stationary
WALNUT2 %>% ur.kpss() %>% summary()
# we need diffrerence (unit root test)
ndiffs(WALNUT2)
# one diff is required
nsdiffs(WALNUT2)
#non seasonal data
dmyts<-diff(WALNUT2)




