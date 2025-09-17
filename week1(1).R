library(gamair)
data(hubble)
with(hubble,plot(x,y,xlab="Distance (Mpc)",ylab="Velocity (km/s)"))
hub.mod <- lm(y~x-1,data=hubble)
summary(hub.mod)
abline(a=0, b=hub.mod$coefficients)

##outliers
points(hubble$x[c(3,15)], hubble$y[c(3,15)], col="red", pch=16)

## omit possible outliers...
hub.mod1 <- lm(y~x-1,data=hubble[-c(3,15),])
summary(hub.mod1)
abline(a=0, b=hub.mod1$coefficients, col="red")

## convert Hubble const to age in years...
##1 megaparsec = 3.09e19 km
hubble.const <- c(coef(hub.mod),coef(hub.mod1))/3.09e19
age <- 1/hubble.const ## in seconds
(age/(60^2*24*365) )/1e9


## residuals...
plot(fitted(hub.mod),residuals(hub.mod),xlab="fitted values",
     ylab="residuals")

plot(fitted(hub.mod1),residuals(hub.mod1),
      xlab="fitted values",ylab="residuals")


##manually calculating standard deviation
sigmasq= 1/(length(hubble$x)-1) * sum((residuals(hub.mod))^2)
sqrt(sigmasq/sum((hubble$x)^2))
stde = sqrt(sigmasq/sum((hubble$x)^2))

##compared to:
summary(hub.mod)$coefficients


## Testing creationist hypothesis
cs.age = 6000 * (60^2*24 * 365)
cs.hubble = (1/cs.age)*3.09e19 ###roughly 163000000

##95% of the data within 2 sigma of mu
pnorm(2)-pnorm(-2)
##99% of the data within 3 sigma of mu
pnorm(3)-pnorm(-3)
## "New physics' threshold: 1 in 3.5 million
1/pnorm(-5)/1e6


round(abs(coef(hub.mod)-cs.hubble)/stde) 

# testing the Creationist hypothesis
t.stat<-abs(coef(hub.mod)-cs.hubble)/summary(hub.mod)$coefficients[2]
pt(-t.stat,df=length(hubble$x)-1)*2 

## reproducing the t and p-values for beta=0
t0 = abs(coef(hub.mod))/summary(hub.mod)$coefficients[2]
p0 = pt(-t0,df=length(hubble$x)-1)*2 
t0
p0
summary(hub.mod)$coefficients

# Confidence interval
q = qt(0.025, df = length(hubble$x)-1)
sigb <- summary(hub.mod1)$coefficients[2]
h.ci <- coef(hub.mod)+c(q,-q)*sigb
h.ci ## Hubble CI

h.ci <- h.ci*60^2*24*365.25/3.09e19 # convert to 1/years
sort(1/h.ci)/1e9 #billions of years
