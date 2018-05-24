library(nlme)
library(ggplot2)
library(gam)
set.seed(10)
prostate = get(load("prostate2.Rdata"))
prostate$svi = as.factor(prostate$svi)
attach(prostate)
options(digits=3)
sample = sample.int(n = nrow(prostate), size = floor(.5*nrow(prostate)), replace = F)
train = prostate[sample,]
test = prostate[-sample,]

# polynomial regression

poly_list = list()
MSE_list = list()
for (i in 1:9){
  model = lm(Cscore~poly(lpsa,i), data=prostate, subset=sample)
  poly_list = c(poly_list, list(model))
  mean = mean((Cscore-predict(model,prostate))[-sample]^2)
  MSE_list = c(MSE_list, list(mean))
}
best_poly_model = poly_list[[which.min(MSE_list)]]

fit = lm(Cscore~poly(lpsa, 3), data = prostate)

lpsalims = range(lpsa)
lpsa.grid = seq(lpsalims[1], lpsalims[2])
prediction = predict(fit, newdata = list(lpsa = lpsa.grid), se = T)
SE.bands = cbind(prediction$fit+2*prediction$se.fit, prediction$fit-2*prediction$se.fit)

plot(lpsa, Cscore, xlim=lpsalims)
title("Third degree polynomial", outer = T)
lines(lpsa.grid, prediction$fit)
matlines(lpsa.grid, SE.bands, lwd=1, col="red", lty=2)

# cubic splines

# Optimising the GAM using backwards selection

# first we can check for concurvity - whether some of our full GAM model terms can be approximated
# by the others

gam.all = gam(Cscore~s(lcavol,5)+s(lweight,5)+s(age,5)+s(lbph,5)+s(lcp,5)+s(lpsa,5), data = train)

gam.nonrest = gam(Cscore~s(lcavol)+s(lweight)+s(age)+s(lbph)+s(lcp)+s(lpsa)+svi, data = train)
