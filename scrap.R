library(nlme)
library(ggplot2)
set.seed(10)
prostate = get(load("prostate2.Rdata"))
prostate$svi = as.factor(prostate$svi)
attach(prostate)
options(digits=3)
sample = sample.int(n = nrow(prostate), size = floor(.5*nrow(prostate)), replace = F)
train = prostate[sample,]
test = prostate[-sample,]
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

