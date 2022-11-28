rm(list=ls())


library(ncf)
set.seed(1234)

### Simulate full data at 100 sites over 2 years
x <- expand.grid(1:20, 1:5)[, 1]
y <- expand.grid(1:20, 1:5)[, 2]
Z <- cbind(
  rmvn.spa(x=x, y=y, p=2, method="gaus"), # true range (p) is 2 under Gaussian covariance function
  rmvn.spa(x=x, y=y, p=2, method="gaus")
)
# Fit ncf
fit0 <- Sncf(x=x, y=y, z=Z, resamp=1E3)
plot(fit0, main='Full Data')
summary(fit0)

### 30% missingness
Z.m30 <- Z
Z.m30[sample(length(x), 0.3*length(x))] <- NA
for (i in 1:nrow(Z)) {
  idx.tmp <- which(is.na(Z.m30[i, ]))
  if (length(idx.tmp>0)) {
    Z.m30[i, idx.tmp] <- mean(Z.m30[i, ], na.rm=TRUE) # replace missing data with row mean
  }
}
fit1 <- Sncf(x=x, y=y, z=Z.m30, resamp=1E3)
plot(fit1, main='30% Missing')

### 70% missingness
Z.m70 <- Z
Z.m70[sample(length(x), 0.7*length(x))] <- NA
for (i in 1:nrow(Z)) {
  idx.tmp <- which(is.na(Z.m70[i, ]))
  if (length(idx.tmp>0)) {
    Z.m70[i, idx.tmp] <- mean(Z.m70[i, ], na.rm=TRUE)
  }
}
# Fit ncf
fit2 <- Sncf(x=x, y=y, z=Z.m70, resamp=1E3)
plot(fit2, main='70% Missing')

### 50% missingness
Z.m50 <- Z
Z.m50[sample(length(x), 0.5*length(x))] <- NA
for (i in 1:nrow(Z)) {
  idx.tmp <- which(is.na(Z.m50[i, ]))
  if (length(idx.tmp>0)) {
    Z.m50[i, idx.tmp] <- mean(Z.m50[i, ], na.rm=TRUE)
  }
}
# Fit ncf
fit3 <- Sncf(x=x, y=y, z=Z.m50, resamp=1E3)
plot(fit3, main='50% Missing')

### 47% missingness
Z.m47 <- Z
Z.m47[sample(length(x), 0.47*length(x))] <- NA
for (i in 1:nrow(Z)) {
  idx.tmp <- which(is.na(Z.m47[i, ]))
  if (length(idx.tmp>0)) {
    Z.m47[i, idx.tmp] <- mean(Z.m47[i, ], na.rm=TRUE)
  }
}
# Fit ncf
fit4 <- Sncf(x=x, y=y, z=Z.m47, resamp=1E3)
plot(fit4, main='47% Missing')
summary(fit4)
