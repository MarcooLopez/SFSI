
oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)

index = which(Y$trial %in% 1:8)     # Use only a subset of data
Y = Y[index,]
M = scale(M[index,])/sqrt(ncol(M))  # Subset and scale markers
G = tcrossprod(M)                   # Genomic relationship matrix
y = scale(Y[,'E1'])            # Subset response variable

# Fit model, whole data
fm1 = fitBLUP(y,K=G)
fm1$varU
fm1$varE
fm1$theta
fm1$h2
cor(y,fm1$u)                 # Prediction accuracy

# Same model different parametrization
fm2 = fitBLUP(y,Z=M)
fm2$varU; fm2$varE; fm2$theta; fm2$h2
cor(y,M%*%fm2$u)             # Prediction accuracy

# Training and testing sets
tst = which(Y$trial %in% 2)
trn = seq_along(y)[-tst]

yNA <- y
yNA[tst] <- NA

# Fit model, split data
fm3 = fitBLUP(yNA,K=G)
plot(y[tst],fm3$u[tst])      # Predicted vs observed values in testing set
cor(y[tst],fm3$u[tst])       # Prediction accuracy in testing set
cor(y[trn],fm3$u[trn])       # Prediction accuracy in training set
fm3$h2                       # Heritability (in training set)

setwd(oldwd)
