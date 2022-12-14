
oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)
X = scale(X[1:400,])/sqrt(ncol(X))   # Subset and scale markers
G = tcrossprod(X)                    # Genomic relationship matrix
y = as.vector(scale(Y[1:400,"YLD"])) # Subset response variable

# Calculate heritability using all data
fm1 = fitBLUP(y,K=G)
h2 = fm1$varU/(fm1$varU + fm1$varE)

# Sparse selection index
fm2 = SSI(y,K=G,h2=h2,nLambda=50)
yHat = fitted(fm2)

plot(fm2)  # Penalization vs accuracy

# Equivalence of the SSI with lambda=0 with G-BLUP
fm3 = SSI(y,K=G,h2=h2,lambda=0,tol=1E-5)

cor(y,fm1$u)        # G-BLUP accuracy
cor(y,fitted(fm3))  # SSI accuracy

# Predicting a testing set using training set
tst = sample(seq_along(y),ceiling(0.3*length(y)))
trn = (seq_along(y))[-tst]

# Calculate heritability in training data
yNA = y
yNA[tst] = NA
fm1 = fitBLUP(yNA,K=G)
(h2 = fm1$varU/(fm1$varU + fm1$varE))

# Sparse selection index
fm2 = SSI(y,K=G,h2=h2,trn=trn,tst=tst)

# Heritability internaly calculated
fm2 = SSI(y,K=G,h2=NULL,trn=trn,tst=tst)
fm2$h2

# Effect of the penalization on the accuracy
plot(fm2)

setwd(oldwd)
