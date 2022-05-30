################################################################################
## Simulatie power PRECISE studie
## Sander van Kuijk
##
## Start: 21/01/2022
## Laatste aanpassing: 30/05/2022
################################################################################

library(MASS)
library(reshape2)
library(nlme)
library(matrixcalc)
library(Matrix)
library(lqmm)

# Simulatieparameters
iters <- 1000

# Steekproefkenmerken
n     <- 320*2 # power = 79.8%
mu    <- c(0.79, 0.48, 0.60, 0.67)
sd    <- c(0.23, 0.28, 0.28, 0.25)
es    <- 0.06

# # Let op! covariance = correlation*sd1*sd2
sig   <- matrix(c(sd[1]^2, 0.16*0.23*0.28, 0.13*0.23*0.28, 0.15*0.23*0.25,
                 0.16*0.28*0.23, sd[2]^2, 0.48*0.28*0.28, 0.41*0.28*0.25,
                 0.13*0.28*0.23, 0.48*0.28*0.28, sd[3]^2, 0.68*0.28*0.25,
                 0.15*0.25*0.23, 0.41*0.25*0.28, 0.68*0.25*0.28, sd[4]^2),
                 nrow = 4, byrow = TRUE)
is.positive.definite(sig)

# Dubbelcheck simulatie data
mean30  <- rep(NA, iters)
mean90  <- rep(NA, iters)
mean180 <- rep(NA, iters)
sd30    <- rep(NA, iters)
sd90    <- rep(NA, iters)
sd180   <- rep(NA, iters)
c3090   <- rep(NA, iters)
c30180  <- rep(NA, iters)
c90180  <- rep(NA, iters)

# Resultaten opslaan
beta    <- rep(NA, iters)
pval    <- rep(NA, iters)
  
# Reproduceerbaarheid
set.seed(7181)

# Simulatie
for(i in 1:iters){
  
  Xi <- mvrnorm(n, mu = mu, Sigma = sig, empirical = TRUE)
  Xi[, 1] <- Xi[, 1]
  Xi[, 2] <- Xi[, 2]
  Xi[, 3] <- Xi[, 3]
  Xi[, 4] <- Xi[, 4]
  
  d  <- data.frame(ID = seq(1:n), Xi)
  dl <- melt(d, id.vars = c("ID", "X1"), measure.vars = c("X2", "X3", "X4"),
           variable.name = "fu_moment", value.name = "qol")
  names(dl)[2] <- "baseline"
  dl <- dl[order(dl$ID), ]
  dl$fu_t <- as.numeric(rep(c(30, 90, 180), n))
  dl$rand <- c(rep(0, 0.5*length(dl$ID)), rep(1, 0.5*length(dl$ID)))
  dl$qol  <- ifelse(dl$rand == 1, dl$qol + 0.5*es, dl$qol - 0.5*es)
  
  mean30[i]  <- mean(dl$qol[dl$fu_t == 30])
  mean90[i]  <- mean(dl$qol[dl$fu_t == 90])
  mean180[i] <- mean(dl$qol[dl$fu_t == 180])
  sd30[i]    <- sd(dl$qol[dl$fu_t == 30])
  sd90[i]    <- sd(dl$qol[dl$fu_t == 90])
  sd180[i]   <- sd(dl$qol[dl$fu_t == 180])
  c3090[i]   <- cor(dl$qol[dl$fu_t == 30], dl$qol[dl$fu_t == 90])
  c30180[i]  <- cor(dl$qol[dl$fu_t == 30], dl$qol[dl$fu_t == 180])
  c90180[i]  <- cor(dl$qol[dl$fu_t == 90], dl$qol[dl$fu_t == 180])

  fit <- lme(qol ~ rand + baseline, data = dl, rand = ~ 1|ID,
             correlation = corAR1(form = ~ fu_t | ID),
             control = list(opt = "optim"))

  beta[i] <- fixef(fit)[2]
  pval[i] <- summary(fit)$tTable[2, 5]
  
}

# Simulatie check
mean(mean30)
mean(mean90)
mean(mean180)
mean(sd30)
mean(sd90)
mean(sd180)
mean(c3090)
mean(c30180)
mean(c90180)

# Resultaten
noquote(paste("Power to detect a difference of ", es, " points is ",
               (sum(pval <= 0.05)/iters)*100, "%.", sep = ""))
round(mean(beta), 3)

### Einde file.