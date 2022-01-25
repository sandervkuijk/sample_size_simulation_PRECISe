################################################################################
## Simulatie power PRECISE studie
## Sander van Kuijk
##
## Start: 21/01/2022
## Laatste aanpassing: 25/01/2022
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
n     <- 392*2
mu    <- c(0.77, 0.26, 0.30, 0.35)
sd    <- c(0.23, 0.31, 0.36, 0.38)
es    <- 0.06

#     b      30     90     180
# b   0.23^2 0.18   0.15   0.21
# 30  0.18   0.31^2 0.71   0.71
# 90  0.15   0.71   0.36^2 0.89
# 180 0.21   0.71   0.89   0.38^2

# Let op! covariance = correlation*sd1*sd2
sig   <- matrix(c(sd[1]^2, 0.18*0.23*0.31, 0.15*0.23*0.36, 0.21*0.23*0.38,
                  0.18*0.31*0.23, sd[2]^2, 0.71*0.31*0.36, 0.71*0.31*0.38,
                  0.15*0.36*0.23, 0.71*0.36*0.31, sd[3]^2, 0.89*0.36*0.38,
                  0.21*0.38*0.23, 0.71*0.38*0.31, 0.89*0.38*0.36, sd[4]^2),
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
mean(mean30)    # 0.26
mean(mean90)    # 0.30
mean(mean180)   # 0.35
mean(sd30)      # 0.31
mean(sd90)      # 0.36
mean(sd180)     # 0.38
mean(c3090)     # 0.71
mean(c30180)    # 0.71
mean(c90180)    # 0.89

# Resultaten
noquote(paste("Power to detect a difference of ", es, " points is ",
               (sum(pval <= 0.05)/iters)*100, "%.", sep = ""))
round(mean(beta), 3)

### Einde file.