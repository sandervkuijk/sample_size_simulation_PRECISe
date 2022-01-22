################################################################################
## Simulatie power PRECISE studie
## Sander van Kuijk
##
## Start: 21/01/2022
## Laatste aanpassing: 22/01/2022
################################################################################

library(MASS)
library(reshape2)
library(nlme)

# Simulatieparameters
iters <- 100

# Steekproefkenmerken
n     <- 600
mu    <- c(20, 30, 40) 
rho   <- c(0.6, 0.3)
sd    <- c(30, 35, 38)
es    <- 6

# Resultaten opslaan
beta <- rep(NA, iters)
pval <- rep(NA, iters)
  
# Simulatie
for(i in 1:iters){
  
  Xi <- mvrnorm(n, mu = mu, Sigma = matrix(c(sd[1]^2, rho, rho[1], sd[2]^2,
                                           rho[1], rev(rho), sd[3]^2), nrow=3))
  d  <- data.frame(ID = seq(1:n), Xi)
  dl <- melt(d, id.vars = "ID", measure.vars = c("X1", "X2", "X3"),
           variable.name = "fu_moment", value.name = "qol")
  dl <- dl[order(dl$ID), ]

  dl$fu_t <- as.numeric(rep(c(1, 3, 6), n))
  dl$rand <- c(rep(0, 0.5*length(dl$ID)), rep(1, 0.5*length(dl$ID)))
  dl$qol  <- ifelse(dl$rand == 1, dl$qol+es, dl$qol)

  fit <- lme(qol ~ rand + fu_t, data = dl, rand = ~ 1|ID,
             correlation = corAR1(form = ~ fu_t | ID))
  beta[i] <- fixef(fit)[2]
  pval[i] <- summary(fit)$tTable[2, 5]
  
}

# Resultaten
noquote(paste("Power to detect a difference of ", es, " points is ",
               (sum(pval <= 0.05)/iters)*100, "%.", sep = ""))
round(mean(beta), 1)
