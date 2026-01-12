seed_id = 3000
set.seed(seed_id)

#install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue, tidyverse, coda)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

setwd("/work/CarlEmilGrum-Nymann#2429/EXAM")

citation("polspline")

##### Preprocessing #####

groupSize <- 4
ntrials <- 10
pi <- 1.6 # used to be 1.4, but the original paper (and Josh' preprint both say 1.6)
ntokens <- 20
vals <- seq(0,ntokens,1)
#vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

rawDat <- read.csv("data/HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game

#- create covariates in raw data matrix
# nation index
rawDat$nation <- c()
rawDat$nation[rawDat$city=="Melbourne"]=1
rawDat$nation[rawDat$city=="Minsk"]=2
rawDat$nation[rawDat$city=="Chengdu"]=3
rawDat$nation[rawDat$city=="Copenhagen"]=4
rawDat$nation[rawDat$city=="Bonn"]=5
rawDat$nation[rawDat$city=="Athens"]=6
rawDat$nation[rawDat$city=="Seoul"]=7
rawDat$nation[rawDat$city=="Samara"]=8
rawDat$nation[rawDat$city=="Zurich"]=9
rawDat$nation[rawDat$city=="St. Gallen"]=9
rawDat$nation[rawDat$city=="Istanbul"]=10
rawDat$nation[rawDat$city=="Nottingham"]=11
rawDat$nation[rawDat$city=="Dnipropetrovs'k"]=12
rawDat$nation[rawDat$city=="Boston"]=13

# create variable for EnvCon. Old data from 
# http://hdr.undp.org/sites/default/files/reports/269/hdr_2009_en_complete.pdf,
# see the suppl material in the published paper for the updated values:
# https://www.sciencedirect.com/science/article/pii/S2666622723000254#sec0014

# # Old pre-published EnvCon-values
# rawDat$EnvCon <- c()
# rawDat$EnvCon[rawDat$city=="Melbourne"]=34.3
# rawDat$EnvCon[rawDat$city=="Minsk"]=25.3
# rawDat$EnvCon[rawDat$city=="Chengdu"]=38.5
# rawDat$EnvCon[rawDat$city=="Copenhagen"]=28.7
# rawDat$EnvCon[rawDat$city=="Bonn"]=31.9
# rawDat$EnvCon[rawDat$city=="Athens"]=34.4
# rawDat$EnvCon[rawDat$city=="Seoul"]=31.6
# rawDat$EnvCon[rawDat$city=="Samara"]=37.5
# rawDat$EnvCon[rawDat$city=="Zurich"]=32.7
# rawDat$EnvCon[rawDat$city=="St. Gallen"]=32.7
# rawDat$EnvCon[rawDat$city=="Istanbul"]=41.9
# rawDat$EnvCon[rawDat$city=="Nottingham"]=34.8
# rawDat$EnvCon[rawDat$city=="Dnipropetrovs'k"]=26.1
# rawDat$EnvCon[rawDat$city=="Boston"]=41.1

# EnvCon-values taken directly from the supplementary materials in the published paper
rawDat$EnvCon <- c()
rawDat$EnvCon[rawDat$city=="Melbourne"]= 55.5
rawDat$EnvCon[rawDat$city=="Minsk"]=61.5 # no data wv5 - from EVS
rawDat$EnvCon[rawDat$city=="Chengdu"]= 73
rawDat$EnvCon[rawDat$city=="Copenhagen"]=70 # no data wv5 - from EVS
rawDat$EnvCon[rawDat$city=="Bonn"]=35.1
rawDat$EnvCon[rawDat$city=="Athens"]=81 # no data wv5 - from EVS
rawDat$EnvCon[rawDat$city=="Seoul"]=75.8
rawDat$EnvCon[rawDat$city=="Samara"]=56.1 # not asked about this specific question - from evs
rawDat$EnvCon[rawDat$city=="Zurich"]=62.6
rawDat$EnvCon[rawDat$city=="St. Gallen"]=62.6
rawDat$EnvCon[rawDat$city=="Istanbul"]=77.7
rawDat$EnvCon[rawDat$city=="Nottingham"]=42.7 # not asked this specific question - from EVS
rawDat$EnvCon[rawDat$city=="Dnipropetrovs'k"]=41
rawDat$EnvCon[rawDat$city=="Boston"]=50.8

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_no_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

missing <- array(0,ngroups)

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Gga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[s,,g] <- colSums(c_no_punish[-s,,g])
    Ga_no_punish[s,,g] <- colMeans(c_no_punish[-s,,g])
    Ggas_no_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# data for punishment condition #
c_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Gga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    Gc_punish[s,,g] <- colSums(c_punish[-s,,g])
    Ga_punish[s,,g] <- colMeans(c_punish[-s,,g])
    Ggas_punish[s,,g] <- colMeans(c_punish[,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Gga <- array(0,c(ntrials,ngroups,2))
Gga[,,1] <- Gga_no_punish
Gga[,,2] <- Gga_punish

Ggas <- array(0,c(groupSize,ntrials,ngroups,2))
Ggas[,,,1] <- Ggas_no_punish
Ggas[,,,2] <- Ggas_punish


Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

Ga <- array(0,c(groupSize,ntrials,ngroups,2))
Ga[,,,1] <- Ga_no_punish
Ga[,,,2] <- Ga_punish

c_choice_index <- c

EnvCon <- array(0, ngroups)
Nation <- array(0, ngroups)
for (g in 1:ngroups) {
  EnvCon[g] <- mean(redDat$EnvCon[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
  Nation[g] <- mean(redDat$nation[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
}

c_win <- c_no_punish[,,!is.na(EnvCon)]
c_keep <- rep()-c_win

Gga_punish <- Gga_punish[,!is.na(EnvCon)]
Gga_no_punish <- Gga_no_punish[,!is.na(EnvCon)]

c <- c[,,!is.na(EnvCon),]
Gga <- Gga[,!is.na(EnvCon),]
Ggas <- Ggas[,,!is.na(EnvCon),]
Gc <- Gc[,,!is.na(EnvCon),]
Ga <- Ga[,,!is.na(EnvCon),]
EnvCon <- EnvCon[!is.na(EnvCon)]
Nation <- Nation[!is.na(Nation)]

#redefine number of groups after removing those without civic scores
ngroups <- length(EnvCon)

# aggregate EnvCon to just 1 number per Nation-index (using the mean here should be unproblematic since all groups within a given nation should have been given the same EnvCon-coefficient)
EnvCon <- aggregate(EnvCon~Nation, FUN=mean)[,2]

nnations <- length(EnvCon)

# calculate the winnings (i.e. apply the multiplication-factor to the sum of each groups contributions)
winnings <-  array(0, ngroups)
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c_win[,,g])*pi)
}

################################################################################
########################### Conditional cooperation model ######################
################################################################################

# JZS priors for partial correlation. Method described here
# https://link.springer.com/article/10.3758/s13423-012-0295-x
# Code available here
# https://github.com/MicheleNuijten/BayesMed/blob/master/R/jzs_corSD.R
# Paper where code is used here (mediation paper)
# https://link.springer.com/article/10.3758/s13428-014-0470-2

#################################################################
#------------------ Winnings analysis ---------------------------
#################################################################

#-------------------  Regress EnvCon on winnings ---------------

# standardise variables

X <- EnvCon
X <- (X-mean(X))/sd(X)

invSigma <- solve(t(X)%*%X) # required for JZS priors

Y <- (winnings-mean(winnings))/sd(winnings)

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX") 

# - run jags code
win.samples <- jags.parallel(data, inits=NULL, params,
                             model.file ="win_corr.txt",
                             n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)

print(win.samples)

# Look at convergence diagnostics
png(filename="winnings_traces.png",
    height = 1400,
    width = 1000)

layout(rbind(c(1,1,1,1,4,2,2,2,2), 
             c(1,1,1,1,4,2,2,2,2),
             c(1,1,1,1,4,2,2,2,2), 
             c(3,3,3,3,4,4,4,4,5),
             
             c(3,3,3,3,7,7,9,9,9), 
             c(3,3,3,3,7,7,9,9,9),
             c(6,6,6,8,8,8,10,10,10), 
             c(6,6,6,8,8,8,10,10,10)))

# Look at convergence diagnostics
win_samples <- as.mcmc(win.samples)

traceplot(win_samples)

dev.off()
#################################################################
#------------------ CC model analysis ---------------------------
#################################################################

#-------------------  Regress EnvCon on belief weights and slope of prefs in CC model ---------------

# standardise covariate
X <- EnvCon
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors
# Ga_old <- Ga
# Ga <- Ggas
data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") 

# - run jags code - original model from skewes and Nockur
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                                  model.file ="CC_corr.txt",
                                  n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)
end_time = Sys.time()
end_time - start_time

png(filename="CC_traces.png",
    height = 1400,
    width = 1000)

layout(rbind(c(1,1,1,2,2,2,3,3,3), 
             c(1,1,1,2,2,2,3,3,3),
             c(4,4,4,5,5,5,6,6,6), 
             c(4,4,4,5,5,5,6,6,6),
             
             c(12,12,12,7,7,7,9,9,9), 
             c(12,12,12,7,7,7,9,9,9),
             c(11,11,11,8,8,8,10,10,10), 
             c(11,11,11,8,8,8,10,10,10)))

# Look at convergence diagnostics
print(CC.samples)

CC_samples <- as.mcmc(CC.samples)

traceplot(CC_samples)

dev.off()


#run means model without correlation
data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("mu_alpha","mu_rho", "mu_omega", "prec_alpha","prec_rho", "prec_omega") 


# - run jags code
start_time = Sys.time()
CC_means.samples <- jags.parallel(data, inits=NULL, params,
                            model.file ="CC_estimating_means_simple_no_hypermu.txt",
                            n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)
end_time = Sys.time()
end_time - start_time

library(coda)

max_plots <- 15
param_names <- varnames(CC_means_samples)

param_groups <- split(
  param_names,
  ceiling(seq_along(param_names) / max_plots)
)

for (g in seq_along(param_groups)) {
  
  n_panels <- length(param_groups[[g]])
  n_row <- ceiling(n_panels / 3)
  n_col <- 3
  
  png(
    filename = paste0("traceplot_page_", g, ".png"),
    width = 1800,
    height = 350 * n_row,
    res = 200
  )
  
  par(mfrow = c(n_row, n_col),
      mar = c(2, 2, 2, 1),
      oma = c(0, 0, 3, 0))
  
  for (p in param_groups[[g]]) {
    traceplot(
      CC_means_samples[, p],
      main = p
    )
  }
  
  mtext(
    paste("Traceplots:", 
          ((g - 1) * max_plots + 1),
          "–",
          min(g * max_plots, length(param_names))),
    outer = TRUE,
    line = 1,
    cex = 1.2
  )
  
  dev.off()
}


# Look at convergence diagnostics
png(filename="CC_means_traces.png",
    height = 1400,
    width = 1400)

layout(rbind(c(1,1,2,2,3,3,4,4), 
             c(1,1,2,2,3,3,4,4),
             c(5,5,6,6,7,7,8,8), 
             c(5,5,6,6,7,7,8,8),
            
             c(9,9,10,10,7,7,9,9), 
             c(9,9,10,7,7,7,9,9),
             c(11,11,11,8,8,8,10,10), 
             c(11,11,11,8,8,8,10,10)))



#print(CC_means.samples)

CC_means_samples <- as.mcmc(CC_means.samples)

traceplot(CC_means_samples)

dev.off()

#################################################################
#------------------ Plotting ---------------------------
#################################################################

# ------ Create empirical and parameter arrays for plots -------
# empirical group winnings data - means and standard deviations
empirical.win <- array(0,c(3,length(EnvCon)))
for (i in 1:length(EnvCon)) {
  empirical.win[1,i] <- mean(winnings[Nation==i]) - sd(winnings[Nation==i]) 
  empirical.win[2,i] <- mean(winnings[Nation==i]) 
  empirical.win[3,i] <- mean(winnings[Nation==i]) + sd(winnings[Nation==i])
}

empirical.initial <- array(0,c(3,length(EnvCon)))
for (i in 1:length(EnvCon)) {
  empirical.initial[1,i] <- mean(colMeans(c[,1,,1])[Nation==i]) - sd(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[2,i] <- mean(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[3,i] <- mean(colMeans(c[,1,,1])[Nation==i]) + sd(colMeans(c[,1,,1])[Nation==i])
}

empirical.contrib <- array(0,c(3,length(EnvCon)))
for (i in 1:length(EnvCon)) {
  empirical.contrib[1,i] <- mean(colSums(colSums(c_win[,,Nation==i]))) - sd(colSums(colSums(c_win[,,Nation==i]))) 
  empirical.contrib[2,i] <- mean(colSums(colSums(c_win[,,Nation==i]))) 
  empirical.contrib[3,i] <- mean(colSums(colSums(c_win[,,Nation==i]))) + sd(colSums(colSums(c_win[,,Nation==i])))
}

#-------------- create high res png for saving ---------------------------------

png(filename="PGG_figure.png",
    height = 1400,
    width = 1000)

# set layout for graphic
layout(rbind(c(1,1,1,2,2,2,4,4,4), 
             c(1,1,1,2,2,2,4,4,4),
             c(1,1,1,3,3,3,4,4,4), 
             c(1,1,1,3,3,3,4,4,4),
             
             c(5,5,5,7,7,7,9,9,9), 
             c(5,5,5,7,7,7,9,9,9),
             c(6,6,6,8,8,8,10,10,10), 
             c(6,6,6,8,8,8,10,10,10)))

# text scale plot parameter
par(cex=1.5)

# set margins - will change when there is no title required
par(mar=c(5,3,5,2))


#----- empirical correlation - EnvCon and group winnings --------------------
plot(c(22,100), c(0,900), type = "n", main = "A: Data - Winnings", 
     xlab = "National EnvCon Coefficient", ylab = "Group winnings",axes=FALSE)
for (i in 1:length(EnvCon)) {
  lines(c(EnvCon[i],EnvCon[i]),c(empirical.win[1,i],empirical.win[3,i]))
  points(EnvCon[i],empirical.win[2,i])
}
axis(1)
axis(2)

#----- posterior of effect - EnvCon and group winnings --------------------
plot(density(win.samples$BUGSoutput$sims.list$betaX),frame=FALSE,lwd=2,ylim=c(0,7),
     cex=2,xlab = "Standardised Effect of EnvCon",ylab = " ",main="B: Winnings")
COL <- adjustcolor(c("red"))
lines(c(win.samples$BUGSoutput$summary[2,3],win.samples$BUGSoutput$summary[2,7]),c(-.1,-.1),col=COL,lwd=2)
points(win.samples$BUGSoutput$summary[2,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2)) # reset margin because no title


#plots based on the correlation model

#----- empirical correlation - EnvCon and initial contribution --------------------
plot(c(22,100), c(0,20), type = "n", main = "C: Data - Initial Contrib.", 
     xlab = "National EnvCon Coefficient", ylab = "Initial Contribution",axes=FALSE)
for (i in 1:length(EnvCon)) {
  lines(c(EnvCon[i],EnvCon[i]),c(empirical.initial[1,i],empirical.initial[3,i]))
  points(EnvCon[i],empirical.initial[2,i])
}
axis(1)
axis(2)

par(mar=c(5,3,1,2)) # reset margin because no title

#---------- Initial Belief -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_alpha),frame=FALSE,lwd=2,ylim=c(0,8),
     cex=2,xlab = "Standardised Effect of EnvCon",ylab = " ",main="D: Initial Belief")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[4,3],CC.samples$BUGSoutput$summary[4,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[4,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2)) # reset margin because no title



#---------- Belief Learning -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_omega),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardised Effect of EnvCon",ylab = " ",main="E: Belief Learning Weight")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[5,3],CC.samples$BUGSoutput$summary[5,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[5,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))


#---------- Preference Slope -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_rho),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardised Effect of EnvCon",ylab = " ",main="F: Conditional Preferences")
COL <- adjustcolor(c("black"))
lines(c(CC.samples$BUGSoutput$summary[6,3],CC.samples$BUGSoutput$summary[6,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[6,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))


#----- empirical correlation - EnvCon and group contributions --------------------
plot(c(22,100), c(0,900), type = "n", main = "H: Data - Contributions", 
     xlab = "National EnvCon Coefficient", ylab = "Group contributions",axes=FALSE)
for (i in 1:length(EnvCon)) {
  lines(c(EnvCon[i],EnvCon[i]),c(empirical.contrib[1,i],empirical.contrib[3,i]))
  points(EnvCon[i],empirical.contrib[2,i])
}
axis(1)
axis(2)


par(mar=c(5,3,1,2)) # reset margin because no title

dev.off()


png(filename="Means_estimation_figure.png",
    height = 1400,
    width = 1000)

# set layout for graphic
layout(rbind(c(1,1,1,1,4,2,2,2,2), 
             c(1,1,1,1,4,2,2,2,2),
             c(1,1,1,1,4,2,2,2,2), 
             c(3,3,3,3,4,4,4,4,5),
             
             c(3,3,3,3,7,7,9,9,9), 
             c(3,3,3,3,7,7,9,9,9),
             c(6,6,6,8,8,8,10,10,10), 
             c(6,6,6,8,8,8,10,10,10)))



#plotting based on the means model

#----- Distribution of estimated national alpha means --------------------
alpha_rows <- grep("^mu_alpha", rownames(CC_means.samples$BUGSoutput$summary))
nnations   <- length(alpha_rows)
plot(
  x = c(
    min(CC_means.samples$BUGSoutput$summary[alpha_rows, 3]),
    max(CC_means.samples$BUGSoutput$summary[alpha_rows, 7])
  ),
  y = c(0.5, nnations + 0.5),
  type = "n",
  xlab = expression("National mean " * alpha),
  ylab = "Nation",
  main = "A: Mational initial belief means",
  cex.main = 2,
  axes = FALSE
)

COL <- adjustcolor("black", alpha.f = 0.8)

for (i in seq_along(alpha_rows)) {
  r <- alpha_rows[i]
  
  lines(
    c(
      CC_means.samples$BUGSoutput$summary[r, 3],  # lower
      CC_means.samples$BUGSoutput$summary[r, 7]   # upper
    ),
    c(i, i),
    lwd = 2,
    col = COL
  )
  
  points(
    CC_means.samples$BUGSoutput$summary[r, 1],   # mean
    i,
    pch = 19,
    col = COL
  )
}

axis(1)
axis(2, at = 1:nnations, labels = 1:nnations, las = 1)
abline(v = 0, col = "gray70", lty = 2)




#----- Distribution of estimated national learning rate means --------------------
omega_rows <- grep("^mu_omega", rownames(CC_means.samples$BUGSoutput$summary))
plot(
  x = c(
    min(CC_means.samples$BUGSoutput$summary[omega_rows, 3]),
    max(CC_means.samples$BUGSoutput$summary[omega_rows, 7])
  ),
  y = c(0.5, nnations + 0.5),
  type = "n",
  xlab = expression("National mean " * omega),
  ylab = "Nation",
  main = "B: National mean sensitivity to others",
  cex.main = 2,
  axes = FALSE
)

COL <- adjustcolor("black", alpha.f = 0.8)

for (i in seq_along(omega_rows)) {
  r <- omega_rows[i]
  
  lines(
    c(
      CC_means.samples$BUGSoutput$summary[r, 3],  # lower
      CC_means.samples$BUGSoutput$summary[r, 7]   # upper
    ),
    c(i, i),
    lwd = 2,
    col = COL
  )
  
  points(
    CC_means.samples$BUGSoutput$summary[r, 1],   # mean
    i,
    pch = 19,
    col = COL
  )
}

axis(1)
axis(2, at = 1:nnations, labels = 1:nnations, las = 1)
abline(v = 0, col = "gray70", lty = 2)

#----- Distribution of estimated national rho means --------------------
rho_rows <- grep("^mu_rho", rownames(CC_means.samples$BUGSoutput$summary))
plot(
  x = c(
    min(CC_means.samples$BUGSoutput$summary[rho_rows, 3]),
    max(CC_means.samples$BUGSoutput$summary[rho_rows, 7])
  ),
  y = c(0.5, nnations + 0.5),
  type = "n",
  xlab = expression("National mean " * rho),
  ylab = "Nation",
  main = "C: National cooperation preference means",
  cex.main = 2,
  axes = FALSE
)

COL <- adjustcolor("black", alpha.f = 0.8)

for (i in seq_along(rho_rows)) {
  r <- rho_rows[i]
  
  lines(
    c(
      CC_means.samples$BUGSoutput$summary[r, 3],  # lower
      CC_means.samples$BUGSoutput$summary[r, 7]   # upper
    ),
    c(i, i),
    lwd = 2,
    col = COL
  )
  
  points(
    CC_means.samples$BUGSoutput$summary[r, 1],   # mean
    i,
    pch = 19,
    col = COL
  )
}

axis(1)
axis(2, at = 1:nnations, labels = 1:nnations, las = 1)
abline(v = 0, col = "gray70", lty = 2)

dev.off()

