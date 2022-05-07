library("cpm")

# Figures 1a) and 1b) don't use the cpm package

# Figure 2a)
set.seed(0)
x <- c(rnorm(200, 0, 1), rnorm(200, 0.5, 1))
resultsStudent <- detectChangePointBatch(x, cpmType = "Student", alpha = 0.05)
resultsMW <- detectChangePointBatch(x, cpmType = "Mann-Whitney", alpha = 0.05)
plot(x, type = "l", xlab = "Observation", ylab = "x", bty = "l")
if (resultsStudent$changeDetected) 
  abline(v = resultsStudent$changePoint, lty = 2)
if (resultsMW$changeDetected) 
  abline(v = resultsMW$changePoint, lty = 2)


# Figure 2b) (run above code first)
plot(resultsMW$Ds, type = "l", xlab = "Observation", ylab = expression(D[t]), bty = "l")
abline(h = resultsMW$threshold, lty = 2)

# Figure 3a)
set.seed(0)
x <- c(rnorm(200, 0, 1), rnorm(200, 1, 1))
resultsStudent <- detectChangePoint(x, cpmType = "Student", ARL0 = 500)
resultsMW <- detectChangePoint(x, cpmType = "Mann-Whitney", ARL0 = 500)
plot(x, type = "l", bty = "l")
if (resultsStudent$changeDetected) 
  abline(v = resultsStudent$detectionTime, col = "red")
if (resultsMW$changeDetected) 
  abline(v = resultsMW$detectionTime, col = "blue")

# Figure 3b) (run above code first)
plot(resultsMW$Ds, xlab = "Observation", ylab = expression(D[t]), type = "l", bty = "l")
lines(resultsMW$thresholds, col = 2)

# Table 1
cpmTypes <- c("Bartlett", "Mood", "Lepage")
changeMagnitudes <- c(1.5, 2.0, 3.0, 0.5, 0.3)
changeLocations <- c(50, 300)
sims <- 50000
ARL0 <- 500
startup <- 20
resultsMean <- resultsSD <- list()
for (cpmType in cpmTypes) {    
   resultsMean[[cpmType]] <- resultsSD[[cpmType]] <- matrix(numeric(length(changeMagnitudes) *
     length(changeLocations)), nrow = length(changeMagnitudes),
     dimnames = list(changeMagnitudes, changeLocations))
   for (cm in 1:length(changeMagnitudes)) {
     for (cl in 1:length(changeLocations)) {
       cat(sprintf("cpm:%s magnitude::%s location:%s\n",
         cpmType, changeMagnitudes[cm], changeLocations[cl]))
       temp <- numeric(sims)
       for (s in 1:sims) {
         x <- c(rnorm(changeLocations[cl], 0, 1), 
           rnorm(2000, 0, changeMagnitudes[cm]))
         temp[s] <- detectChangePoint(x, cpmType, 
           ARL0 = ARL0, startup = startup)$detectionTime
       }
       resultsMean[[cpmType]][cm, cl] <- 
         mean(temp[temp > changeLocations[cl]]) -
           changeLocations[cl]
       resultsSD[[cpmType]][cm, cl] <- 
         sd(temp[temp > changeLocations[cl]])
    }
  }
}
lapply(resultsMean, round, digits = 1)
lapply(resultsSD, round, digits = 1)

# Table 2
theta0 <- 0.4
theta1 <- 0.6
ARL0s <- c(100, 500, 1000)
startup <- 20
changePoint <- 100
sims <- 1000
results <- matrix(numeric(length(ARL0s) * sims), ncol = sims)
for (i in 1:length(ARL0s)) {    
 for (s in 1:sims) {
   x <- c(rbinom(changePoint, 1, theta0), rbinom(500, 1, theta1))
   results[i, s] <- detectChangePoint(x, "FET", 
     ARL0 = ARL0s[i], startup = startup, lambda = 0.1)$detectionTime
 }
}
# average delays
apply(results, 1, function(z) mean(z[z >= changePoint]) - changePoint)
# proportion of false positives
apply(results, 1, function(z) length(z[z < changePoint])/sims)
# output results as a table
library("xtable") # used to produce LaTeX markup
M <- cbind(ARL0s,
 apply(results, 1, function(z) mean(z[z >= changePoint]) - changePoint),
 apply(results, 1, function(z) length(z[z < changePoint])/sims))
colnames(M) <- c("ARL0", "Delay", "False Positives")
rownames(M) <- NULL
xtable(M, label = NULL)

# Figure 4
set.seed(0)
x <- c(rt(50, 5), rt(50, 5) + 3, rt(50, 5))
res <- processStream(x, cpmType = "Mann-Whitney", ARL0 = 500, startup = 20)

plot(x, type = "l", xlab = "Observation", ylab = "", bty = "l")
abline(v = res$detectionTimes)
abline(v = res$changePoints, lty = 2)


# Figures 5
data("ForexData", package = "cpm")
x <- ForexData[, 3]
diffs <- diff(log(x))
results <- processStream(diffs, cpmType = "Mood", ARL0 = 5000, startup = 20)

# Figure 5a)
plot(x, type = "l", xlab = "Observation", ylab = "log(Exchange Rates)", bty = "l")

# Figure 5b) (run above code first)
plot(diffs, type = "l", xlab = "Observation", ylab = "Differences", bty = "l")
abline(v = results$changePoints, lty = 2)

# Figure 5c) (run above code first)
ewma <- numeric(length(diffs))
ewma[1] <- diffs[1]^2
lambda <- 0.995
for (i in 2:length(ewma)) {
  ewma[i] <- lambda * ewma[i - 1] + (1 - lambda) * diffs[i]^2
}
plot(ewma, type = "l", xlab = "Observation", ylab = "", bty = "l")
abline(v = results$changePoints, lty = 2)


# Manipulating CPMs as S4 objects]

detectiontimes <- numeric()
changepoints <- numeric()

cpm <- makeChangePointModel(cpmType = "Mood", ARL0 = 5000, startup = 20)
i <- 0
while (i < length(diffs)) {
  i <- i + 1
  cpm <- processObservation(cpm, diffs[i])
  if (changeDetected(cpm)) {
    cat(sprintf("Change detected at observation %d\n", i))
    detectiontimes <- c(detectiontimes, i)
    Ds <- getStatistics(cpm)
    tau <- which.max(Ds)
    if (length(changepoints) > 0) 
      tau <- tau + changepoints[length(changepoints)]
    changepoints <- c(changepoints, tau)
    cpm <- cpmReset(cpm)
    i <- tau
  }
}
changepoints
