setwd("/home/graciela/nextcloud/saco.csic.es/1 graciela/projects/Plan nacional/CORREDORAS/data/chronology_2025/BSM24")
# setwd es la función que dirije a una 
# localización concreta, hay que escribir la ruta 
# donde tengamos guardado el archivo .csv de edades

## Installing and loading needed packages.

install.packages(c("clam","rbacon", "dplyr", "tidyverse"))
libraries <- c("clam","rbacon", "dplyr", "tidyverse")
lapply (libraries, library, character.only = TRUE)

#we define the parameters for Rbacon
thick<-20
d.by<-1
core.name="bsm24_v2"
acc.mean=10
acc.shape=2
mem.strength=0.5
mem.mean=0.5

# mem.strength * (1 - mem.mean) #should be smaller than 1

#Bacon(core.name, acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean, mem.strength=mem.strength, thick=thick, d.by=d.by, ask=FALSE, suggest=FALSE, plot.pdf=TRUE, depths.file=FALSE, normal=FALSE, rotate.axes=TRUE)

x <- Bacon(
  coredir = "/home/graciela/nextcloud/saco.csic.es/1 graciela/projects/Plan nacional/CORREDORAS/data/chronology_2025/BSM24/",
  core.name,
  acc.mean=acc.mean,
  acc.shape=acc.shape,
  mem.mean=mem.mean,
  mem.strength=mem.strength,
  thick=thick,
  d.by=d.by,
  ask=FALSE,
  suggest=FALSE,
  plot.pdf=TRUE,
  depths.file=FALSE,
  normal=FALSE,
  rotate.axes=TRUE)

##ACCUMULATION RATE PLOTS

accrate.age.ghost() #estimates yr/cm in an age scale
accrate.depth.ghost() #estimates yr/cm in a depth scale

##CALIBRATING ORIGINAL DATES

#functions to compute weighted mean and standard deviation (from SDMTools)
wt.mean <- function(x, wt) {
  s = which(is.finite(x*wt)); wt = wt[s]; x = x[s] #remove NA info
  sum(wt * x)/sum(wt)
}

wt.sd <- function(x, wt) {
  s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
  xbar = wt.mean(x,wt) #get the weighted mean
  sqrt(sum(wt *(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2))))
}

#extracting calibration object from info
calibration <- info$calib$probs

#every object in calibration is a matrix. The first column is the age, and the second column is the probability
calibration[[1]]

#getting ages and probs from the matrix
ages <- calibration[[1]][, 1]
probs <- calibration[[1]][, 2]

#calibration curve
plot(ages, probs, type = "l")

#weighted mean of first date
wt.mean(x = ages, wt = probs)

#weighted standard deviation of first date
wt.sd(x = ages, wt = probs)

#the important part

#weighted mean for every date
calibration.mean <- lapply(
  calibration,
  FUN = function(x) wt.mean(x = x[, 1], wt = x[, 2])
) %>%
  unlist()

#weighted sd for every date
calibration.sd <- lapply(
  calibration,
  FUN = function(x) wt.sd(x = x[, 1], wt = x[, 2])
) %>%
  unlist()

#to the dates data frame
df <- read.table(
  "/home/graciela/nextcloud/saco.csic.es/1 graciela/projects/Plan nacional/CORREDORAS/data/chronology_2025/BSM24/bsm24_v2/bsm24_v2.csv",
  header = TRUE,
  sep = ","
) %>%
  dplyr::mutate(
    calibration.mean = calibration.mean,
    calibration.sd = calibration.sd
  )

#Saving the dataframe to a table
setwd("/home/graciela/nextcloud/saco.csic.es/1 graciela/projects/Plan nacional/CORREDORAS/data/chronology_2025/BSM24/bsm24_v2/")
write.table(df, file = "bsm24_v2_cal.csv", row.names = FALSE, col.names = TRUE, sep=",")
