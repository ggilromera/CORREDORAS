rm(list = ls())
setwd("~/Dropbox/Projects/108_CORREDORAS/src/")
## load package
require(rEDM)
require(viridis)

## load data
d_sp = read.csv("../data/BSM_age_pollen_temp_fire/bsm_raw_pollen.csv")
d_age = d_sp[,1]
d_temp = read.csv("../data/BSM_age_pollen_temp_fire/BSM_Temp.csv")

d_sp = d_sp[,colSums(d_sp,na.rm=TRUE)>0][,-1]


# for demostration purposes - keep only the first 10 species
d_sp = d_sp[,1:10]



## interpolate datasets
# species
timesteps = seq(min(d_age), max(d_age), length = nrow(d_sp))

d_sp_int = matrix(nrow = nrow(d_sp), ncol = ncol(d_sp), data = NA)

for(i in 1:ncol(d_sp)) {
  mod = approxfun(d_age, d_sp[,i], yleft = 0, yright = 0)
  d_sp_int[,i] = pmax(mod(timesteps), 0)
}

# check species
matplot(d_age, d_sp, type = "l", col = viridis(ncol(d_sp)))
matplot(timesteps, d_sp_int, type = "l", col = viridis(ncol(d_sp)))

# temp
d_temp$WAPLS_C2 = as.numeric(d_temp$WAPLS_C2)
mod = approxfun(d_temp$age, d_temp$WAPLS_C2, yleft = d_temp$WAPLS_C2[1], yright = d_temp$WAPLS_C2[nrow(d_temp)])
d_temp_int = pmax(mod(timesteps), 0)

plot(d_temp$age, d_temp$WAPLS_C2, type = "l")
lines(timesteps, d_temp_int, col = 2)

## run CCM on all species pairs
CCM_out = matrix(nrow = ncol(d_sp), ncol = ncol(d_sp), data = NA)

print(paste("num_sp =", ncol(d_sp)))
for(i in 1:(ncol(d_sp)-1)) {
  for(j in (i+1):ncol(d_sp)) {
    
    dtmp = data.frame(A = d_sp[,i], B = d_sp[,j])
    dtmp[is.na(dtmp)] = 0
    dtmp = data.frame(time = 1:nrow(dtmp), dtmp)
    
    # get embedding dimensions
    Emax = 10
    simplex_out = data.frame(rho = rep(NA, Emax))
    for(ii in 1:Emax) {
      tmp = Simplex(dataFrame = dtmp, E = ii, columns = "A", target = "A", lib = paste("1", nrow(dtmp)), pred = paste("1", nrow(dtmp)))
      simplex_out$rho[ii] = cor(tmp$Observations, tmp$Predictions, use = "pairwise.complete.obs")
    }
    EA = (1:Emax)[which.max(simplex_out$rho)]
    
    simplex_out = data.frame(rho = rep(NA, Emax))
    for(ii in 1:Emax) {
      tmp = Simplex(dataFrame = dtmp, E = ii, columns = "B", target = "B", lib = paste("1", nrow(dtmp)), pred = paste("1", nrow(dtmp)))
      simplex_out$rho[ii] = cor(tmp$Observations, tmp$Predictions, use = "pairwise.complete.obs")
    }
    EB = (1:Emax)[which.max(simplex_out$rho)]
    
    # does B causally force A? use embedding dimension from B
    niter = 100
    ccm_out_BcauseA = numeric(niter)
    for(ii in 1:niter) {
      ccm_out_BcauseA[ii] = CCM(dataFrame = dtmp, E = EB, columns = "A", target = "B", libSizes = nrow(d_sp)-EB, sample = 1)[1,2]
    }
    
    # does A causally force B? use embedding dimension from A
    ccm_out_AcauseB = numeric(niter)
    for(ii in 1:niter) {
      ccm_out_AcauseB[ii] = CCM(dataFrame = dtmp, E = EA, columns = "B", target = "A", libSizes = nrow(d_sp)-EA, sample = 1)[1,2]
    }
    
    
    # check significance based on cross-validation
    q_BcauseA = quantile(ccm_out_BcauseA, c(0.025, 0.5, 0.975))
    q_AcauseB = quantile(ccm_out_AcauseB, c(0.025, 0.5, 0.975))
    
    if(sign(q_BcauseA[1]) == sign(q_BcauseA[3])) { # check whether 95% CI crosses zero
      # does j (B) causally force i (A)?
      CCM_out[i,j] = q_BcauseA[2]
    } else {
      CCM_out[i,j] = 0
    }
    
    if(sign(q_AcauseB[1]) == sign(q_AcauseB[3])) { # check whether 95% CI crosses zero
      # does i (A) causally force j (B)?
      CCM_out[j,i] = q_AcauseB[2]
    } else {
      CCM_out[j,i] = 0
    }
    
    print(paste("i = ", i, "; j = ", j, sep = ""))
  }
}

# Shows average information transfer from column (j) to row (i):
# roughly interpreted as correlation coefficient describing the "strength of the effect of column on row"
CCM_out



## Apply EDM to get time-varying interactions
# describes the effects of species 2:10 on species 1
focal_species = 1
CCM_out[focal_species,]

smap_data = d_sp[,c(focal_species, which(!is.na(CCM_out[1,]) & CCM_out[1,]>0))]
smap_data = data.frame(Time = 1:nrow(smap_data), smap_data)

thetalst = c(0, 0.001, 0.01, 0.1, 0.2, 0.5, 1, 2, 5) # nonlinearity parameter
smap_full = NULL
for(i in 1:length(thetalst)) {
  smap_full[[i]] = SMap(dataFrame=smap_data, theta = thetalst[i], embedded = TRUE, target = colnames(d_sp)[focal_species], columns = colnames(smap_data)[-1],
       lib = paste("1", nrow(smap_data)), pred = paste("1", nrow(smap_data)))
}
rho_out = sapply(smap_full, function(x) cor(x$predictions$Observations, x$predictions$Predictions, use = "pairwise.complete.obs")) # get fits
thetalst[which.max(rho_out)] # which theta? 0 mean linear process, >>0 means nonlinear

smap_out = smap_full[[which.max(rho_out)]]

# plot of effects of each species on the focal species (Abies) over time
matplot(timesteps, smap_out$coefficients[-1,-c(1,2)], type = "l",
        xlab = "time", ylab = "effect on focal species", col = viridis(ncol(smap_data)))
legend("topright", colnames(smap_out$coefficients[-c(1:2)]),#
       lty = 1,
       col = viridis(ncol(smap_data)))

