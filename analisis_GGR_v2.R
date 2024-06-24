
library(Hmisc)
library(igraph) 
library(ggplot2)
library(dplyr)
library(gridExtra)

dev.off()

### Network creation ###
# Data preparation


edad<-read.table('edad.txt',header=F)
datos<-read.table('data.txt',header=F) 
bsmchar<-read.table("BSM_microchar.csv", header=TRUE, sep=",")
agebsm <- as.numeric(bsmchar$Age.BP)
scoresbsm <- as.numeric(bsmchar$Z.SCORES)
datos=t(datos)# Transpose the 'datos' matrix.

# Remove columns from 'datos' where the sum of values is zero.
for (i in dim(datos)[2]:1) {
  if (sum(datos[, i]) == 0) {
    datos = datos[, -i]
    print(i)
  }
}

S = dim(datos)[2]

# Calculate Spearman correlation coefficients and p-values for the data in 'datos'.
R <- rcorr(as.matrix(datos), type = "spearman")$r
P <- rcorr(as.matrix(datos), type = "spearman")$P

# Set correlation values to zero if the p-value is greater than or equal to 0.05.
for (i in 1:S) {
  for (j in 1:S) {
    if (i != j) {
      if (P[i, j] >= 0.05) {
        R[i, j] = 0
      }
    }
    if (i == j) {
      R[i, j] = 0
    }
  }
}

# Create an undirected graph from the absolute values of the correlation matrix 'R'.
g <- graph.adjacency(abs(R), mode = "undirected", weighted = TRUE)

# Plot the graph with a specific layout and vertex settings.
plot(g, layout = layout_with_lgl, vertex.size = 10, vertex.label = NA)

##CENTRALITY
# Calculate degree centrality for each node in the graph
degree_centrality <- degree(g)

# Print degree centrality values for all nodes
print(degree_centrality)

# Visualize the degree centrality distribution
hist(degree_centrality, main = "Degree Centrality Distribution", xlab = "Degree Centrality", col = "lightblue", border = "black")

# Find nodes with the highest degree centrality
max_degree <- max(degree_centrality)
most_central_nodes <- which(degree_centrality == max_degree)
cat("Nodes with the highest degree centrality:", most_central_nodes, "\n")


# MOVING WINDOW

w = 20  # Define the window size.

# Create an empty array 'res' to store results.
res <- array(0, c(3, dim(datos)[1] - (w - 1)))

# Iterate through moving windows.
for (mw in 1:(dim(datos)[1] - (w - 1))) {
  # Calculate the mean of values in the current window.
  res[1, mw] <- rowMeans(edad[mw:(mw + w - 1)])
  
  # Extract the data window.
  win <- datos[mw:(mw + w - 1),]
  
  # Remove columns where the sum of values is zero.
  for (i in dim(win)[2]:1) {
    if (sum(win[, i]) == 0) {
      win = win[, -i]
    }
  }
  res[2,mw]<-dim(win)[2]
  
  # Calculate Spearman correlation coefficients and p-values for the window.
  Rw <- rcorr(as.matrix(win), type = "spearman")$r
  Pw <- rcorr(as.matrix(win), type = "spearman")$P
  
  # Set correlation values to zero if the p-value is greater than or equal to 0.05 or if the correlation is negative.
  for (i in 1:dim(Rw)[1]) {
    for (j in 1:dim(Rw)[2]) {
      if (i != j) {
        if (Pw[i, j] >= 0.05 || Rw[i, j] < 0) {
          Rw[i, j] = 0
        }
      }
      if (i == j) {
        Rw[i, j] = 0
      }
    }
  }
  
  # Calculate the density of the filtered correlation matrix.
  res[3, mw] <- sum(Rw > 0) / (res[2, mw] * res[2, mw] - res[2, mw])
}

#Basic plots for:
#number of species holding at least significant correlation in the network
plot(res[2,], type="l") 

# connectance measured as the potential number of correlations between taxa 
#in the data window actually present and positive after filtering. 
plot(res[3,], type="l")

#xlim <- range(0,12000)
#xmin <-0
#xmax <-12000
#x.tks <-seq(xmin, 12000,500)
#altitude.anomaly.redon <- as.numeric(Redon$Detrended.Altitude.Anomaly..m.)
#Redon <- cbind(age.redon, altitude.anomaly.redon)
#Redon <- as.data.frame(Redon)


# CI - Null Model Analysis

# Number of simulations to perform (a parameter you've set)
sim = 154

# Create an array to store the results of the null model analysis
# The array has two rows and a number of columns equal to the number of moving windows analyzed
CI <- array(0, c(2, sim))

# Perform 'sim' simulations to create null models
for (k in 1:sim) {
  # Randomly sample 'w' rows (data points) from your original data
  sel <- sample(1:dim(datos)[1], size = w)
  net_sel <- datos[sel,]
  
  # Remove columns in 'net_sel' where the sum of values is zero
  for (j in dim(net_sel)[2]:1) {
    if (sum(net_sel[, j]) == 0) {
      net_sel = net_sel[,-j]
    }
  }
  
  # Calculate and store the number of variables (columns) in the null model
  CI[1, k] <- dim(net_sel)[2]
  
  # Calculate Spearman correlation coefficients and p-values for the null model
  Rw <- rcorr(as.matrix(net_sel), type = "spearman")$r
  Pw <- rcorr(as.matrix(net_sel), type = "spearman")$P
  
  # Filter the null model's correlations in the same way as before
  for (i in 1:dim(Rw)[1]) {
    for (j in 1:dim(Rw)[2]) {
      if (i != j) {
        if (Pw[i, j] >= 0.05 || Rw[i, j] < 0) {
          Rw[i, j] = 0
        }
      }
      if (i == j) {
        Rw[i, j] = 0
      }
    }
  }
  
  # Calculate and store the complexity (density) of the null model
  CI[2, k] <- sum(Rw > 0) / (CI[1, k] * (CI[1, k] - 1))
}
#Estimate the boundaries of the quantiles

q2.25<-quantile(CI[1,],0.025)
q2.75<-quantile(CI[1,],0.975)
CI_interval_r2<-cbind(q2.25,q2.75)

q3.25<-quantile(CI[2,],0.025)
q3.75<-quantile(CI[2,],0.975)
CI_interval_r2<-cbind(q2.25,q2.75)

#Creating the confidence plot

#combined_plot2 <- grid.arrange(pres2, confidence_plot2, ncol = 1)
#print(combined_plot2)

##These data are not available yet##
#Basa de la Mora: charcoal and chironomid derived T
#charcoal
plotcharbsm <- bsmchar %>%
  ggplot(aes(x=age.BSM_microchar, y=z.scores.BSM_microchar)) +
  xlab ("age") +
  ylab ("Microcharcoal anomalies in sed.rate") +
  theme_classic() +
  geom_line() +
  geom_area(aes(), fill='#FF8C00')+
  scale_x_reverse(limits = c(12.238, -0.1), breaks=c(
    1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5,9, 9.5, 10, 10.5, 11))+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank())

#Temperature plot for BSM
BSM_Temp <- read.csv2("BSM_Temp.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
age.BSM_Temp <- as.numeric(BSM_Temp$age)
age.BSM_Temp <- age.BSM_Temp/1000
WAPLS.BSM_Temp <- as.numeric(BSM_Temp$WAPLS_C2)
BSM_Temp <- cbind(age.BSM_Temp, WAPLS.BSM_Temp)
BSM_Temp <- as.data.frame(BSM_Temp)

chiroT <- BSM_Temp %>%
  ggplot(aes(x=age.BSM_Temp, y=WAPLS.BSM_Temp))+
  xlab("") +
  ylab("Temperature (ºC)")+
  theme_classic() +
  geom_line() + scale_x_reverse(limits = c(12.238, 1.307), breaks=c(
    1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5,9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5))+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank())
