
# TITLE: TraMiner - A short example to begin with (official docume --------

# Web: http://mephisto.unige.ch/pub/TraMineR/doc/TraMineR-Users-Guide.pdf, chapter2



# Package Handling --------------------------------------------------------


library(TraMineR)
library(cluster)



# Data IO -----------------------------------------------------------------

data(mvad)
head(mvad)


# Defining a vector containing the legends for the states to appear in the graphics and creatinga sequence object which will be used as argument to the next functions 
mvad.labels <- c("employment", "further education", "higher education", "joblessness", "school", "training")
mvad.scode <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, 17:86, states = mvad.scode, labels = mvad.labels, xtstep = 6)
class(mvad.seq)


# Plots -------------------------------------------------------------------

par(mfrow = c(2,2))
# Drawing in a single figure
seqiplot(mvad.seq, withlegend = F, title = "Index plot (10 first sequences)",border = NA)

# The sequence frequency plot of the 10 most frequent sequences with bar width proportional to the frequencies
seqfplot(mvad.seq, withlegend = F, border = NA, title = "Sequence frequency plot")

# The state distirbution by time points
seqdplot(mvad.seq, withlegend = F, border = NA, title = "State distribution plot")

# The legend as a separate graphic since several plots use the same color codes for the states
seqlegend(mvad.seq, fontsize = 1.3)

# Plot the entropy of the state distribution of each time point
seqHtplot(mvad.seq, title = "Entropy index")

# Compute, summarize and plot the histogram of the sequence turbulences
Turbulence <- seqST(mvad.seq)
summary(Turbulence)
hist(Turbulence, col = "cyan", main = "Sequence turbulence")




# Typology of the trajectories --------------------------------------------

# cluster library needed 
# build Ward hierarchical clustering of sequences from the optimal matching distanced
# and retrieve for each individual sequence the cluster membership of the 4 class solution

# Compute the optimal matching distances using substitution costs based on transition rates
# observed in the data 
submat <- seqsubm(mvad.seq, method = "TRATE")
dist.om1 <- seqdist(mvad.seq, method = "OM", indel = 1, sm = submat)

clusterward1 <- agnes(dist.om1, diss = TRUE, method = "ward")
plot(clusterward1)
cl1.4 <- cutree(clusterward1, k = 4)
cl1.4fac <- factor(cl1.4, labels = paste("Type", 1:4)) 

# Plot the state distirbution at each time point within each cluster
seqdplot(mvad.seq, group = cl1.4fac, border = NA)

# Plot the sequence frequencies within each cluster
seqfplot(mvad.seq, group = cl1.4fac, border = NA)




# Event sequence analysis -------------------------------------------------

# Instead of focusing on sequences of states, we will look at
# sequences of transitions or events

# 1 - Define the sequences of transitions
mvad.seqe <- seqecreate(mvad.seq)

# 2 - Look for frequent event subsequences and plot the 15 most frequent ones
fsubseq <- seqefsub(mvad.seqe, pMinSupport = 0.05)
plot(fsubseq[1:15], col = "tomato")

# 3 - Determine the most disriminating transitions between clusters
# and plot the frequencies by cluster of the 6 first ones
discr <- seqecmpgroup(fsubseq, group = cl1.4fac)
plot(discr[1:6])








