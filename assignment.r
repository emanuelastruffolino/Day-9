###############################
######## Assignment 9 #########
###############################

rm(list = ls())
#setwd("/Users/emanuelastruffolino/Desktop/SequenceCourse/Assignments")
#getwd()
library (TraMineRextras)

#  Discrepancy between sequences and representative sequences

#1. Consider the state sequence object biofam.seq and the matrix dOM of 
  # pairwise OM dissimilarities based on the properties matrix considered 
  # in the previous assignments.

data(biofam)
mycol<-brewer.pal(,"RdBu")
biofam.lab<-c("Parent", "Left", "Married","Left+Marr","Child","Left+Child","Left+Marr+Child","Divorced")
biofam.shortlab<-c("P", "L", "M", "LM", "C","LC","LMC", "D")
biofam.seq <- seqdef(biofam[,10:25],states=biofam.shortlab,labels=biofam.lab)
summary (biofamseq)

properties <- matrix(c(# left, married, child, divorced
  0, 0, 0, 0,  # parent
  1, 0, 0, 0,  # left
  0, 1, .5, 0, # marr
  1, 1, 0, 0,  # left+marr
  0, 0, 1, 0,  # child
  1, 0, 1, 0,  # left+child
  1, 1, 1, 0,  # left+marr+child
  .5, 1, .5, 1 # divorced
), 8, 4, byrow=TRUE)
scm <- as.matrix(dist(properties))
indel <- .5*max(scm)
distOM <- seqdist(biofam.seq, method="OM", indel=indel, sm=scm, full.matrix = FALSE)

#2. Plot the single representative of the whole dataset for each of the 
  # ‘frequency’ and ‘centrality’ criteria.

seqrplot(biofam.seq, dist.matrix=distOM, criterion="dist", nrep=1)
seqrplot(biofam.seq, dist.matrix=distOM, criterion="freq", nrep=1)
seqlegend(biofam.seq, ncol=4, position="center")

#3. Plot the densest-neighborhood representative of the whole dataset by 
  # setting successively the neighborhood radius tsim as .1, .2, .3, .4.

seqrplot(biofam.seq, dist.matrix = distOM, tsim=.1,criterion = "density", nrep = 1, border = NA)
seqrplot(biofam.seq, dist.matrix = distOM, tsim=.2,criterion = "density", nrep = 1, border = NA)
seqrplot(biofam.seq, dist.matrix = distOM, tsim=.3,criterion = "density", nrep = 1, border = NA)
seqrplot(biofam.seq, dist.matrix = distOM, tsim=.4,criterion = "density", nrep = 1, border = NA)

#or better

tsimval <- c(.1,.2,.3,.4)
for (i in 1:4){
  seqrplot(biofam.seq, dist.matrix=distOM, criterion="density", tsim=tsimval[i], nrep=1)
}

#4. Plot the set of representatives using a neighborhood radius of .15 and 
  # setting successively the minimum wanted coverage (trep) as .25, .4, .5, .75, 1

trepval <-c(25, .4, .5, .75, 1)
for (i in 1:5){
seqrep(biofam.seq, dist.matrix = distOM, criterion = "density",
       trep = trepval [i], tsim = 0.15)
}

#5. Consider the 6-group solution of the PAM clustering of the set of biofam sequences 
   # and represent the clusters with representative plots using a neighborhood radius of 
   # .15 and a minimum coverage of 70%.

library(WeightedCluster)
weight<-attr(biofam.seq, "weight")
clu.pam <- wcKMedoids(distOM, k = 6, weight = weight)
clu.pam6 <- factor(clu.pam$clustering)
seqrplot(biofam.seq, dist.matrix=distOM, criterion="density", tsim=.15, trep=.7, group=clu.pam6, border=NA)

#6. Examine the quality measures of the representatives. How do you explain that the 
   # overall Q is negative when there is a single representative?

seqrep <- seqrep.grp(biofam.seq, mdist=distOM, group=clu.pam6, criterion="density", tsim=.15, trep=.7, ret = "rep")

# here I don't know how to compute some kind of loop: I've looked at the solution but I didn't get it... 

#7. Compute the discrepancy (pseudo-variance) of the whole set of data and the within 
   # dicrepancy in each cluster. (Tip: Currently, dissvar.grp accepts only a matrix as 
   # first argument. If you have computed the dissimilarity matrix with full.matrix=FALSE, 
   # pass the dissimilarity matrix as as.matrix(dOM))

dissvar(distOM)
dissvar.grp(as.matrix(distOM), group=clu.pam6)
