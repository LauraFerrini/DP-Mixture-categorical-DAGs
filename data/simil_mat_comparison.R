######################
# some exploratory data analysis on cluster membership 

load("S100K_round1.RData")
out1 = out_mcmc
load("S100K_round2.RData")
out2 = out_mcmc
source("clustering_functions.R")
# vi1    =  minVI(out1$simil_mat)$cl
# vi2    =  minVI(out2$simil_mat)$cl

clus_hat1 = from_simil_to_clust(out1$simil_mat)
clus_hat2 = from_simil_to_clust(out2$simil_mat)

ordered_individuals1 = unlist(sapply(1:length(table(clus_hat1)),
                                     function(k) c(which(clus_hat1 == k))))


simil_ordered1 = out1$simil_mat[ordered_individuals1,ordered_individuals1]
simil_ordered2 = out2$simil_mat[ordered_individuals1,ordered_individuals1]


plot(out1$simil_mat, out2$simil_mat)
library(fields)
#par(mfrow=c(1,2))
colori = colorRampPalette(c('white','black'))
#image(simil_ordered1)
#image(simil_ordered2)

image.plot(simil_ordered1, col = colori(100), zlim = c(0,1), cex.sub = 1.5, xlab = "", 
           axes = F, horizontal = F, legend.shrink = 1)
# With Variation of Information
library(mcclust.ext)
vi1 = minVI(out1$simil_mat)$cl
vi2 = minVI(out2$simil_mat)$cl
sum(vi1 == 1)
sum(vi1 == 2)
ordered_ind1_vi = unlist(sapply(1:length(table(vi1)),
                                     function(k) c(which(vi1 == k))))


simil_ordered1_vi = out1$simil_mat[ordered_ind1_vi,ordered_ind1_vi]
simil_ordered2_vi= out2$simil_mat[ordered_ind1_vi,ordered_ind1_vi]

dev.off()
pdf(file = "psm_round1.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7.6) # The height of the plot in inches
image.plot(simil_ordered1_vi, col = colori(100), zlim = c(0,1), cex.sub = 3, 
           axes = F, horizontal = F, legend.shrink = 1, xlab = " Individuals (i)",
           ylab = "Individuals (i')", cex.axis = 4 )




########
n = length(vi1)
grid_labs = c(1, round(n/4), n/2, round(3*n/4), n)
pos_labs  = grid_labs/n; pos_labs[1] = 0

library(fields)


dev.off()
pdf(file = "psm_round1.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7.4) # The height of the plot in inches
colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_ordered1_vi, col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "subjects (i)",
           ylab = "subjects (i')", axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)
round(simil_ordered1_vi, 2)
dev.off()
########

dev.off()
pdf(file = "psm_round2.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7.2) # The height of the plot in inches
image.plot(simil_ordered2_vi, col = colori(100), zlim = c(0,1), cex.sub = 1.5, 
           axes = F, horizontal = F, legend.shrink = 1, 
           main = "", xlab = " Individuals (i)",
           ylab = "Individuals (i')")
dev.off()
