######################
# Posterior similarity Matrix with individuals ordered accordingly to estimated point 
# clustering structure, obtained via Variation of Information criterion
library(mcclust.ext)
library(fields)

load("data/out_breastcancer1.RData")
str(out_mcmc)

vi = minVI(out_mcmc$simil_mat)$cl

sum(vi == 1) # n_1
sum(vi == 2) # n_2
ordered_ind1_vi = unlist(sapply(1:length(table(vi)),
                                     function(k) c(which(vi == k))))

simil_ordered1_vi = out_mcmc$simil_mat[ordered_ind1_vi,ordered_ind1_vi]

n = length(vi)
grid_labs = c(1, round(n/4), n/2, round(3*n/4), n)
pos_labs  = grid_labs/n; pos_labs[1] = 0

pdf(file = "psm_breastcancer_def.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7.4) # The height of the plot in inches
colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_ordered1_vi, col = colori(100), zlim = c(0,1), cex.sub = 1,
           xlab = "subjects (i)", ylab = "subjects (i')", 
           axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)
dev.off()


