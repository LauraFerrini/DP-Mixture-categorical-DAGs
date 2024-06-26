from_simil_to_clust = function(simil_probs){
  
  simil_mat = 1*(simil_probs > 0.5)
  
  clust_ind = c()
  
  for(i in ncol(simil_mat):1){
    
    clust_ind[simil_mat[i,] == 1] = i
    
  }
  
  clust_ind = as.factor(clust_ind)
  
  levels(clust_ind) = 1:(length(levels(clust_ind)))
  
  return(clust_ind)
  
}


#######################
## Results for paper ##
#######################

load("out_leukemia_inits_K2_alpha1_n250.RData")

simil_probs = out_mcmc$simil_probs
graph_probs = out_mcmc$graph_probs

est_causal_1 = out_mcmc$est_causal$`response y = 1`
est_causal_2 = out_mcmc$est_causal$`response y = 2`
est_causal_3 = out_mcmc$est_causal$`response y = 3`


n = nrow(X)

clus_hat = from_simil_to_clust(simil_probs)

table(clus_hat)

ordered_individuals = unlist(sapply(1:length(table(clus_hat)), function(k) c(which(clus_hat == k))))

simil_ordered = simil_probs[ordered_individuals,ordered_individuals]


#######################
## Similarity matrix ##
#######################

grid_labs = c(1, round(n/4), n/2, round(3*n/4), n)
pos_labs  = grid_labs/n; pos_labs[1] = 0

library(fields)

pdf("simil_map_ordered.pdf", width = 7.1, height = 6.3)

colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_ordered, col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "subjects (i)", ylab = "subjects (i')", axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)

dev.off()


######################
## Graph estimation ##
######################

## Subject-specific graphs

prob_1 = graph_probs[,,1]
prob_6 = graph_probs[,,6]

rownames(prob_1) = rownames(prob_6) = colnames(X)
colnames(prob_1) = colnames(prob_6) = colnames(X)

pdf("probs_two.pdf", width = 9.8, height = 4.8)

par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(0.5,0.5,0.1,0.5))
set.panel(1,2)

colori = colorRampPalette(c('white','black'))

image.plot(t(prob_1), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)

image.plot(t(prob_6), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)

dev.off()



####################
## Causal effects ##
####################

## Sort individuals according to the "estimated" cluster allocations

est_causal_1_sort = est_causal_1[ordered_individuals,]
est_causal_2_sort = est_causal_2[ordered_individuals,]
est_causal_3_sort = est_causal_3[ordered_individuals,]

colnames(X)[1:3]

# to set a range (equal for all heat maps)

min(c(est_causal_1, est_causal_2, est_causal_3))
max(c(est_causal_1, est_causal_2, est_causal_3))

pdf("causal_1.pdf", width = 10, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.4))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_1_sort, col = colori(100), zlim = c(-0.5,0.5), cex.sub = 1, xlab = "subjects (i)", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, expression(AKT), outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -0.002)

dev.off()


pdf("causal_2.pdf", width = 10, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.4))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_2_sort, col = colori(100), zlim = c(-0.5,0.5), cex.sub = 1, xlab = "subjects (i)", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, expression(AKT.p308), outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -0.002)

dev.off()


pdf("causal_3.pdf", width = 10, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.4))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_3_sort, col = colori(100), zlim = c(-0.5,0.5), cex.sub = 1, xlab = "subjects (i)", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, expression(AKT.p473), outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -0.002)

dev.off()



###################################
## Comparison with no clustering ##
###################################

## For comparison without groups (no DP)

load("out_leukemia_no_mixture.RData")

probs_no_cluster = out_mcmc$graph_probs

rownames(probs_no_cluster) = colnames(X)
colnames(probs_no_cluster) = colnames(X)


## Graph estimate

library(fields)

pdf("probs_no_cluster.pdf", width = 5.6, height = 4.8)

par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(0.5,0.5,0.1,0.4))
set.panel(1,1)

colori = colorRampPalette(c('white','black'))

image.plot(t(probs_no_cluster), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)

dev.off()


####################
## Causal effects ##
####################

## Sort individuals according to the "estimated" cluster allocations

est_causal_1_no = out_mcmc$est_causal$`response y = 1`
est_causal_2_no = out_mcmc$est_causal$`response y = 2`
est_causal_3_no = out_mcmc$est_causal$`response y = 3`

est_causal_1_no_sort = est_causal_1_no
est_causal_2_no_sort = est_causal_2_no
est_causal_3_no_sort = est_causal_3_no


pdf("causal_1_no_cluster.pdf", width = 2.5, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,0.5,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_1_no_sort, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -1)

dev.off()


pdf("causal_2_no_cluster.pdf", width = 2.5, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,0.5,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_2_no_sort, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -1)

dev.off()


pdf("causal_3_no_cluster.pdf", width = 2.5, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,0.5,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_3_no_sort, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -1)

dev.off()







################################
## Diagnsotics of convergence ##
################################

load("out_leukemia_5_inits_K2_alpha1_n250.RData")

simil_probs_1 = out_mcmc$simil_probs
graph_probs_1 = out_mcmc$graph_probs

est_causal_1_1 = out_mcmc$est_causal$`response y = 1`
est_causal_2_1 = out_mcmc$est_causal$`response y = 2`
est_causal_3_1 = out_mcmc$est_causal$`response y = 3`

load("out_leukemia_6_inits_K2_alpha1_n250.RData")

simil_probs_2 = out_mcmc$simil_probs
graph_probs_2 = out_mcmc$graph_probs

est_causal_1_2 = out_mcmc$est_causal$`response y = 1`
est_causal_2_2 = out_mcmc$est_causal$`response y = 2`
est_causal_3_2 = out_mcmc$est_causal$`response y = 3`



par(mfrow = c(1,2))

image(simil_probs_1)
image(simil_probs_2)

image(graph_probs_1[,,1])
image(graph_probs_2[,,1])

sub = which(clus_hat_1 == clus_hat_2)


table(clus_hat_1)
table(clus_hat_2)


probs_1_1 = graph_probs_1[,,sub]
probs_1_2 = graph_probs_2[,,sub]



df_cfr_1 = data.frame(causal_1 = c(est_causal_1_1), causal_2 = c(est_causal_1_2))
df_cfr_2 = data.frame(causal_1 = c(est_causal_2_1), causal_2 = c(est_causal_2_2))
df_cfr_3 = data.frame(causal_1 = c(est_causal_3_1), causal_2 = c(est_causal_3_2))

library(ggplot2)

pdf(file = "causal_conv_1.pdf", width = 4, height = 3.6)

ggplot(df_cfr_1, aes(causal_1, causal_2)) +
  geom_abline(col = "grey") +
  geom_point(aes(causal_1, causal_2, col = "grey70"), size = 1.3, shape = 16) +
  labs(x = "causal effect (chain I)", y = "causal effect (chain II)") +
  scale_color_manual(values = c("gray30")) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.4,0.4), "cm")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position = "none")

dev.off()


pdf(file = "causal_conv_2.pdf", width = 4, height = 3.6)

ggplot(df_cfr_2, aes(causal_1, causal_2)) +
  geom_abline(col = "grey") +
  geom_point(aes(causal_1, causal_2, col = "grey70"), size = 1.3, shape = 16) +
  labs(x = "causal effect (chain I)", y = "causal effect (chain II)") +
  scale_color_manual(values = c("gray30")) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.4,0.4), "cm")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position = "none")

dev.off()


pdf(file = "causal_conv_3.pdf", width = 4, height = 3.6)

ggplot(df_cfr_3, aes(causal_1, causal_2)) +
  geom_abline(col = "grey") +
  geom_point(aes(causal_1, causal_2, col = "grey70"), size = 1.3, shape = 16) +
  labs(x = "causal effect (chain I)", y = "causal effect (chain II)") +
  scale_color_manual(values = c("gray30")) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.4,0.4), "cm")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position = "none")

dev.off()




simil_ordered_1 = simil_probs_1[ordered_individuals,ordered_individuals]
simil_ordered_2 = simil_probs_2[ordered_individuals,ordered_individuals]







grid_labs = c(1, round(n/4), n/2, round(3*n/4), n)
pos_labs  = grid_labs/n; pos_labs[1] = 0


library(fields)

pdf("simil_map_1.pdf", width = 7, height = 6.3)

colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_probs_1, col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "subjects (i)", ylab = "subjects (i')", axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)

dev.off()

pdf("simil_map_2.pdf", width = 7, height = 6.3)

colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_probs_2, col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "subjects (i)", ylab = "subjects (i')", axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)

dev.off()

library(fields)

pdf("simil_map_ordered_1.pdf", width = 7, height = 6.3)

colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_ordered_1, col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "subjects (i)", ylab = "subjects (i')", axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)

dev.off()


pdf("simil_map_ordered_2.pdf", width = 7, height = 6.3)

colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_ordered_2, col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "subjects (i)", ylab = "subjects (i')", axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)

dev.off()





###############################
## Comparison between groups ##
###############################

prop = c(6.5,13.8,26.8,23.4,3.4,2.3,5.4,3.8,2.7,1.9,1.9,8)/100

sum(prop)
par(mfrow = c(3,6))

labs_subtypes = data$FAB

labs_subtypes[labs_subtypes == ""] = "Unknown"

table(labs_subtypes)

tab = table(clus_hat,labs_subtypes)

PP = t(tab)/colSums(tab)

library(xtable)

xtable(round(t(PP), 2))


cfr_clus = data.frame(clus_hat = clus_hat, subtype = labs_subtypes)

library(ggplot2)
library(ggmosaic)
library(dplyr)

pdf("cfr_clusters_subtype.pdf", width = 9, height = 4)

ggplot(data = cfr_clus) +
  geom_mosaic(aes(x = product(clus_hat, subtype), fill = clus_hat)) + 
  labs(y = "Estimated cluster", x = "Subtype", title = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.4,0.4,0.4), "cm")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) +
  scale_fill_manual(values = c("gray30", "gray70")) +
  guides(x =  guide_axis(angle = 90)) +
  theme(legend.position = "none")

dev.off()

boxplot(X[,1]~clus_hat)
boxplot(X[,2]~clus_hat)
boxplot(X[,3]~clus_hat)
boxplot(X[,4]~clus_hat)
boxplot(X[,5]~clus_hat)
boxplot(X[,6]~clus_hat)
boxplot(X[,7]~clus_hat)
boxplot(X[,8]~clus_hat)
boxplot(X[,9]~clus_hat)
boxplot(X[,10]~clus_hat)
boxplot(X[,11]~clus_hat)
boxplot(X[,12]~clus_hat)



##############################
## Comparison with group M4 ##
##############################

load("out_leukemia_subtype_M4.RData")

probs_no_cluster = out_mcmc_M4$graph_probs

rownames(probs_no_cluster) = colnames(X_M4)
colnames(probs_no_cluster) = colnames(X_M4)


## Graph estimate

library(fields)

pdf("probs_subtype_M4.pdf", width = 5.6, height = 4.8)

par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(0.5,0.5,0.1,0.4))
set.panel(1,1)

colori = colorRampPalette(c('white','black'))

image.plot(t(probs_no_cluster), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)

dev.off()


## Causal effects

est_causal_1_M4 = out_mcmc_M4$est_causal$`response y = 1`
est_causal_2_M4 = out_mcmc_M4$est_causal$`response y = 2`
est_causal_3_M4 = out_mcmc_M4$est_causal$`response y = 3`

est_causal_1_M4_sort = est_causal_1_M4
est_causal_2_M4_sort = est_causal_2_M4
est_causal_3_M4_sort = est_causal_3_M4


pdf("causal_1_M4_cluster.pdf", width = 2.8, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,0.5,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_1_M4_sort, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, text = "AKT (Subtype M4)", outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -1)

dev.off()

pdf("causal_2_M4_cluster.pdf", width = 2.8, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,0.5,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_2_M4_sort, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, text = "AKT.p308 (Subtype M4)", outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -1)

dev.off()

pdf("causal_3_M4_cluster.pdf", width = 2.8, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,0.5,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(est_causal_3_M4_sort, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = NA, las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, text = "AKT.p473 (Subtype M4)", outer = F)
abline(h = 1.030)
abline(v = 1.002)
abline(h = -0.03)
abline(v = -1)

dev.off()


setwd("C:/Users/Utente/Dropbox/2020 BNP mixture DAG models/Castelletti_Consonni_codes/Application Leukemia")

load("out_leukemia_inits_K2_alpha1_n250.RData")

simil_probs = out_mcmc$simil_probs
graph_probs = out_mcmc$graph_probs

est_causal_1 = out_mcmc$est_causal$`response y = 1`
est_causal_2 = out_mcmc$est_causal$`response y = 2`
est_causal_3 = out_mcmc$est_causal$`response y = 3`

# Patients with subtype M4

causal_patient_DP_M4_1 = est_causal_1[data$FAB == "M4",]
causal_patient_DP_M4_2 = est_causal_2[data$FAB == "M4",]
causal_patient_DP_M4_3 = est_causal_3[data$FAB == "M4",]

pdf("causal_1_DP_M4.pdf", width = 6, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(causal_patient_DP_M4_1, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "subjects", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, text = "AKT (DP mixture)", outer = F)
abline(h = 1.030)
abline(v = 1.008)
abline(h = -0.03)
abline(v = -0.008)

dev.off()

pdf("causal_2_DP_M4.pdf", width = 6, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(causal_patient_DP_M4_2, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "subjects", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, text = "AKT.p308 (DP mixture)", outer = F)
abline(h = 1.030)
abline(v = 1.008)
abline(h = -0.03)
abline(v = -0.008)

dev.off()

pdf("causal_3_DP_M4.pdf", width = 6, height = 5.5)

par(mar = c(4.5,5.5,3.5,2), oma = c(0,0,0,0), mai = c(1.02,1,0.82,0.5))
set.panel(1,1)

colori = colorRampPalette(c('cyan4','white','brown'))

image.plot(causal_patient_DP_M4_3, col = colori(100), zlim = c(-0.5,0.5), legend.only = FALSE, cex.sub = 1, xlab = "subjects", axes = F, horizontal = F, legend.shrink = 1)
abline(h = -0.03)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2)
mtext(line = 1, side = 3, text = "AKT.p473 (DP mixture)", outer = F)
abline(h = 1.030)
abline(v = 1.008)
abline(h = -0.03)
abline(v = -0.008)

dev.off()


# Or comparison through scatterplots

plot(causal_patient_6_1, est_causal_1_M4, xlim = c(-0.35,0.35), ylim = c(-0.35,0.35))
plot(causal_patient_6_2, est_causal_2_M4, xlim = c(-0.35,0.35), ylim = c(-0.35,0.35))
plot(causal_patient_6_3, est_causal_3_M4, xlim = c(-0.35,0.35), ylim = c(-0.35,0.35))


##############################################
## Compare level of variables across groups ##
##############################################

ls()

summary(clus_hat)

dim(X)


out_test = matrix(NA, nrow = 3, ncol = q)

rownames(out_test) = c("Cluster 1", "Cluster 2", "p value")

colnames(out_test) = colnames(X)

for(j in 1:q){
  
  t.tmp = t.test(x = X[clus_hat == 1,j], y = X[clus_hat == 2,j])
  
  out_test[3,j] = t.tmp$p.value
  
  out_test[1:2,j] = t.tmp$estimate
  
}

round(out_test, 3)

library(xtable)

xtable(out_test)


set = c(2,4,10,11)

colnames(X)[set]

X_sorted = X[c(which(clus_hat == 1), which(clus_hat == 2)),]

X_sub_df = data.frame(AKT.p308 = X_sorted[,2],
                      BAD = X_sorted[,4],
                      BCLXL = X_sorted[,10],
                      CCND1 = X_sorted[,11],
                      Cluster = sort(clus_hat),
                      subject = 1:n)


library(ggplot2)

plot1 = ggplot(data = X_sub_df) +
  geom_point(aes(y = AKT.p308, x = subject, col = cluster), size = 1) +
  scale_x_continuous(breaks = c(1,50,100,150,200,250)) +
  ylim(-2, 4) +
  theme_bw() +
  theme(legend.position = "none")
  
plot2 = ggplot(data = X_sub_df) +
  geom_point(aes(y = BAD, x = subject, col = cluster), size = 1) +
  scale_x_continuous(breaks = c(1,50,100,150,200,250)) +
  ylim(-3, 3) +
  theme_bw() +
  theme(legend.position = "none")

plot3 = ggplot(data = X_sub_df) +
  geom_point(aes(y = BCLXL, x = subject, col = cluster), size = 1) +
  scale_x_continuous(breaks = c(1,50,100,150,200,250)) +
  ylim(-2, 6) +
  theme_bw() +
  theme(legend.position = "none")

plot4 = ggplot(data = X_sub_df) +
  geom_point(aes(y = CCND1, x = subject, col = cluster), size = 1) +
  scale_x_continuous(breaks = c(1,50,100,150,200,250)) +
  ylim(-3, 5) +
  theme_bw() +
  theme(legend.position = "none")

library(gridExtra)

pdf("cfr_groups_variables.pdf", height = 5, width = 6.5)

grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

dev.off()


## Alternative ways to estimate the clustering

library(mcclust)

devtools::install_github("sarawade/mcclust.ext")

library(mcclust.ext)

?minbinder
?minVI

binder.clust = minbinder(simil_probs)

binder.clust.ext = minbinder.ext(simil_probs)

vi.clust.ext = minVI(simil_probs)

table(binder.clust$cl, clus_hat)
table(binder.clust.ext$cl, clus_hat)
table(vi.clust.ext$cl, clus_hat)
