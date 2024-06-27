library(mcclust.ext)
library(fields)

load("data/out_breastcancer1.RData")
#############################
## Subject-specific graphs ##
#############################

graph_probs = out_mcmc$graph_probs

vi    =  minVI(out_mcmc$simil_mat)$cl
# draw at random 2 individuals
set.seed(1234)
i_cl1 = sample(which(vi == 1), 1)
set.seed(5678)
i_cl2 = sample(which(vi == 2), 1)

prob_1 = graph_probs[,,i_cl1]
prob_2 = graph_probs[,,i_cl2]

X = read.csv("data/data.csv")

ncol(X)
X = X[,-1]
X = X[, c(3,2,1,4:21)]
colnames(X)

new_names = unname(sapply(colnames(X), function(x) gsub("_", " ", x)))
colnames(X) = new_names
n = nrow(X)
q = ncol(X)

rownames(prob_1) = colnames(prob_1) = colnames(X)
colnames(prob_2) = rownames(prob_2) = colnames(X)


colori = colorRampPalette(c('white','black'))
pdf(file = "data/PPI_cl1.pdf",   # The directory you want to save the file in
     width = 9, # The width of the plot in inches
     height = 7.2) # The height of the plot in inches
par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(1,1,0.1,1))
image.plot(t(prob_1), col = colori(100), zlim = c(0,1), cex.sub = 1.5, xlab = "", 
           axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 20, cex = 1, cex.axis = 1)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1, cex.axis = 1)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)
dev.off()

pdf(file = "data/PPI_cl2.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 7.2) # The height of the plot in inches
par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(1,1,0.1,1))
image.plot(t(prob_2), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2, cex.axis = 1)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2, cex.axis = 1)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)
dev.off()

