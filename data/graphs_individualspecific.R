load("S100k_round1.RData")

library(mcclust.ext)
library(fields)

#############################
## Subject-specific graphs ##
#############################

graph_probs = out_mcmc$graph_probs
dim(graph_probs)
vi    =  minVI(out_mcmc$simil_mat)$cl

i_cl1 = sample(which(vi == 1), 1)
i_cl2 = sample(which(vi == 2), 1)


prob_1 = graph_probs[,,i_cl1]
prob_2 = graph_probs[,,i_cl2]

nrow(prob_1)
setwd("/Users/laura/Desktop/2024 - clustering categorical models/real data - cardiotoxicity dataset for breast cancer patients ")

X = read.csv("data.csv")

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

# pdf("probs_two.pdf", width = 9.8, height = 4.8)

par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(1,1,0.1,1))
#set.panel(1,2)

colori = colorRampPalette(c('white','black'))
pdf(file = "estimated_graph_cl1.pdf",   # The directory you want to save the file in
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

pdf(file = "estimated_graph_cl2.pdf",   # The directory you want to save the file in
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
prob_2[13, 1, 2]

# dati i dag con la posterior probability of edge inclusion di ciasucna freccia recupera
# stima dei parametri 

#########################
#####•# plot DAGs #######
#########################
library(network)
dag_i_cl1 = (prob_1 > 0.3)*1
dag_i_cl2 = (prob_2 > 0.5)*1

D1 = network(dag_i_cl1, label = colnames(X))
D2 = network(dag_i_cl2, label = colnames(X))

vertex_col = rep("white", ncol(dag_i_cl1))
# vertex_col[1:4] = "grey60" # per i nodi di interevento 
out_net = plot.network(D1, displaylabels = TRUE, label = colnames(X),vertex.col = vertex_col,
                       mode = "circle",
                       #coord = coord,
                       label.pos = 6,
                       usecurve = T, edge.curve = 0, vertex.cex = 1, edge.col = "gray40",
                       label.cex = 0.75, edge.lwd = 0.01, arrowhead.cex = 0.8)

out_net = plot.network(D2, displaylabels = TRUE, label = colnames(X), vertex.col = vertex_col,
                       mode = "circle",
                       #coord = coord,
                       label.pos = 6,
                       usecurve = T, edge.curve = 0, vertex.cex = 1, edge.col = "gray40",
                       label.cex = 0.75, edge.lwd = 0.01, arrowhead.cex = 0.8)
