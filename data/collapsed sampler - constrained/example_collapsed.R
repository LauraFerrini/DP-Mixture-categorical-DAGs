setwd("/Users/laura/Desktop/2024 - clustering categorical models/collapsed sampler - constrained") 

source("generate_data_b.R")
source("GIBBS_collapsed_rcpp_constrained.R")

# Generate two DAGs of q nodes
q = 4

adj.matr1 = matrix( c(0,0,0,0,
                      0,0,0,0,
                      0,1,0,0,
                      0,1,0,0), nrow = q, ncol = q, byrow = T)

adj.matr2 = matrix( c(0,0,0,0,
                      0,0,0,0,
                      1,0,0,0,
                      1,0,0,0), nrow = q, ncol = q, byrow = T)

# Generate n_k observations from the two dags 

Y1_data = gen_data_from_dag(123, 200, adj.matr1, delta = TRUE, alpha = 0.2, b = 0)
Y2_data = gen_data_from_dag(1, 200, adj.matr2, delta = TRUE, alpha = 0.2, b = 0)

Y1 = Y1_data$Y
Y2 = Y2_data$Y

Y = rbind(Y1, Y2)

# Suppose to have some constraints on the adjacency matrices: 
# For example 2 -> 1 is not allowed
a.constr = matrix(0, q, q)
a.constr[2,1] = NA

########################
#### MCMC arguments ####
########################

S = 2500
burn_in = 500
a_pi = 1

# b_pi = 2*(q-2)/3
b_pi = 2*q

a_alpha = 1
b_alpha = 3

# ne = NULL per il moneto senza constraints sul numero massimo di neighboors per node
# se si vuole includere va modificata la funzione move constrained
a = 1

I.cal  = sapply(1:ncol(Y), function(j) length(unique(Y[,j]))) # gives the number of levels for each var
Y = as.matrix(Y)
t0 = proc.time()
out_mcmc = Gibbs_collapsed(Y = Y, S = S, burn_in = burn_in, a_pi = a_pi,
                                 b_pi = b_pi, a_alpha = a_alpha, a = a,
                                 b_alpha = b_alpha, A_constr = a.constr)
t1 = proc.time() - t0
out_mcmc$graph_probs[2,1, ]


n = nrow(Y)
simil_probs = out_mcmc$simil_mat
graphs_probs = out_mcmc$graph_probs

library(fields)
colori = colorRampPalette(c('white','black'))
par(mar = c(4,4,1,2), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(2,0.5,0))
image.plot(t(simil_probs), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "i'", ylab = "i", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)

axis(1, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)
axis(2, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)

