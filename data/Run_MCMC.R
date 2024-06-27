library(ggplot2)
library(ggmosaic)
library(dplyr)

source("data/collapsed sampler - constrained/GIBBS_collapsed_rcpp_constrained.R")

X = read.csv("data/data.csv")
head(X)

X = X[,-1]
X = X[, c(3,2,1,4:21)] # just some reversal of the columns
colnames(X)
n = nrow(X)
q = ncol(X)
# set the constraints
A_constr = matrix(0,q,q)

A_constr[1,] = NA
A_constr[,2] = NA

S = 100000
burn_in = 10000
# pi ~ Beta(a_pi, b_pi)
a_pi = 1
b_pi = 2*q
# Prior hyper-prameters on the concentration param of DP prior 
a_alpha = 3; b_alpha = 1 

# ne = NULL -> no constraints on the maximum number of neighborhoods per node
a = 1
I.cal  = sapply(1:ncol(X), function(j) length(unique(X[,j]))) # gives the number of levels for each var
X = as.matrix(X)
t0 = proc.time()

set.seed(1234)
out_mcmc = Gibbs_collapsed(Y = X, S = S, burn_in = burn_in, a_pi = a_pi,
                           b_pi = b_pi, a_alpha = a_alpha, a = a,
                           b_alpha = b_alpha, A_constr = A_constr)
t1 = proc.time() - t0

# the output of the mcmc algorithm has been stored in the file "out_breastcancer1.RData"

save(out_mcmc, file = "out_breastcancer1.RData")
## Recover posterior distribution of DAG parameters (and so causal effects) with the file ""






