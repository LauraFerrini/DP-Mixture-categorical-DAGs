#setwd("/Users/laura/Desktop/2024 - clustering categorical models")
setwd("./real data - cardiotoxicity dataset for breast cancer patients ")

X = read.csv("data.csv")

setwd("./collapsed sampler - constrained") 
source("GIBBS_collapsed_rcpp_constrained.R")

library(ggplot2)
library(ggmosaic)
library(dplyr)
X = X[,-1]
X = X[, c(3,2,1,4:21)]
colnames(X)
n = nrow(X)
q = ncol(X)
A_constr = matrix(0,q,q)

A_constr[1,] = NA
A_constr[,2] = NA
S = 100000
burn_in = 10000
a_pi = 1

# b_pi = 2*(q-2)/3
b_pi = 2*q

# a_alpha = 1
# b_alpha = 3

a_alpha = 3; b_alpha = 1 

# ne = NULL per il moneto senza constraints sul numero massimo di neighboors per node
# se si vuole includere va modificata la funzione move constrained
a = 1

I.cal  = sapply(1:ncol(X), function(j) length(unique(X[,j]))) # gives the number of levels for each var
X = as.matrix(X)
t0 = proc.time()
out_mcmc = Gibbs_collapsed(Y = X, S = S, burn_in = burn_in, a_pi = a_pi,
                           b_pi = b_pi, a_alpha = a_alpha, a = a,
                           b_alpha = b_alpha, A_constr = A_constr)
t1 = proc.time() - t0

# the output of the mcmc prodìcedure has been stored in the file ""

# save(out_mcmc, file = "S100K_round2.RData")
## Recover posterior distribution of DAG parameters (and so causal effects)






