## Recover posterior distribution of DAG parameters (and so causal effects)
library(bnlearn)
library(gtools)

load("data/out_breastcancer1.RData")
load("data/theta_function.R")
X = read.csv("data/breast_cancer.csv")
colnames(X)
n = nrow(X)
q = ncol(X)

X = data.frame(X)
X[,1:ncol(X)] = lapply(X[,1:ncol(X)], as.factor)
colnames(X) = as.character(paste0("X",1:q))

S = 100000
burn_in = 10000
seed = 1234



save(Theta_chain, file = "Theta_chain.RData")

