## Recover posterior distribution of DAG parameters (and so causal effects)
library(bnlearn)
library(gtools)
load("data/out_breastcancer1.RData")
X = read.csv("data/data.csv")
X = X[,-1]
X = X[, c(3,2,1,4:21)]
colnames(X)
n = nrow(X)
q = ncol(X)

X = data.frame(X)
X[,1:ncol(X)] = lapply(X[,1:ncol(X)], as.factor)

colnames(X) = as.character(paste0("X",1:q))

S = 100000
burn_in = 10000
Theta_chain = vector(mode = "list", length = S - burn_in)

bn.tmp = empty.graph(as.character(paste0("X",1:q)))
r = 1
set.seed(1234)

for(s in (burn_in + 1):S){
  
  Theta = list()
  xi    = out_mcmc$Xi[,s]
  
  for(k in unique(xi)){
    
    Xk = X[xi == k, TRUE]
    Dk = out_mcmc$DAG[[s]][,,k]
    
    colnames(Dk) = rownames(Dk) = as.character(paste0("X",1:q))
    amat(bn.tmp) = Dk
    
    thetak = suppressWarnings(bn.fit(x = bn.tmp, data = Xk, keep.fitted = TRUE, 
                                     debug = FALSE, 
                                     method = "bayes-sample"))
    
    Theta[[k]] = thetak
    
  }
  
  Theta_chain[[r]] = Theta
  r = r + 1
  if(s%%100 == 0){print(s)}
  
}


save(Theta_chain, file = "Theta_chain.RData")

