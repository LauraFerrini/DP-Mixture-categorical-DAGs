# bnlearn version required: to be downloaded from github https://github.com/mldv/bnlearn
# the following code requires to use
# bn.fit(..., method = "bayes-sample") which is adapted from bnlearn package. 
library(devtools)

library(bnlearn)
library(gtools)

sample_theta <- function(burn_in, S, Xi_chain, DAG_chain, X, q, seed) {
  
  Theta_chain = vector(mode = "list", length = S - burn_in)
  bn.tmp = empty.graph(as.character(paste0("X",1:q)))
  
  X = data.frame(X)
  X[,1:ncol(X)] = lapply(X[,1:ncol(X)], as.factor)
  
  colnames(X) = as.character(paste0("X", 1:q))
  
  r = 1
  
  for(s in (burn_in + 1):S){
    
    Theta = list()
    
    set.seed(seed)
    
    xi = Xi_chain[,s]
    
    for(k in unique(xi)){
      
      Xk = X[xi == k, TRUE]
      Dk = DAG_chain[[s]][,,k]
      
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
  return(Theta_chain)
}