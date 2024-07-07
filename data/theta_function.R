# insert lines for downloading new version of bnlearn 
library(devtools)
install_github("mldv/bnlearn")
theta_chain <- function(burn_in, S, out_mcmc, X, q, seed) {
  
  Theta_chain = vector(mode = "list", length = S - burn_in)
  bn.tmp = empty.graph(as.character(paste0("X",1:q)))
  r = 1
  
  for(s in (burn_in + 1):S){
    
    Theta = list()
    xi    = out_mcmc$Xi[,s]
    
    set.seed(seed)
    
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
  return(Theta_chain)
}