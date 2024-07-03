library(pcalg)
library(gRbase)

gen_data_from_dag = function(i, n, dag, alpha){
  
  # This function randomly generates a categorical binary dataset from a given DAg
  
  ###########
  ## INPUT ##
  ###########
  
  # i     : a numerical seed, required for reproducibility
  # n     : the number of observations
  # dag   : (q, q) adjacency matrix of a dag 
  # alpha
  
  ############
  ## OUTPUT ##
  ############
  
  # Y   : the (n,q) generated categorical (binary) dataset
  
  # If not specified otherwise, cut-offs are set equal to zero for each variable
  
  q = ncol(dag)
  
  set.seed(i)
  
  B = dag*matrix(runif(q*q, 1.1, 2), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE); diag(B) = 1
  
  Sigma_cond = diag(rep(1, q))
  
  Sigma = solve(t(B))%*%Sigma_cond%*%solve(B); mu = c(rep(0, q))
  mu    = runif(n = q, -b, b)
  
  library(mvtnorm)
  
  Z = data.frame(rmvnorm(n, mu, Sigma))
  Y = Z
  
  for(j in 1:q){
    
    gamma_j = runif(1, quantile(Z[,j], alpha), quantile(Z[,j], 1-alpha))
    Y[,j][Z[,j] >= gamma_j] = 1; Y[,j][Z[,j] < gamma_j] = 0
      
  }
  
  return(out_data = list(Y = Y, Z = Z, Sigma = Sigma))
}




