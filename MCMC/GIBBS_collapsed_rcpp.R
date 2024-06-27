source("MCMC/sample_from_baseline.R")
source("MCMC/prob_ik_nonempty_function_new_rcpp.R")
source("MCMC/marg_dag.R")
source("MCMC/normalize_weights.R")


library(abind)
library(plyr)
library(prodlim)

Gibbs_collapsed <- function(Y, S, burn_in, a_pi, b_pi, a_alpha, b_alpha, ne, a) {
  
  ###############
  #### INPUT ####
  ###############
  
  # Y           : (n, q) data matrix (nb this is a numeric matrix to use the rcpp function)
  # S           : number of MCMC iterations
  # burn_in     : the burn-in 
  # a           : the common prior hyper-parameter of Thetas in the Dirichlet
  # a_alpha     : hyper-paramerameter of the concentration parameter alpha0 of the DP 
  # b_alpha     : hyper-paramerameter of the concentration parameter alpha0 of the DP 
  # ne          : maximum number of node-neighbors in the DAG -> parameter of the move() function 
  # a_pi        : hyper-parameter of the Beta prior on probability of edge inclusion pi -> sample_baseline_dags()
  # b_pi        : hyper-parameter of the Beta prior on probability of edge inclusion pi
  
  ################
  #### OUTPUT ####
  ################
  
  # DAG         : a list of length S, where each element is a (q,q,K) array containing the D_k, forall k
  # Xi          : an (n, S) matrix containing observations' assignment to clusters at each iteration of the MCMC
  # alpha_0     : an (S, 1) vecor collecting posterior draws of the concentration parameter alpha_0 of the DP
  # simil_mat   : an (n, n) posterior similarity matrix between individuals
  # graph_probs : a (q,q,n) array with n (q,q) matrices collecting subject-specific posterior probabilities of edge inclusion
  
  n = nrow(Y)
  q = ncol(Y)
  colnames(Y) = as.character(paste0("X",1:q))
  I.cal  = sapply(1:ncol(Y), function(j) length(unique(Y[,j]))) # gives the number of levels for each var
  I.size = prod(I.cal)
  
  Xi_chain = matrix(NA, n, S)
  
  A_chain = vector(mode = "list", length = S)
  
  alpha_0_chain = rep(NA, S)
  
  n_base    = S
  burn_base = burn_in
  
  #0. Get draws from the baseline of the Dags
  
  out_baseline = sample_baseline_dags(S = n_base, burn = burn_base, 
                                      q = q, a_pi = a_pi, b_pi = b_pi)$DAG_chain
  
  #####################
  # Set initial values#
  #####################
  
  K_inits = 2    ## number of clusters
  
  alpha_0 = 0.5 ## precision parameter
  alpha_0_chain[[1]] = alpha_0
  
  ## DAGs D_1, ..., D_K (empty)
  A_0 = array(0, c(q, q, K_inits))
  A_chain[[1]] = A_0
  
  ## Update of cluster allocators
  
  set.seed(123)
  
  xi = sample(K_inits, n, replace = TRUE)
  
  while(length(table(xi)) < K_inits){
    
    xi = sample(K_inits, n, replace = TRUE)
    
  }
  
  #xi = c(rep(1, nrow(Y)/2), rep(2, nrow(Y)/2))
  
  Xi_chain[,1] = xi
  
  
  Dags  = A_0
  r = table(xi)
  
  #############################
  ###### MCMC iterations ######
  #############################
  
  pb = txtProgressBar(min = 2, max = S, style = 3)
  
  for (t in 2:S){
    
    Dag_star   = out_baseline[,,sample(n_base - burn_base, 1)]
    Dags  = abind(Dags, Dag_star)
    
    K_star = dim(Dags)[3] # there will always be one empty cluster, which is the last one
    # sampled from the baseline of the dags 
    
    logProbs = matrix(nrow = n, ncol = K_star) 
    
    # I.cal da inserire come argomento anche in funzione per la predictive
    
    for(k in 1:(K_star - 1)){
      
      
      fa.k.list = lapply(1:q, function(j) fa(j, Dags[,,k]))
      pa.k.list = lapply(1:q, function(j) pa(j, Dags[,,k]))
      
     
      logProbs[,k] = sapply(1:n, function(i) prob_ik_nonempty(Y = Y[xi == k,, drop = FALSE], pa.list = pa.k.list, fa.list = fa.k.list, 
                                                              member = (xi[i] == k),
                                                              yi = Y[i,], a = a, I.cal = I.cal)) 
      
    }
    
    logProbs[, K_star] = log(alpha_0/I.size)
    
    Probs = t(sapply(1:nrow(logProbs), function(i) normalize_weights(logProbs[i,])))
    
    xi_star = sapply(1:n, function(i) sample(1:(K_star), size = 1, prob = Probs[i,]))
    
    labs = as.integer(names(table(xi_star)))
    K_star = length(labs) # effective new number of clusters
    
    Dags  = array(Dags[,,labs], c(q, q, K_star))
    
    xi_star = as.factor(xi_star); levels(xi_star) = 1:K_star   # update labels
    
    r = table(xi_star)
    xi = c(xi_star)
    K = dim(Dags)[3] # number of non-empty clusters after having assigned each obs to a cluster 
    
    
    ###############################
    ## Update of alpha_0 given K ##
    ###############################
    
    eta = rbeta(1, alpha_0 + 1, n)
    
    alpha_0 = c(rgamma(1, shape = a_alpha + K, 
                       rate = b_alpha - log(eta)), 
                rgamma(1, shape = a_alpha + K - 1, 
                       rate = b_alpha - log(eta)))[sample(c(1,2), 1, prob = c(a_alpha + K - 1, n*(b_alpha - log(eta))))]
    
    alpha_0_chain[t] = alpha_0
    
    
    ################################
    ## Update DAGs D_1, ..., D_K  ##
    ################################
    
    set = 1:K
    
    ## Potremmo recuperare gli N_k da quelli già costruiti per aggionrare xi;
    ## là andrebbero costruiti in una lista fuori dal ciclo for su k
    ## similmente per i pa e fa set con una lista di liste
    
    ## 1.1 Propose new DAG
    for (k in set){
      
      Dag = Dags[,,k]
      Dag_move = move(A = Dag, q = q, ne = ne)
      
      Dag_star      = Dag_move$A_new           # adjacency matrix of the proposed DAG
      nodes_star    = Dag_move$nodes           # nodes (u,v) involved in the local move leading to Dag_star
      type.operator = Dag_move$type.operator   # the type of the applied operator in the local move (1: Insert; 2: Delete; 3: Reverse edge u -> v)
      
      Y_k =(Y[xi == k, , drop = FALSE])
      N_k = plyr::count(Y_k)
      
      ## 1.2 Compute (log)priors p(Dag), p(Dag_star) and (log)prior ratio
      
      logprior.star = lgamma(sum(Dag_star) + a_pi) + 
        lgamma(q*(q-1)/2 - sum(Dag_star) + b_pi - 1)
      
      logprior = lgamma(sum(Dag) + a_pi) + 
        lgamma(q*(q-1)/2 - sum(Dag) + b_pi - 1)
      
      logprior.ratio = logprior.star - logprior
      
      ## 1.3 Compute the marginal likelihood of Dag and Dag_star
      
      ## We can distinguish among the three types of operators (1,2,3) applied in the local move
      
      if(type.operator == 1 | type.operator == 2){
        
        # (1) Insert a directed edge or
        # (2) Delete a directed edge
        
        marg_star = marg_j(j = nodes_star[2], dag = Dag_star, N = N_k, I.cal, I.size, a = a)
        marg      = marg_j(j = nodes_star[2], dag = Dag, N = N_k, I.cal, I.size, a = a)
        
      }else{
        
        # (3) Reverse a directed edge
        
        marg_star = marg_j(j = nodes_star[1], dag = Dag_star, N = N_k, I.cal, I.size, a = a) +
          marg_j(j = nodes_star[2], dag = Dag_star, N = N_k, I.cal, I.size, a = a)
        
        marg = marg_j(j = nodes_star[1], dag = Dag, N = N_k, I.cal, I.size, a = a) +
          marg_j(j = nodes_star[2], dag = Dag, N = N_k, I.cal, I.size, a = a)
        
      }
      
      ## 1.4 Compute the acceptance ratio
      
      ratio = min(0, marg_star - marg + logprior.ratio)
      
      ## 1.5 Accept/Reject Dag_star and update Dag
      
      if(log(runif(1)) < ratio){
        Dag = Dag_star
      }
      
      ## 1.6 Store the accepted DAG in the chain
      
      Dags[,,k] = Dag # fin qui dovrebbe essere ok 
      
      
    }
    
    A_chain[[t]]     = Dags
    Xi_chain[,t]     = xi
    
    setTxtProgressBar(pb, t)
    close(pb)
    
  }
  
  ################################
  #### End of MCMC iterations ####
  ################################
  # similarity matrix
  # posterior prob of edge inclusion 
  
  # 1. Similarity matrix
  simil_mat = matrix(0, nrow = n, ncol = n)
  
  for(t in (burn_in + 1):S){
    
    simil_mat = simil_mat + (matrix(Xi_chain[,t], nrow = n, ncol = n) == t(matrix(Xi_chain[,t], nrow = n, ncol = n)))*1
    
  }
  
  simil_probs = simil_mat/(S-burn_in)
  
  #2. Construct subject-specific matrices with posterior probabilities of edge inclusion
  graph_probs = array(NA, c(q, q, n))
  
  for(i in 1:n){
    
    probs_i = sapply((burn_in + 1): S, function(t) c(A_chain[[t]][,,Xi_chain[i,t]]))
    graph_probs[,,i] = matrix(rowMeans(probs_i), q, q)
    
  }
  
  return(list(DAG = A_chain, Xi = Xi_chain, alpha_0 = alpha_0_chain,
              simil_mat = simil_probs, graph_probs = graph_probs))
  
  
}

