mcmc_pooled = function(Y, S, a, a_pi, b_pi, verbose = FALSE, A_constr = NULL){
  
  # This function implements our main MCMC scheme for joint posterior inference on
  # DAGs, DAG-parameters and causal effects
  
  ###########
  ## INPUT ##
  ###########
  
  # Y : an (n,q) categorical dataset
  # S : the number of MCMC iterations
  # a : common prior hyperparameter of the DAG-Dirichlet prior
  # a_pi, b_pi : hyperparameters of the prior on pi (probability of edge inclusion)
  # ne : a maximum number of node-neighbors imposing a sparsity constraint on the DAG space
  # v_set : a set of nodes for which the causal effect of do(Xv = x) on y = 1 is required; if null, all nodes are included
  
  # Required auxiliary functions:
  source("MCMC/move_dag.R")
  source("MCMC/marg_dag.R")
  
  
  n = nrow(Y)
  q = ncol(Y)
  
  if(is.null(A_constr)){
    
    A_constr = matrix(0, q, q)
    
  }
  
  bn_tmp = empty.graph(as.character(paste0("X",1:q))) # needed to convert adjacency matrices into bn objects
  
  
  Y = as.data.frame(Y)
  Y[,1:ncol(Y)] = lapply(Y[,1:ncol(Y)], as.factor)
  
  colnames(Y) = as.character(paste0("X",1:q))
  
  # Store space for parameters
  
  Dags = array(NA, dim = c(q, q, S))
  
  # matrix collecting the adjacency matrices of the DAGs
  # (each column is the by-column-vectorized adjacency matrix of the accepted DAG) ??
  
  Theta = list()
  
  # list collecting draws from the posterior of node-parameters
  
  #Out_causal = matrix(NA, S, length(vset))
  
  # matrix collecting draws from the posterior of causal effect coefficients
  
  # Inits values
  
  Dag = matrix(0, q, q); colnames(Dag) = rownames(Dag) = 1:q; Dags[,,1] = Dag
  
  N = plyr::count(Y) # gives the counts of observations having that specific combinations of levels
  
  I.cal  = sapply(1:ncol(Y), function(j) length(unique(Y[,j]))) # this should count
  # the levels assumed by each variable 
  
  I.cal.size  = sapply(1:ncol(Y), function(j) length(unique(Y[,j]))) # da rimuovere
  
  I.size = prod(I.cal) 
  
  # MCMC iterations
  
  pb = txtProgressBar(min = 2, max = S, style = 3)
  
  for(t in 1:S){
    
    
    ###################
    ## 1. Update DAG ##
    ###################
    
    ## 1.1 Propose new DAG
    
    Dag_move = move(A = Dag, q = q, A_constr = A_constr)
    
    Dag_star      = Dag_move$A_new           # adjacency matrix of the proposed DAG
    nodes_star    = Dag_move$nodes           # nodes (u,v) involved in the local move leading to Dag_star
    type.operator = Dag_move$type.operator   # the type of the applied operator in the local move (1: Insert; 2: Delete; 3: Reverse edge u -> v)
    
    ## 1.2 Compute (log)priors p(Dag), p(Dag_star) and (log)prior ration
    
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
      
      marg_star = marg_j(j = nodes_star[2], dag = Dag_star, N = N, I.cal, I.size, a = a)
      marg      = marg_j(j = nodes_star[2], dag = Dag, N = N, I.cal, I.size, a = a)
      
    }else{
      
      # (3) Reverse a directed edge
      
      marg_star = marg_j(j = nodes_star[1], dag = Dag_star, N = N, I.cal, I.size, a = a) +
        marg_j(j = nodes_star[2], dag = Dag_star, N = N, I.cal, I.size, a = a)
      
      marg = marg_j(j = nodes_star[1], dag = Dag, N = N, I.cal, I.size, a = a) +
        marg_j(j = nodes_star[2], dag = Dag, N = N, I.cal, I.size, a = a)
      
    }
    
    ## 1.4 Compute the acceptance ratio
    
    ratio = min(0, marg_star - marg + logprior.ratio)
    
    ## 1.5 Accept/Reject Dag_star and update Dag
    
    if(log(runif(1)) < ratio){
      Dag = Dag_star
    }
    
    ## 1.6 Store the accepted DAG in the chain
    
    Dags[,,t] = Dag
    
    
    ############################################
    ## 2. Sample theta given the accepted DAG ##
    ############################################
    
    ## Sample from the collection of posterior Dirichlet distributions using function pa_conf_summary
    
    ## Theta[[t]] = lapply(1:q, function(j) pa_conf_summary(I.cal = I.cal, v = j, N = N, a = a, dag = Dag, I.cal.size = I.cal.size, out = "post"))
    
    ## Sample from the collection of posterior Dirichlet distributions using function bn.fit with method = "bayes-sample"
    
    colnames(Dag) = rownames(Dag) = as.character(paste0("X",1:q))
    
    amat(bn_tmp) = Dag
    
    theta = bn.fit(x = bn_tmp, data = Y, keep.fitted = TRUE, debug = FALSE, method = "bayes-sample")
    
    Theta[[t]] = theta
    
    
    ###############################
    ## 3. Compute causal effects ##
    ###############################
    
    
    #Out_causal[t,] = c(sapply(vset, function(vv) gammav(theta, y = "X1", v = sprintf("X%s", vv))))
    
    #setTxtProgressBar(pb, t)
    #close(pb)
    
  }
  
  return(out = list(Y = Y, Dags = Dags, Theta = Theta))
  
}


