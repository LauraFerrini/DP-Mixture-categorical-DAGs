library(parallel)
library(klaR)
library(poLCA)

# simulations over samples

nc   = detectCores()   # detect cores
clus = makeCluster(min(nc,maxcore)) # create cluster

parallel::clusterExport(clus, ls(), envir = .GlobalEnv) # export vars to nodes

# memory.size(50000)

out = parLapply(clus, 1:N, function(l){
  
  library(mvtnorm)
  library(parallel)
  library(klaR)
  library(poLCA)
  
  ## Simulate data
  
  source("simulations/generate_data.R")
  
  set.seed(l)
  dag_1 = (t(as(randomDAG(n = q, prob = w), "matrix")) != 0)*1
  Y1_data = gen_data_from_dag(l, n_k, dag_1, alpha = alpha)
  
  set.seed(l+1)
  dag_2 = (t(as(randomDAG(n = q, prob = w), "matrix")) != 0)*1
  Y2_data = gen_data_from_dag(l+1, n_k, dag_2, alpha = alpha)
  
  Y_data  = list(Y1_data = Y1_data, Y2_data = Y2_data, dag_1 = dag_1, dag_2 = dag_2)
  
  Y1 = Y1_data$Y
  Y2 = Y2_data$Y
  
  Y = rbind(Y1, Y2)
  Y = as.matrix(Y)
  
  
  ## Method 1 : Collapsed sampler: dag based
  
  a_pi = 1
  b_pi = 2*q
  
  a_alpha = 3
  b_alpha = 1
  
  #ne = q/2
  a  = 1
  
  I.cal  = sapply(1:ncol(Y), function(j) length(unique(Y[,j]))) # gives the number of levels for each var
  
  source("MCMC/GIBBS_collapsed_rcpp.R")
  
  set.seed(123)
  out_mcmc_collapsed = Gibbs_collapsed(Y = Y, S = S, burn_in = burn, a_pi = a_pi,
                                       b_pi = b_pi, a_alpha = a_alpha, a = a,
                                       b_alpha = b_alpha, A_constr = NULL)
  
  
  out_simil_dag = out_mcmc_collapsed$simil_mat
  out_probs_dag = out_mcmc_collapsed$graph_probs
  
  
  ## Method 2: Oracle -> only for structure learning
  
  xi = c(rep(1, nrow(Y1)), rep(2, nrow(Y2)))
  
  source("MCMC/Gibbs_collapsed_oracle.R")
  
  set.seed(123)
  out_oracle = Gibbs_collapsed_ORACLE(Y, S, burn_in = burn, a_pi = a_pi, b_pi = b_pi, ne = NULL, 
                                      a = a, xi = xi, A_constr = NULL)
  
  out_probs_dag_oracle = out_oracle$graph_probs
  
  
  ## Method 3: No dags -> the update of dags is replaced by proposing always an empty dag 
  
  source("MCMC/Gibbs_collapsed_nodags.R")
  
  set.seed(123)
  out_mcmc_nodags = Gibbs_collapsed_no_dags(Y = Y, S = S, burn_in = burn,
                                            a_alpha = a_alpha, a = a,
                                            b_alpha = b_alpha, ne = NULL)
  
  out_simil_nodag = out_mcmc_nodags$simil_mat
  
  
  ## Method 4: K Modes
  
  set.seed(123)
  out_kmodes = kmodes(Y, 2)
  
  
  ## Method 5: Latent class model
  
  Y_polca = data.frame(Y + 1)
  set.seed(123)
  out_polca =  poLCA(cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)~ 1, data = Y_polca, nclass = 2, nrep = 5, verbose = F)
  
  
  return(out = list(data    = Y_data,
                    xi_true = xi,
                    out_simil_dag = out_simil_dag,
                    out_probs_dag = out_probs_dag,
                    out_probs_dag_oracle = out_probs_dag_oracle,
                    out_simil_nodag      = out_simil_nodag,
                    out_kmodes = out_kmodes$cluster,
                    out_polca = out_polca$predclass
  ))
  
})
