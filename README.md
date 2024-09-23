# Bayesian nonparametric mixtures of categorical directed graphs for personalized causal inference

This repository contains the R codes implementing DP mixture of categorical DAGs.

## MCMC
The folder MCMC contains the R codes for the implementation of the main MCMC algorithms for the DP mixture of categorical DAGs.  
In particular: 
  * GIBBS_collapsed_rcpp.R       : contains the main MCMC algorithm for posterior inference
  * move_dag.R                   : implements the proposal distribution of DAGs
  * sample_from_baseline.R       : samples from the baseline over the space of DAGs
  * marg_dag.R                   : computes the marginal likelihood
  * prob_ik_nonempty_function.R  : computes the probability of allocating an individual to a non empty cluster
  * gamma_causal.R               : computes the causal effects at subject-specific level

  * Gibbs_collapsed_nodags.R     : implements the **no DAG** version of our MCMC algoriithm  
  * Gibbs_collapsed_oracle.R     : implements the **ORACLE** version of our MCMC algorithm 
  * theta_function.R             : draws from the posterior of DAG parameters
  * GIBBS_joint_rcpp.R           : contains the collapsed MCMC and retrieves the DAG parameters 
  
## Data
The folder data contains the R codes for the the analysis of cardiac side effects induced by anticancer treatments on breast cancer patients. 

 * breast_cancer.csv    : contains the data used in the analysis. 
 * Run_MCMC.R           : implements the MCMC on breast cancer patients and it produces plots for the analysis on the results
