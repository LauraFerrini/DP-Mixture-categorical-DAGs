# DP Mixture of categorical DAGs

This repository contains the R codes implementing our DP mixture of categorical D.
Specifically:

## Folder MCMC
The folder MCMC contains the R codes for the implementation of the MCMC algorithms for the DP mixture of categorical DAGs.  \\
In particular: 
\begin{itemize}
  \item GIBBS_collapsed_rcpp.R : contains the main MCMC algorithm for posterior inference
  
  \item move_dag.R                   : performs one move from a DAG to an adjacent DAG (implements the proposal distribution over the space of DAGs)
  \item sample_from_baseline.R       : samples from the baseline over the space of DAGs
  \item marg_dag.R                   : computes the marginal likelihood
  \item prob_ik_nonempty_function.R  : computes the probability of allocating an individual to a non empty cluster
  \item gamma_causal.R               : computes the causal effects at subject-specific level

  \item Gibbs_collapsed_nodags.R     : implements the \textbf{no DAG} version of our MCMC algoriithm  
  \item Gibbs_collapsed_oracle.R     : implements the \textbf{ORACLE} verion of our MCMC algorithm 
  \item theta_function.R             : draws from the posterior of DAG parameters
  \item GIBBS_joint_rcpp.R           : contains the collapsed MCMC and retrieves the DAG parameters 
  
\end{itemize}


## Folder data
The folder data contains the R codes for the quantification of cardiac side effects of anti-cancer therapies (antracyclines) on breast cancer patients. 

\begin{itemize}
  \item breast_cancer.csv  : contains the data used in the analysis. Original source: 
  \item Run_MCMC           : implements the MCMC on breast cancer patients and it produces plots for the analysis on the results
\end{itemize}


## Folder simulations 
The folder simulations contains all the R codes for the replicability of our simulation study. In particular:

\begin{itemize}
  \item generate_data.R     : generates the data
  \item parallel_sim.R      : contains the function used to parallelize our simulation study
  \item sim.R               : contains all the analysed scenarios and produces all the plots in the paper
  
\end{itemize}
