#########################
## Auxiliary functions ##
#########################

library(gRbase)

## find parents of node j in dag

pa = function(j, dag){
  ifelse(all(dag[,j] == 0),
         return(NULL),
         return(as.numeric(which(dag[,j] != 0))))
}

## find the family of node j in dag

fa = function(j, dag){
  return(as.numeric(c(j, which(dag[,j] != 0))))
}


#########################################
## Compute the DAG marginal likelihood ##
#########################################

## In our MCMC scheme we compare two DAGs (Dag and Dag_star) which differ locally (by the insertion, removal, reversal of a single edge u -> j)
## Because the marginal likelihood admits a node-by-node factorization this will simplify in the ratio between Dag and Dag_star except for component j
## We can implement a formula which computes the marginal likelihood relative to component j only

marg_j = function(j, a, N, I.cal, I.size, dag){

  ## Notice that there are configurations that are never observed in the data and therefore not contained in table N
  ## We need to include for these configurations the corresponding zero frequencies in the formulas below
  ## Therefore I distinguish between observed and unobserved configurations/frequencies
  
  # Observed frequencies for each level of Yj and its parents (i.e. variables in fa(j,dag))
  
  pa.j = pa(j, dag)
  fa.j = fa(j, dag)
  
  N.j.obs   = aggregate(N$freq, by = as.list(N[fa.j]), FUN = sum)$x  # frequncies of the combinations of leveles that we observe
  N.j.zeros = rep(0, prod(I.cal[fa.j]) - length(N.j.obs)) # replicate 0 frequencies for all the combinations of levels that are not observed
  
  N.j = c(N.j.obs, N.j.zeros)
  
  if(length(pa.j) == 0){
    
    return((lgamma(a) - lgamma(a + sum(N.j))) + 
             sum(lgamma(a/length(N.j) + N.j)) - length(N.j)*lgamma(a/length(N.j)))
    
  }else{
    
    N.paj.obs   = aggregate(N$freq, by = as.list(N[pa.j]), FUN = sum)$x
    N.paj.zeros = rep(0, prod(I.cal[pa.j]) - length(N.paj.obs))
    
    N.pa = c(N.paj.obs, N.paj.zeros)
    
    return((length(N.pa)*lgamma(a/length(N.pa)) - sum(lgamma(a/length(N.pa) + N.pa))) + 
             sum(lgamma(a/length(N.j) + N.j)) - length(N.j)*lgamma(a/length(N.j)))
    
  }
  
}
