## Function to compute posterior probability of allocation in cluster k for subject i
library(Rcpp)

cppFunction(
  "LogicalVector compareToRow(NumericMatrix x, NumericVector y) {
  const int nr = x.nrow();
  const int nc = x.ncol();
  LogicalVector ret(nr, true);
  for (int j=0; j < nr; ++j) {
    for (int k=0; k < nc; ++k) {
      if (x(j, k) != y[k]) {
        ret[j] = false;
        break;
      }
    }
  }
  return ret;
}")




prob_ik_nonempty = function(Y, pa.list, fa.list, member, yi, a, I.cal){
  
  # Y : matrix collecting individuals assigned to cluster under evaluation
  # pa.list : list of node-by-node parents of DAG
  # fa.list : list of node-by-node families of DAG
  # member : logical (TRUE if i belongs to the cluster, FALSE otherwise)
  # yi : observation to evaluate
  
  p_ik = log(dim(Y)[1] - member) + sum(sapply(1:q, function(j) 
    log(a/prod(I.cal[fa.list[[j]]]) + sum(compareToRow(Y[, fa.list[[j]], drop = FALSE], yi[fa.list[[j]]])) - member) - 
      log(a/prod(I.cal[pa.list[[j]]]) + sum(compareToRow(Y[, pa.list[[j]], drop = FALSE], yi[pa.list[[j]]])) - member)))
  
  return(p_ik)
  
}


