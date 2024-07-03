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


# Y = as.matrix(Y)
# str(Y)
# 
# xi = c(rep(1, nrow(Y)/2), rep(2, nrow(Y)/2))
# # k = 1
# # i = 345
# # Y_k = Y[xi ==1, ]
# # yi = Y[i, ]
# # str(Y_k)
# # D1 = adj.matr1
# fa.list = lapply(1:q, function(j) fa(j, adj.matr1))
# pa.list = lapply(1:q, function(j) pa(j, adj.matr1))
# # 
# i = 1
# k = 1
# member = xi[i] == k
# j = 1
# 
# n_i_faj_list = sapply(1:q, function(j) sum(compareToRow(Y_k[, fa.k.list[[j]], drop = FALSE], yi[fa.k.list[[j]]])))
# 
# n_i_paj_list = sapply(1:q, function(j) sum(compareToRow(Y_k[, pa.k.list[[j]], drop = FALSE], yi[pa.k.list[[j]]])))
# 
# sapply(1:q, function(j) sum(apply(Y_k[, pa.k.list[[j]], drop = F], 1, function(i) all(i == yi[pa.k.list[[j]]]))) )
# 
# p_ik = log(dim(Y_k)[1] - member) + sum(sapply(1:q, function(j)
#   log(a/prod(I.cal[fa.k.list[[j]]]) + sum(compareToRow(Y_k[, fa.k.list[[j]], drop = FALSE], yi[fa.k.list[[j]]])) - member) -
#     log(a/prod(I.cal[pa.k.list[[j]]]) +sum(compareToRow(Y_k[, pa.k.list[[j]], drop = FALSE], yi[pa.k.list[[j]]])) - member)))
# 
# # check che venga uguale a prima: giusto
# 
# 
# colnames(Y_k) = as.character(paste0("X",1:q))
# N = plyr::count(Y_k)
# str(N)
# Nk.list.fa = lapply(1:q, function(j) aggregate(N$freq, by = as.list(N[fa.k.list[[j]]]), FUN = sum))
# 
# set_pa = which(lapply(1:q, function(j) length(pa.k.list[[j]])) > 0)
# 
# Nk.list.pa = vector("list", length = q) # anche al di fuori di tutto il codice; lo sovrascriviamo ogni volta
# Nk.list.pa[set_pa] = lapply(set_pa, function(j) aggregate(N$freq, by = as.list(N[pa.k.list[[j]]]), FUN = sum))
# 
# 
# prob_ik_nonempty(N, Nk.list.fa, Nk.list.pa, pa.k.list, fa.k.list, TRUE, yi)
