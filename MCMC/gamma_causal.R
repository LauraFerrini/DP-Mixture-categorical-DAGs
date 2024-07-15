library(bnlearn)
gammav = function(theta, y, v, v_k, v_h) {
  
  # e.g. v_k = 1 and v_h = 0 for E(Y | X = v_k) - E(Y | X = v_h)
  
  # This function computes the causal effect of do(Xv = x) on y
  # from the DAG-parameter theta (a collection of node-conditional probabilities)
  
  res = 0
  pa = theta[[v]]$parents
  eventStr = sprintf("(%s == 1)", y)
  
  if (length(pa)== 0) {
    cp1 = sprintf("cpquery(theta, event = %s, evidence = ((%s == %s)))", eventStr, v, v_k)
    cp2 = sprintf("cpquery(theta, event = %s, evidence = ((%s == %s)))", eventStr, v, v_h)
    res = eval(parse(text=cp1)) - eval(parse(text=cp2))
  } else { # 2 or more parents
    xxx = theta[[v]]$prob
    tab = apply(xxx, seq(2, length(dim(xxx))), sum)
    if (length(pa)==1) {
      df = data.frame(pa = names(tab))
      paname = pa
    } else {
      df = as.data.frame.table(tab)[,1:length(pa)]
      #paname = names(df[i,1:ncol(df)])
      paname = names(df[,1:ncol(df)])
    }
    for (i in 1:nrow(df)) {
      evidenceStr = ""
      for (j in 1:ncol(df)) {
        val = df[i, j]
        n = paname[j]
        if (evidenceStr == "") {
          evidenceStr = sprintf("(%s == %s)", n, val)
        } else {
          evidenceStr = sprintf("%s & (%s == %s)", evidenceStr, n, val)
        }
      }
      cp1 = sprintf("cpquery(theta, event = %s, evidence = ((%s == %s) & %s))", eventStr, v, v_k, evidenceStr)
      cp2 = sprintf("cpquery(theta, event = %s, evidence = ((%s == %s) & %s))", eventStr, v, v_h, evidenceStr)
      cp3 = sprintf("cpquery(theta, event = (%s), evidence = TRUE)", evidenceStr)
      res = res + (eval(parse(text=cp1)) - eval(parse(text=cp2))) * eval(parse(text=cp3))
    }
  }
  return(res)
}


individual_causal <- function(n, S, burnin, Xi_chain, Theta_chain, y, v, v_k, v_h) {
  
  res_causal = matrix(nrow = n, ncol= (S-burnin))
  r = 1
  for (s in (burnin + 1):S){
    
    xis = Xi_chain[,s]
    causal = sapply(1:length(unique(xis)), function(k) gammav(Theta_chain[[s]][k][[1]],
                                                              y = y, v = v, v_k = v_k, v_h = v_h) )
    res_causal[,r] = causal[xis]
    r = r + 1
    
    if(s%%100 == 0){print(s)}
    
  }
  return(res_causal)
}

