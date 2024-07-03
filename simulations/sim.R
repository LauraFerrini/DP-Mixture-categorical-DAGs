## Simulation setting with two clusters

rm(list = ls())

# memory.size(500000)

######

q = 10
w = 3/(2*q)

alpha = 0.4

n_k = 100
#n_k = 200
#n_k = 500

S    = 12000
burn = 2000

N = 40 # changed as 40 the number of replications

######

library("parallel")

maxcore = min(detectCores(), N)

maxcore = 10

tm   = proc.time()
out  = source("simulations/parallel_sim.R")
time = proc.time() - tm


filesaved = paste0("simulations/q", q, "_nk", n_k, "_alpha", alpha*100,".RData")
save.image(file = filesaved)


##################################################
#### this can be deleted ????????


#########
## OLD ## 
#########

str(out)

n = n_k*2

library(pcalg)

shd(dag2essgraph(as(round(out$value[[2]]$graph_probs[,,n]), "graphNEL")), dag2essgraph(as(out$value[[2]]$dags[,,2], "graphNEL")))
shd(dag2essgraph(as(round(out$value[[2]]$graph_probs[,,1]), "graphNEL")), dag2essgraph(as(out$value[[2]]$dags[,,1], "graphNEL")))


shd_to_EG = function(j){
  
  shds = c()
  
  for(i in 1:n){
    
    dag_est  = as(round(out$value[[j]]$graph_probs[,,i]), "graphNEL")
    dag_true = as(out$value[[j]]$dags[,,out$value[[j]]$xi[i]], "graphNEL")
    
    shds[i] = shd(dag2essgraph(dag_est),dag2essgraph(dag_true))
    
  }
  
  return(shds)
  
}

n = n_k*2

shds_1 = sapply(1:N, function(j) shd_to_EG(j))

shds_2 = sapply(1:N, function(j) shd_to_EG(j))

shds_3 = sapply(1:N, function(j) shd_to_EG(j))

shds_4 = sapply(1:N, function(j) shd_to_EG(j))

par(mfrow = c(1,4))
boxplot(c(shds_4), ylim = c(0,30), outline = FALSE)
boxplot(c(shds_3), ylim = c(0,30), outline = FALSE)
boxplot(c(shds_2), ylim = c(0,30), outline = FALSE)
boxplot(c(shds_1), ylim = c(0,30), outline = FALSE)

n = n_k*K

shds_5 = sapply(1:N, function(j) shd_to_EG(j))

shds_6 = sapply(1:N, function(j) shd_to_EG(j))

shds_7 = sapply(1:N, function(j) shd_to_EG(j))

shds_8 = sapply(1:N, function(j) shd_to_EG(j))

par(mfrow = c(1,4))
boxplot(c(shds_5), ylim = c(0,30), outline = FALSE)
boxplot(c(shds_6), ylim = c(0,30), outline = FALSE)
boxplot(c(shds_7), ylim = c(0,30), outline = FALSE)
boxplot(c(shds_8), ylim = c(0,30), outline = FALSE)
