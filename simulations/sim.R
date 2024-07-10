## Simulation setting with two clusters
library("parallel")
library(tidyr)
library(ggplot2)
library(mcclust.ext)
rm(list = ls())

# memory.size(500000)

######

q = 10      # number of nodes
w = 3/(2*q) # probability of edges 
S    = 12000
burn = 2000

N = 40  # number of replications

#maxcore = min(detectCores(), N)
maxcore = 10

################################################
###### SCENARIO 1: alpha = 0.2, n_k = 100 ######
################################################

alpha = 0.2
n_k = 100
out  = source("simulations/parallel_sim.R")

filesaved = paste0("simulations/q", q, "_nk", n_k, "_alpha", alpha*100,".RData")
save.image(file = filesaved)

################################################
###### SCENARIO 2: alpha = 0.2, n_k = 200 ######
################################################
n_k = 200
out  = source("simulations/parallel_sim.R")
filesaved = paste0("simulations/q", q, "_nk", n_k, "_alpha", alpha*100,".RData")
save.image(file = filesaved)

################################################
###### SCENARIO 3: alpha = 0.2, n_k = 500 ######
################################################

n_k = 500
out  = source("simulations/parallel_sim.R")
filesaved = paste0("simulations/q", q, "_nk", n_k, "_alpha", alpha*100,".RData")
save.image(file = filesaved)

################################################
###### SCENARIO 4: alpha = 0.4, n_k = 100 ######
################################################
alpha = 0.4
n_k = 100
out  = source("simulations/parallel_sim.R")
filesaved = paste0("simulations/q", q, "_nk", n_k, "_alpha", alpha*100,".RData")
save.image(file = filesaved)

################################################
###### SCENARIO 5: alpha = 0.4, n_k = 200 ######
################################################
n_k = 200
out  = source("simulations/parallel_sim.R")
filesaved = paste0("simulations/q", q, "_nk", n_k, "_alpha", alpha*100,".RData")
save.image(file = filesaved)

################################################
###### SCENARIO 4: alpha = 0.4, n_k = 500 ######
################################################
n_k = 500
out  = source("simulations/parallel_sim.R")
filesaved = paste0("simulations/q", q, "_nk", n_k, "_alpha", alpha*100,".RData")
save.image(file = filesaved)

########################################################################################
############## PLOT: comparison of the retrieved clustering structure ##################
########################################################################################

Vi = function(simil_mat, xi_true, n_k, method){
  if (method == "bnp"){
    min_VI  =  minVI(simil_mat)$cl
    vi = (vi.dist(min_VI,   xi_true))/log(2*n_k)
  }
  else{
    vi = vi.dist(simil_mat, xi_true)/ log(2*n_k) 
  }
  return(vi)
}

################################################
###### SCENARIO 1: alpha = 0.2, n_k = 100 ######
################################################

load("simulations/q10_nk100_alpha20.RData")
out = out$value

vi_1_dags    = sapply(1:N, function(i) Vi(out[[i]]$out_simil_dag, xi, n_k, "bnp") )
vi_1_no_dags = sapply(1:N, function(i) Vi(out[[i]]$out_simil_nodag, xi, n_k, "bnp") )
vi_1_polca = sapply(1:N, function(i) Vi(out[[i]]$out_polca, xi, n_k, "polca"))
vi_1_kmodes = sapply(1:N, function(i) Vi(out[[i]]$out_kmodes, xi, n_k, "kmodes"))


################################################
###### SCENARIO 2: alpha = 0.2, n_k = 200 ######
################################################

load("simulations/q10_nk200_alpha20.RData")
out = out$value

vi_2_dags    = sapply(1:N, function(i) Vi(out[[i]]$out_simil_dag, xi, n_k, "bnp") )
vi_2_no_dags = sapply(1:N, function(i) Vi(out[[i]]$out_simil_nodag, xi, n_k, "bnp") )
vi_2_polca   = sapply(1:N, function(i) Vi(out[[i]]$out_polca, xi, n_k, "polca"))
vi_2_kmodes  = sapply(1:N, function(i) Vi(out[[i]]$out_kmodes, xi, n_k, "kmodes"))

################################################
###### SCENARIO 3: alpha = 0.2, n_k = 500 ######
################################################
load("simulations/q10_nk500_alpha20.RData")
out = out$value

vi_3_dags    = sapply(1:N, function(i) Vi(out[[i]]$out_simil_dag, xi, n_k, "bnp") )
vi_3_no_dags = sapply(1:N, function(i) Vi(out[[i]]$out_simil_nodag, xi, n_k, "bnp") )
vi_3_polca   = sapply(1:N, function(i) Vi(out[[i]]$out_polca, xi, n_k, "polca"))
vi_3_kmodes  = sapply(1:N, function(i) Vi(out[[i]]$out_kmodes, xi, n_k, "kmodes"))

################################################
###### SCENARIO 4: alpha = 0.4, n_k = 100 ######
################################################
load("simulations/q10_nk100_alpha40.RData")
out = out$value

vi_4_dags    = sapply(1:N, function(i) Vi(out[[i]]$out_simil_dag, xi, n_k, "bnp") )
vi_4_no_dags = sapply(1:N, function(i) Vi(out[[i]]$out_simil_nodag, xi, n_k, "bnp") )
vi_4_polca   = sapply(1:N, function(i) Vi(out[[i]]$out_polca, xi, n_k, "polca"))
vi_4_kmodes  = sapply(1:N, function(i) Vi(out[[i]]$out_kmodes, xi, n_k, "kmodes"))

################################################
###### SCENARIO 5: alpha = 0.4, n_k = 200 ######
################################################

load("simulations/q10_nk200_alpha40.RData")
out = out$value

vi_5_dags    = sapply(1:N, function(i) Vi(out[[i]]$out_simil_dag, xi, n_k, "bnp") )
vi_5_no_dags = sapply(1:N, function(i) Vi(out[[i]]$out_simil_nodag, xi, n_k, "bnp") )
vi_5_polca   = sapply(1:N, function(i) Vi(out[[i]]$out_polca, xi, n_k, "polca"))
vi_5_kmodes  = sapply(1:N, function(i) Vi(out[[i]]$out_kmodes, xi, n_k, "kmodes"))

################################################
###### SCENARIO 6: alpha = 0.4, n_k = 500 ######
################################################
load("simulations/q10_nk500_alpha40.RData")
out = out$value

vi_6_dags    = sapply(1:N, function(i) Vi(out[[i]]$out_simil_dag, xi, n_k, "bnp") )
vi_6_no_dags = sapply(1:N, function(i) Vi(out[[i]]$out_simil_nodag, xi, n_k, "bnp") )
vi_6_polca   = sapply(1:N, function(i) Vi(out[[i]]$out_polca, xi, n_k, "polca"))
vi_6_kmodes  = sapply(1:N, function(i) Vi(out[[i]]$out_kmodes, xi, n_k, "kmodes"))

##########################
#### Plot alpha = 0.2 ####
##########################
VI_all_alpha20 = data.frame(nk_100_dags = vi_1_dags,
                    nk_100_no_dag = vi_1_no_dags,
                    nk_100_polca = vi_1_polca,
                    nk_100_kmodes = vi_1_kmodes,
                    nk_200_dags = vi_2_dags,
                    nk_200_no_dag = vi_2_no_dags,
                    nk_200_polca = vi_2_polca,
                    nk_200_kmodes = vi_2_kmodes,
                    nk_500_dags = vi_3_dags,
                    nk_500_no_dag = vi_3_no_dags,
                    nk_500_polca = vi_3_polca,
                    nk_500_kmodes = vi_3_kmodes)


df_long = gather(VI_all_alpha20, key = "variable", value = "value")
cols = rep(c(rep("steelblue3", N) , rep("salmon2", N), rep("gold2", N),
             rep("palegreen2", N)), 3)
cols = as.factor(cols)
df_long$cols = cols
x_breaks = c("nk_100_no_dag", "nk_200_no_dag", "nk_500_no_dag")  # Positions where ticks will be placed
x_labels = c(100, 200, 500)  # Corresponding labels for the ticks

pdf(file = "simulations/vi_alpha20.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6.5) # The height of the plot in inches
ggplot(df_long, aes(x = variable, y = value, fill = cols)) +
  geom_boxplot() + 
  scale_fill_manual(values = levels(df_long$cols), labels = c("LCA", "K-modes", "No DAG", "DAG mixture"))+
  labs(x = expression(n[k]), y = "Variation of Information", title = expression(alpha ~ "= 0.2"), fill ="Method") +
  scale_x_discrete(breaks = x_breaks, labels =x_labels) +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.text=element_text(size=14),
        legend.title = element_text(size=16),
        plot.title = element_text(size = 16))
dev.off()


##########################
#### Plot alpha = 0.4 ####
##########################

VI_all_alpha40 = data.frame(nk_100_dags = vi_4_dags,
                            nk_100_no_dag = vi_4_no_dags,
                            nk_100_polca = vi_4_polca,
                            nk_100_kmodes = vi_4_kmodes,
                            nk_200_dags = vi_5_dags,
                            nk_200_no_dag = vi_5_no_dags,
                            nk_200_polca = vi_5_polca,
                            nk_200_kmodes = vi_5_kmodes,
                            nk_500_dags = vi_6_dags,
                            nk_500_no_dag = vi_6_no_dags,
                            nk_500_polca = vi_6_polca,
                            nk_500_kmodes = vi_6_kmodes)


df_long = gather(VI_all_alpha40, key = "variable", value = "value")
cols = rep(c(rep("steelblue3", N) , rep("salmon2", N), rep("gold2", N),
             rep("palegreen2", N)), 3)
cols = as.factor(cols)
df_long$cols = cols
x_breaks = c("nk_100_no_dag", "nk_200_no_dag", "nk_500_no_dag")  # Positions where ticks will be placed
x_labels = c(100, 200, 500)  # Corresponding labels for the ticks

pdf(file = "simulations/vi_alpha40.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6.5) # The height of the plot in inches
ggplot(df_long, aes(x = variable, y = value, fill = cols)) +
  geom_boxplot() + 
  scale_fill_manual(values = levels(df_long$cols), labels = c("LCA", "K-modes", "No DAG", "DAG mixture"))+
  labs(x = expression(n[k]), y = "Variation of Information", title = expression(alpha ~ "= 0.2"), fill ="Method") +
  scale_x_discrete(breaks = x_breaks, labels =x_labels) +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        plot.title = element_text(size = 16))
dev.off()



#################################################################################
############## PLOT: comparison of the retrieved DAG structure ##################
#################################################################################

shd_oracle_fun <- function(N, out, q) {
  
  shd_oracle =  matrix(0, nrow = N, ncol = 2)
  
  for(r in 1:N){
    
    A_true1 = dag2cpdag(as(out[[r]]$data$dag_1, "graphNEL" ))
    A_true2 = dag2cpdag(as(out[[r]]$data$dag_2, "graphNEL" ))
      
    A_est1 = round(out[[r]]$out_probs_dag_oracle[,,1])
    A_est2 = round(out[[r]]$out_probs_dag_oracle[,,2])
      
    colnames(A_est1) = rownames(A_est1) = 1:q
    colnames(A_est2) = rownames(A_est2) = 1:q
    A_est1 = dag2cpdag(as(A_est1, "graphNEL") )
    A_est2 = dag2cpdag(as(A_est2, "graphNEL") )
      
    shd_oracle[r, 1]   = pcalg::shd(A_est1, A_true1 )
    shd_oracle[r, 2] = pcalg::shd(A_est2, A_true2 )
    
  }
  return(shd_oracle)
}


shd_function <- function(arr_true, q, xi_i, a_est) {
  
  colnames(a_est) = rownames(a_est) = 1:q
  return(pcalg::shd(dag2cpdag(as(a_est, "graphNEL")), dag2cpdag(as(arr_true[,,xi_i], "graphNEL"))) )
}

shd_mean_function <- function(n_k, N, out, q, xi) {
  
  shd_individual = matrix(0, nrow = 2*n_k, ncol = N)
  idx = 1
  for (r in 1:N){
    shd_individual[, r] = sapply(1:(2*n_k), function(i) shd_function(abind(out[[r]]$data$dag_1,
                                                                            out[[r]]$data$dag_2, along = 3),
                                                                           q, xi[i], round(out[[r]]$out_probs_dag[,,i])) )
    }
  
  shd_means = sapply(1:N, function(r) c(mean(shd_individual[1:n_k,r]), mean(shd_individual[(n_k+1):(2*n_k),r])) )
  return(shd_means)
}

################################################
###### SCENARIO 1: alpha = 0.2, n_k = 100 ######
################################################
load("simulations/q10_nk100_alpha20.RData")
out = out$value

shd1_oracle = shd_oracle_fun(N = N, out = out, q = q)
shd1 = t(shd_mean_function(n_k, N, out, q, xi))
################################################
###### SCENARIO 2: alpha = 0.2, n_k = 200 ######
################################################
load("simulations/q10_nk200_alpha20.RData")
out = out$value

shd2_oracle = shd_oracle_fun(N = N, out = out, q = q)
shd2 = t(shd_mean_function(n_k, N, out, q, xi))
################################################
###### SCENARIO 3: alpha = 0.2, n_k = 500 ######
################################################
load("simulations/q10_nk500_alpha20.RData")
out = out$value

shd3_oracle = shd_oracle_fun(N = N, out = out, q = q)
shd3 = t(shd_mean_function(n_k, N, out, q, xi))

################################################
###### SCENARIO 4: alpha = 0.4, n_k = 100 ######
################################################
load("simulations/q10_nk100_alpha40.RData")
out = out$value

shd4_oracle = shd_oracle_fun(N = N, out = out, q = q)
shd4 = t(shd_mean_function(n_k, N, out, q, xi))
################################################
###### SCENARIO 5: alpha = 0.4, n_k = 200 ######
################################################
load("simulations/q10_nk200_alpha40.RData")
out = out$value

shd5_oracle = shd_oracle_fun(N = N, out = out, q = q)
shd5 = t(shd_mean_function(n_k, N, out, q, xi))

################################################
###### SCENARIO 6: alpha = 0.4, n_k = 500 ######
################################################
load("simulations/q10_nk500_alpha40.RData")

out = out$value

shd6_oracle = shd_oracle_fun(N = N, out = out, q = q)
shd6 = t(shd_mean_function(n_k, N, out, q, xi))

#########################
### PLOT: alpha = 0.2 ###
#########################

shd_alpha20 = data.frame( shd_oracle_nk100 = as.vector(shd1_oracle),
                           shd_nk100 = as.vector(shd1),
                           shd_oracle_nk200 = as.vector(shd2_oracle),
                           shd_nk200 = as.vector(shd2),
                           shd_oracle_nk500 = as.vector(shd3_oracle),
                           shd_nk500 = as.vector(shd3) )
png("alpha20_shd.png", width = 600, height = 400)    
boxplot( shd_alpha20, col = rep(c("lightblue",
                                   "khaki"),3), ylab = "SHD", xaxt="n", main = expression(alpha ~ "= 0.2"),
         ylim = c(0, 15))

axis(1,at=1.3+2*(0:2),labels=c(100,200, 500),tick=FALSE)
axis(1,at=4,labels= expression(n[k]), tick = F, line = 2)


legend("topright", legend = c("Oracle", "DP mixture"),
       col = c("lightblue",
               "khaki"), lty = 1, cex = 0.8)
dev.off()

#########################
### PLOT: alpha = 0.4 ###
#########################

shd_alpha40 = data.frame( shd_oracle_nk100 = as.vector(shd4_oracle),
                          shd_nk100 = as.vector(shd4),
                          shd_oracle_nk200 = as.vector(shd5_oracle),
                          shd_nk200 = as.vector(shd5),
                          shd_oracle_nk500 = as.vector(shd6_oracle),
                          shd_nk500 = as.vector(shd6) )
png("alpha40_shd.png", width = 600, height = 400)    
boxplot( shd_alpha40, col = rep(c("lightblue",
                                  "khaki"),3), ylab = "SHD", xaxt="n", main = expression(alpha ~ "= 0.2"),
         ylim = c(0, 15))

axis(1,at=1.3+2*(0:2),labels=c(100,200, 500),tick=FALSE)
axis(1,at=4,labels= expression(n[k]), tick = F, line = 2)

legend("topright", legend = c("Oracle", "DP mixture"),
       col = c("lightblue",
               "khaki"), lty = 1, cex = 0.8)
dev.off()

