## Simulation setting with two clusters
library("parallel")
library(tidyr)
library(ggplot2)
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

load("simulations/q10_nk100_alpha20_b0.RData")
xi_true = c(rep(1, n_k), rep(2, n_k))

vi_1_dags    = sapply(1:N, function(i) Vi(out_simil_dag[[i]], xi_true, n_k, "bnp") )
vi_1_no_dags = sapply(1:N, function(i) Vi(out_simil_nodag[[i]], xi_true, n_k, "bnp") )
vi_1_polca = sapply(1:N, function(i) Vi(pred_cl_polca[[i]], xi_true, n_k, "polca"))
vi_1_kmodes = sapply(1:N, function(i) Vi(pred_cl_kmodes[[i]], xi_true, n_k, "kmodes"))


################################################
###### SCENARIO 2: alpha = 0.2, n_k = 200 ######
################################################

load("simulations/q10_nk200_alpha20_b0.RData")
xi_true = c(rep(1, n_k), rep(2, n_k))

vi_2_dags    = sapply(1:N, function(i) Vi(out_simil_dag[[i]], xi_true, n_k, "bnp") )
vi_2_no_dags = sapply(1:N, function(i) Vi(out_simil_nodag[[i]], xi_true, n_k, "bnp") )
vi_2_polca   = sapply(1:N, function(i) Vi(pred_cl_polca[[i]], xi_true, n_k, "polca"))
vi_2_kmodes  = sapply(1:N, function(i) Vi(pred_cl_kmodes[[i]], xi_true, n_k, "kmodes"))


################################################
###### SCENARIO 3: alpha = 0.2, n_k = 500 ######
################################################
load("simulations/q10_nk500_alpha20_b0.RData")
xi_true = c(rep(1, n_k), rep(2, n_k))

vi_3_dags    = sapply(1:N, function(i) Vi(out_simil_dag[[i]], xi_true, n_k, "bnp") )
vi_3_no_dags = sapply(1:N, function(i) Vi(out_simil_nodag[[i]], xi_true, n_k, "bnp") )
vi_3_polca   = sapply(1:N, function(i) Vi(pred_cl_polca[[i]], xi_true, n_k, "polca"))
vi_3_kmodes  = sapply(1:N, function(i) Vi(pred_cl_kmodes[[i]], xi_true, n_k, "kmodes"))


################################################
###### SCENARIO 4: alpha = 0.4, n_k = 100 ######
################################################
load("simulations/q10_nk100_alpha40_b0.RData")
xi_true = c(rep(1, n_k), rep(2, n_k))

vi_4_dags    = sapply(1:N, function(i) Vi(out_simil_dag[[i]], xi_true, n_k, "bnp") )
vi_4_no_dags = sapply(1:N, function(i) Vi(out_simil_nodag[[i]], xi_true, n_k, "bnp") )
vi_4_polca   = sapply(1:N, function(i) Vi(pred_cl_polca[[i]], xi_true, n_k, "polca"))
vi_4_kmodes  = sapply(1:N, function(i) Vi(pred_cl_kmodes[[i]], xi_true, n_k, "kmodes"))

################################################
###### SCENARIO 5: alpha = 0.4, n_k = 200 ######
################################################

load("simulations/q10_nk200_alpha40_b0.RData")
xi_true = c(rep(1, n_k), rep(2, n_k))

vi_5_dags    = sapply(1:N, function(i) Vi(out_simil_dag[[i]], xi_true, n_k, "bnp") )
vi_5_no_dags = sapply(1:N, function(i) Vi(out_simil_nodag[[i]], xi_true, n_k, "bnp") )
vi_5_polca   = sapply(1:N, function(i) Vi(pred_cl_polca[[i]], xi_true, n_k, "polca"))
vi_5_kmodes  = sapply(1:N, function(i) Vi(pred_cl_kmodes[[i]], xi_true, n_k, "kmodes"))


################################################
###### SCENARIO 6: alpha = 0.4, n_k = 500 ######
################################################
load("simulations/q10_nk500_alpha40_b0.RData")
xi_true = c(rep(1, n_k), rep(2, n_k))

vi_6_dags    = sapply(1:N, function(i) Vi(out_simil_dag[[i]], xi_true, n_k, "bnp") )
vi_6_no_dags = sapply(1:N, function(i) Vi(out_simil_nodag[[i]], xi_true, n_k, "bnp") )
vi_6_polca   = sapply(1:N, function(i) Vi(pred_cl_polca[[i]], xi_true, n_k, "polca"))
vi_6_kmodes  = sapply(1:N, function(i) Vi(pred_cl_kmodes[[i]], xi_true, n_k, "kmodes"))

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
#### Plot alpha = 0.2 ####
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
k = 1


for(r in 1:rep){
  for (l in 1:N){
    
    A_true1 = dag2cpdag(as(out3[[r]][[l]]$data$Y1_data$dag, "graphNEL" ))
    A_true2 = dag2cpdag(as(out3[[r]][[l]]$data$Y2_data$dag, "graphNEL" ))
    
    A_est1 = round(out3[[r]][[l]]$out_probs_dag_oracle[,,1])
    A_est2 = round(out3[[r]][[l]]$out_probs_dag_oracle[,,2])
    
    colnames(A_est1) = rownames(A_est1) = 1:q
    colnames(A_est2) = rownames(A_est2) = 1:q
    A_est1 = dag2cpdag(as(A_est1, "graphNEL") )
    A_est2 = dag2cpdag(as(A_est2, "graphNEL") )
    
    shd3_oracle[k, 1]   = pcalg::shd(A_est1, A_true1 )
    shd3_oracle[k, 2] = pcalg::shd(A_est2, A_true2 )
    k = k + 1
  }
}

################################################
###### SCENARIO 1: alpha = 0.2, n_k = 100 ######
################################################


################################################
###### SCENARIO 2: alpha = 0.2, n_k = 200 ######
################################################


################################################
###### SCENARIO 3: alpha = 0.2, n_k = 500 ######
################################################


################################################
###### SCENARIO 4: alpha = 0.4, n_k = 100 ######
################################################


################################################
###### SCENARIO 5: alpha = 0.4, n_k = 200 ######
################################################



################################################
###### SCENARIO 6: alpha = 0.4, n_k = 500 ######
################################################


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
