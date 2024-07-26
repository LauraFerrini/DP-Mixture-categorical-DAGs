library(ggplot2)
library(ggmosaic)
library(dplyr)

source("MCMC/GIBBS_joint_rcpp.R")

X = read.csv("data/breast_cancer.csv")
head(X)
n = nrow(X)
q = ncol(X)
# set the constraints of the adjacency matrix
A_constr = matrix(0,q,q)

colnames(A_constr) = rownames(A_constr) = colnames(X)

A_constr[1,] = NA
A_constr[,c(2,11,12,21)] = NA; A_constr[2,c(11,12,21)] = 0; A_constr[12,11] = 0

risk_factors = c(3,4,8,9,10,15,16,17,18,19,20)
therapies = c(6,7,13,14)

A_constr[therapies, risk_factors] = NA

S = 100000
burn_in = 10000
# pi ~ Beta(a_pi, b_pi)
a_pi = 1
b_pi = 2*q
# Prior hyper-prameters on the concentration param of DP prior 
a_alpha = 3; b_alpha = 1 

# ne = NULL -> no constraints on the maximum number of neighborhoods per node
a = 1
#I.cal  = sapply(1:ncol(X), function(j) length(unique(X[,j]))) # gives the number of levels for each var
X = as.matrix(X)
t0 = proc.time()

set.seed(1234)
out_mcmc = Gibbs_joint(Y = X, S = S, burn_in = burn_in, a_pi = a_pi,
                           b_pi = b_pi, a_alpha = a_alpha, a = a,
                           b_alpha = b_alpha, A_constr = A_constr)
t1 = proc.time() - t0

out_mcmc$DAG
# the output of the mcmc algorithm has been stored in the file "out_breastcancer1.RData"

save(out_mcmc, file = "data/out_breastcancer1.RData")

## recover posterior distribution of causal effects of interest 


######################################
############# PLOTS ##################
######################################
#load("data/out_breastcancer1.RData")

###################################
### Posterior Similarity matrix ###
###################################


library(mcclust.ext)
library(fields)

vi = minVI(out_mcmc$simil_mat)$cl
table(vi)

ordered_ind1_vi = unlist(sapply(1:length(table(vi)),
                                function(k) c(which(vi == k))))
simil_ordered1_vi = out_mcmc$simil_mat[ordered_ind1_vi,ordered_ind1_vi]

n = length(vi)
grid_labs = c(1, round(n/4), n/2, round(3*n/4), n)
pos_labs  = grid_labs/n; pos_labs[1] = 0

pdf(file = "data/psm_breastcancer_def.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7.4) # The height of the plot in inches
colori = colorRampPalette(c('white','black'))
par(mar = c(5,5,1,1), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3.5,1,0), mfrow = c(1,1), mai = c(1,1,0.5,0.5))
image.plot(simil_ordered1_vi, col = colori(100), zlim = c(0,1), cex.sub = 1,
           xlab = "subjects (i)", ylab = "subjects (i')", 
           axes = F, horizontal = F, legend.shrink = 1, border = "black", lwd = 2)
axis(1, at = pos_labs, lab = grid_labs, las = 2)
axis(2, at = pos_labs, lab = grid_labs, las = 2)
abline(h = 1.0025)
abline(v = 1.0025)
abline(h = -0.002)
abline(v = -0.002)
dev.off()


###################################
######### Individual PPI ##########
###################################
library(fields)
graph_probs = out_mcmc$graph_probs
vi    =  minVI(out_mcmc$simil_mat)$cl

# find the two individuals with the smallest posterior probability of being assigned
# to the same clusters 
idx = which(out_mcmc$simil_mat == min(out_mcmc$simil_mat), arr.ind = TRUE)
set.seed(123)
idx = idx[sample(1:nrow(idx), 1), ]
i_cl1 = unname(idx[which(vi[idx]==1)])
i_cl2 = unname(idx[which(vi[idx]==2)])
prob_1 = graph_probs[,,i_cl1]
prob_2 = graph_probs[,,i_cl2]

X = read.csv("data/breast_cancer.csv")
colnames(X)
new_names = unname(sapply(colnames(X), function(x) gsub("_", " ", x)))
colnames(X) = new_names
n = nrow(X)
q = ncol(X)

rownames(prob_1) = colnames(prob_1) = colnames(X)
colnames(prob_2) = rownames(prob_2) = colnames(X)


colori = colorRampPalette(c('white','black'))
pdf(file = "data/PPI_cl1.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 7.2) # The height of the plot in inches
par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(1,1,0.1,1))
image.plot(t(prob_1), col = colori(100), zlim = c(0,1), cex.sub = 1.5, xlab = "", 
           axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 20, cex = 1, cex.axis = 1)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1, cex.axis = 1)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)
dev.off()

pdf(file = "data/PPI_cl2.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 7.2) # The height of the plot in inches
par(oma = c(3,3,1,1.5), mar = c(5,4,0.3,1), mai = c(1,1,0.1,1))
image.plot(t(prob_2), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "", axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2, cex.axis = 1)
axis(2, at = seq(0, 1, l = q), lab = colnames(X), las = 2, srt = 35, cex = 1.2, cex.axis = 1)
abline(h = 1.030)
abline(v = 1.030)
abline(h = -0.03)
abline(v = -0.03)
dev.off()

####################################
########## Spider Plots ############
####################################
library(ggplot2)
library(scales)
library(ggradar)
library(dplyr)

X = read.csv("data/breast_cancer.csv")

X[] <- lapply(X, factor)
str(X)
X_dummy = model.matrix(~., data = X)
X_dummy = X_dummy[,-1]

fun_for_names <- function(my_string) {
  my_string = gsub("_", " ", my_string)
  if (substr(my_string, nchar(my_string), nchar(my_string)) == "1") {
    # Drop the last character
    my_string <- substr(my_string, 1, nchar(my_string) - 1)
  }
  return(my_string)
}

new_names = unname(sapply(colnames(X_dummy), function (x) fun_for_names(x)))
colnames(X_dummy) = new_names

prop = apply(X_dummy, 2, mean)
prop = as.data.frame(t(prop))
dim(prop)

vi      =  minVI(out_mcmc$simil_mat)$cl
sum(vi==1)

X_cl1 = X_dummy[which(vi==1), ]
X_cl2 = X_dummy[which(vi==2), ]
dev.off()
prop_all = rbind(prop, apply(X_cl1, 2, mean), apply(X_cl2, 2, mean))
prop_all = cbind(data = c("pooled", "cluster 1", "cluster 2"), prop_all)

pdf('data/spider_cl1.pdf', pointsize=10, width=9, height=9)
p = prop_all[1:2, ] %>%
  ggradar(group.line.width = 1,
          group.point.size = 1.5, group.colours = c("#1b9e77", "grey")) +
  theme(legend.position = "bottom", legend.title = element_text(size = 17))
p
dev.off()
pdf('data/spider_cl2.pdf', pointsize=10, width=9, height=9)
p = prop_all[c(1,3), ] %>%
  ggradar(group.line.width = 1,
          group.point.size = 1.5,  group.colours = c("#d95f02", "grey")) +
  theme(legend.position = "bottom", legend.title = element_text(size = 17))
p
dev.off()

#################################################
############## Causal effects ###################
#################################################

## Retrieve the causal effects of interest 
source("MCMC/gamma_causal.R")

Xi_burned = out_mcmc$Xi[,(burn_in+1):S]

set.seed(123)
res_ACprev = individual_causal(n = nrow(X), 
                               Xi_chain = Xi_burned,
                               out_mcmc$theta_chain, y = "X1", v = "X13", v_k = 1, v_h = 0)
set.seed(123)
res_AC = individual_causal(n = nrow(X),
                           Xi_chain = Xi_burned,
                           out_mcmc$theta_chain, y = "X1", v = "X6", v_k = 1, v_h = 0)

set.seed(123)
res_antiHER2= individual_causal(n = nrow(X), 
                                Xi_chain = Xi_burned,
                                out_mcmc$theta_chain, y = "X1", v = "X7", v_k = 1, v_h = 0)

out_causal = list(res_AC, res_antiHER2)
save(out_causal, file ="out_causal.RData")

#################################################################
## Causal effects if we would have neglected the heterogeneity ##
#################################################################
source("MCMC/MCMC_pooled.R")
X = read.csv("data/breast_cancer.csv")
n = nrow(X)
q = ncol(X)
# set the constraints of the adjacency matrix as before
A_constr = matrix(0,q,q)

colnames(A_constr) = rownames(A_constr) = colnames(X)

A_constr[1,] = NA
A_constr[,c(2,11,12,21)] = NA; A_constr[2,c(11,12,21)] = 0; A_constr[12,11] = 0

risk_factors = c(3,4,8,9,10,15,16,17,18,19,20)
therapies = c(6,7,13,14)

A_constr[therapies, risk_factors] = NA

S = 100000
burn_in = 10000
# pi ~ Beta(a_pi, b_pi)
a_pi = 1
b_pi = 2*q
# Prior hyper-prameters on the concentration param of DP prior 
a_alpha = 3; b_alpha = 1 

# ne = NULL -> no constraints on the maximum number of neighborhoods per node
a = 1

X = as.matrix(X)
out_pooled = mcmc_pooled(X, S, burn_in, a, a_pi, b_pi, A_constr, joint = TRUE)
AC_pooled = sapply(out_pooled$Theta, function(theta) gammav(theta, y = "X1", v = "X6",
                                                            v_k = 1, v_h =0))
antiHER2_pooled = sapply(out_pooled$Theta, function(theta) gammav(theta, y = "X1", v = "X7",
                                                            v_k = 1, v_h =0))
####################################
##### Plots of causal effects ######
####################################

library(mcclust.ext)
library(ggplot2)
source("out_causal.RData")
str(out_causal)

post_means = lapply(out_causal, function(i) rowMeans(i))
#post_means_nocluster = unlist(lapply(out_causal, function(i) mean(i)))

vi      =  minVI(out_mcmc$simil_mat)$cl

causal_cl1 = lapply(post_means, function (i) i[which(vi == 1)])
causal_cl2 = lapply(post_means, function (i) i[which(vi == 2)])

ggplot_fun = function(causal_cl1, causal_cl2, n, #post_means_nocluster,
                      names_x){
  ggplot(data.frame(post.means = c(causal_cl1, causal_cl2),
                    cluster    = as.factor(c(rep(1, length(causal_cl1)),rep(2, length(causal_cl2)) )),
                    idx = 1:n),
         aes(x= idx, y = post.means, color = cluster)) +
    geom_point() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          legend.text=element_text(size=14),
          legend.title = element_text(size=16),
          plot.title = element_text(size = 16)) +
    scale_color_brewer(palette = "Dark2") +
    labs(title= names_x,
         x = "i", y = "Posterior mean of causal effect") + 
    #geom_hline(yintercept = post_means_nocluster[j], color = "darkblue", lwd = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") 
  
}

par(mfrow=c(2,3))

####################
## Anthracyclines ##
####################

pdf(file = "CE_AC.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 7.5)
ggplot_fun(causal_cl1[[1]], causal_cl2[[1]], nrow(X), "AC")
dev.off()


########################
## AntiHER2 therapies ##
########################

pdf(file = "CE_antiHER2.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 7.5)
ggplot_fun(causal_cl1[[2]], causal_cl2[[2]], nrow(X), "AntiHER2 Thearpies")
dev.off()

