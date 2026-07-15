library(ggplot2)
library(dplyr)
library(bnlearn)
library(gtools)
#library(ggmosaic)
source("GIBBS_collapsed_speedy.R")

X = read.csv("breast_cancer.csv")
head(X)

sapply(1:ncol(X), function(j) table(X[,j]))
q = ncol(X)
n = nrow(X)


# set the constraints of the adjacency matrix
A_constr = matrix(0,q,q)

colnames(A_constr) = rownames(A_constr) = colnames(X)

A_constr[1,] = NA
A_constr[,c(2,11,12,19)] = NA; A_constr[2,c(11,12,19)] = 0; A_constr[12,11] = 0

risk_factors = c(3,4,8,9,10,13:18)
therapies = c(6,7)

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
I.cal  = sapply(1:ncol(X), function(j) length(unique(X[,j]))) # gives the number of levels for each var
X = as.matrix(X)

t0 = proc.time()

set.seed(1234)
out_mcmc = Gibbs_collapsed_speedy(Y = X, S = S, burn_in = burn_in, a_pi = a_pi,
                       b_pi = b_pi, a_alpha = a_alpha, a = a,
                       b_alpha = b_alpha, A_constr = A_constr)
t1 = proc.time() - t0

save(out_mcmc, file = "out_breastcancer.RData")

## Causal effects ------- 
## Sample Theta ---- 
source("theta_function.R")

Dags = out_mcmc$DAG
Xi = out_mcmc$Xi

Theta_all = sample_theta(burn_in = burn_in, S = S, Xi_chain = Xi, DAG_chain = Dags,
                         X = X, q = q, seed = 123)

## Retrieve the causal effects of interest 
source("gamma_causal.R")

Xi_burned = out_mcmc$Xi[,(burn_in + 1):S]
set.seed(123)
res_AC = individual_causal(n = nrow(X),
                           Xi_chain = Xi_burned,
                           Theta_all, y = "X1", v = "X6", v_k = 1, v_h = 0)

set.seed(123)
res_antiHER2= individual_causal(n = nrow(X), 
                                Xi_chain = Xi_burned,
                                Theta_all, y = "X1", v = "X7", v_k = 1, v_h = 0)

out_causal = list(AC = res_AC, antiHER2 = res_antiHER2)
save(out_causal, file ="out_causal.RData")

#################################################################
## Causal effects if we would have neglected the heterogeneity ##
#################################################################

source("MCMC_pooled.R")
X = read.csv("breast_cancer.csv")
n = nrow(X)
q = ncol(X)
# set the constraints of the adjacency matrix as before
A_constr = matrix(0,q,q)

colnames(A_constr) = rownames(A_constr) = colnames(X)

A_constr[1,] = NA
A_constr[,c(2,11,12,19)] = NA; A_constr[2,c(11,12,19)] = 0; A_constr[12,11] = 0

risk_factors = c(3,4,8,9,10,13:18)
therapies = c(6,7)

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

set.seed(123)

out_pooled = mcmc_pooled(Y = X, S = S, burn_in = burn_in,
                         a = a, a_pi = a_pi, b_pi = b_pi, A_constr = A_constr, 
                         joint = TRUE)

out_pooled_burned = out_pooled$Theta[(burn_in+1):S]
AC_pooled = sapply(out_pooled_burned, function(theta) gammav(theta, y = "X1", v = "X6",
                                                             v_k = 1, v_h =0))

antiHER2_pooled = sapply(out_pooled_burned, function(theta) gammav(theta, y = "X1", v = "X7",
                                                                   v_k = 1, v_h =0))

out_causal_pooled = list(AC = AC_pooled, antiHER2 = antiHER2_pooled)
save(out_causal_pooled, file ="out_causal_pooled.RData")


######################################
############# PLOTS ##################
######################################

load("out_breastcancer.RData")

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

pdf(file = "psm.pdf",   # The directory you want to save the file in
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

X = read.csv("breast_cancer.csv")
colnames(X)
new_names = unname(sapply(colnames(X), function(x) gsub("_", " ", x)))
colnames(X) = new_names
n = nrow(X)
q = ncol(X)

rownames(prob_1) = colnames(prob_1) = colnames(X)
colnames(prob_2) = rownames(prob_2) = colnames(X)


colori = colorRampPalette(c('white','black'))
pdf(file = "PPI_cl1.pdf",   # The directory you want to save the file in
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


pdf(file = "PPI_cl2.pdf",   # The directory you want to save the file in
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

X = read.csv("breast_cancer.csv")

X[] = lapply(X, factor)
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

plot_spider_base = function(profile_mat, file_name, main,
                            cols = c("grey40", "#1b9e77")) {
  profile_mat = as.matrix(profile_mat)
  labels = colnames(profile_mat)
  K = ncol(profile_mat)
  angles = seq(0, 2 * pi, length.out = K + 1)[-(K + 1)]
  x = sin(angles)
  y = cos(angles)
  
  pdf(file_name, pointsize = 10, width = 9, height = 9)
  par(mar = c(1, 1, 3, 1))
  plot(NA, xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35),
       axes = FALSE, xlab = "", ylab = "", asp = 1, main = main)
  
  for (r in seq(0.25, 1, by = 0.25)) {
    polygon(r * c(x, x[1]), r * c(y, y[1]), border = "grey85")
    text(0, r, paste0(100 * r, "%"), col = "grey50", cex = 0.75, pos = 3)
  }
  
  segments(0, 0, x, y, col = "grey85")
  text(1.15 * x, 1.15 * y, labels = labels, cex = 0.65)
  
  for (g in seq_len(nrow(profile_mat))) {
    values = as.numeric(profile_mat[g, ])
    lines(c(values * x, values[1] * x[1]),
          c(values * y, values[1] * y[1]),
          col = cols[g], lwd = 2)
    points(values * x, values * y, col = cols[g], pch = 16, cex = 0.6)
  }
  
  legend("bottom", legend = rownames(profile_mat), col = cols,
         lwd = 2, pch = 16, bty = "n", horiz = TRUE)
  dev.off()
}

prop_AC = apply(X_dummy[as.character(X$AC) == "1", , drop = FALSE], 2, mean)
prop_antiHER2 = apply(X_dummy[as.character(X$antiHER2) == "1", , drop = FALSE], 2, mean)

prop_AC_all = rbind(pooled = unlist(prop), `AC = 1` = prop_AC)
plot_spider_base(
  profile_mat = prop_AC_all,
  file_name = "spider_AC_conditioned.pdf",
  main = "Variables profile conditional on AC = 1",
  cols = c("grey40", "#1b9e77")
)

prop_antiHER2_all = rbind(pooled = unlist(prop), `antiHER2 = 1` = prop_antiHER2)
plot_spider_base(
  profile_mat = prop_antiHER2_all,
  file_name = "spider_antiHER2_conditioned.pdf",
  main = "Variables profile conditional on antiHER2 = 1",
  cols = c("grey40", "#d95f02")
)

vi      =  minVI(out_mcmc$simil_mat)$cl
sum(vi==1)

X_cl1 = X_dummy[which(vi==1), ]
X_cl2 = X_dummy[which(vi==2), ]
dev.off()
prop_all = rbind(prop, apply(X_cl1, 2, mean), apply(X_cl2, 2, mean))
prop_all = cbind(data = c("pooled", "cluster 1", "cluster 2"), prop_all)

pdf('spider_cl1.pdf', pointsize=10, width=9, height=9)
p = prop_all[1:2, ] %>%
  ggradar(group.line.width = 1,
          group.point.size = 1.5, group.colours = c("#1b9e77", "grey")) +
  theme(legend.position = "bottom", legend.title = element_text(size = 17)) +
  scale_fill_discrete(labels = c("pooled", expression(hat(c)[1]),  values = c("grey", "#1b9e77")))
p
dev.off()

pdf('spider_cl2.pdf', pointsize=10, width=9, height=9)
p = prop_all[c(1,3), ] %>%
  ggradar(group.line.width = 1,
          group.point.size = 1.5,  group.colours = c("#d95f02", "grey")) +
  theme(legend.position = "bottom", legend.title = element_text(size = 17))
p
dev.off()

#################################################
############## Causal effects ###################
#################################################


####################################
##### Plots of causal effects ######
####################################

library(mcclust.ext)
library(ggplot2)

load("out_causal.RData")
load("out_causal_pooled.RData")
str(out_causal)

post_means = lapply(out_causal, function(i) rowMeans(i))
post_means_pooled = unlist(lapply(out_causal_pooled, function(i) mean(i)))

vi      =  minVI(out_mcmc$simil_mat)$cl

causal_cl1 = lapply(post_means, function (i) i[which(vi == 1)])
causal_cl2 = lapply(post_means, function (i) i[which(vi == 2)])

ggplot_fun = function(causal_cl1, causal_cl2, n, post_means_pooled,
                      names_x, xlab = "subject (i)", ylab = expression(hat(gamma)[i])){
  ggplot(data.frame(post.means = c(causal_cl1, causal_cl2),
                    cluster    = as.factor(c(rep(1, length(causal_cl1)),rep(2, length(causal_cl2)) )),
                    idx = 1:n),
         aes(x= idx, y = post.means, color = cluster)) +
    geom_point() + 
    ylim(-0.013, 0.05) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size = 16, margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16),
          plot.title = element_text(size = 16)) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = names_x, x = xlab, y = ylab, color = expression(hat(c))) + 
    geom_hline(yintercept = post_means_pooled, color = "darkblue", lwd = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") 
}

par(mfrow=c(2,3))

####################
## Anthracyclines ##
####################

pdf(file = "AC.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 7)
ggplot_fun(causal_cl1[[1]], causal_cl2[[1]], nrow(X), post_means_pooled[1], "AC")
dev.off()

########################
## AntiHER2 therapies ##
########################

pdf(file = "CE_antiHER2.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 7)
ggplot_fun(causal_cl1[[2]], causal_cl2[[2]], nrow(X), post_means_pooled[2], "antiHER2")
dev.off()
