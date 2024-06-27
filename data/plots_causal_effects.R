load("out_causal_effects-2.RData")
load("S100K_round1.RData")
load("out_causal_nocluster-2.RData")

library(mcclust.ext)
library(ggplot2)

str(out_causal)
str(out_causal_nocluster)
post_means = lapply(out_causal, function(i) rowMeans(i))
post_means_nocluster = unlist(lapply(out_causal, function(i) mean(i)))

int_var = names(out_causal)
vi      =  minVI(out_mcmc$simil_mat)$cl
causal_cl1 = lapply(post_means, function (i) i[which(vi == 1)])

causal_cl2 = lapply(post_means, function (i) i[which(vi == 2)])


plot_fun = function(j, causal_cl1, causal_cl2, post_means_nocluster){
  plot(c(causal_cl1[[j]], causal_cl2[[j]]),
       col = c(rep("salmon", length(causal_cl1[[1]])),
               rep("dodgerblue", length(causal_cl2[[1]]))),
       ylab = paste("Posterior means of causal effects of ", 
                    int_var[j], " on CTRCD"), 
       ylim = c(-0.1, 0.1))
  abline(h = post_means_nocluster[j], col = "darkblue", lwd = 2)
}
int_var = unname(sapply(int_var, function(x) gsub("_", " ", x)))
nomi_estesi = c("Previous antracycline treatment", "Previous anti-HER2 therapies",
                "Hypertension", "Diabetes mellitus", "Left Ventricular Ejection Fraction", 
                "Left Ventricular Ejection Fraction", "Antracycline", "antiHER2")
par(mfrow=c(2,3))
lapply(1:length(causal_cl1), function (j) plot_fun(j, causal_cl1, causal_cl2, post_means_nocluster))
ggplot_fun = function(j, causal_cl1, causal_cl2, post_means_nocluster, names_x){
  ggplot(data.frame(post.means = c(causal_cl1[[j]], causal_cl2[[j]]),
                    cluster    = as.factor(c(rep(1, length(causal_cl1[[j]])),rep(2, length(causal_cl2[[j]])) )),
                    idx = 1:length(vi)),
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
    labs(title=paste(nomi_estesi[j]," - ", names_x[j]),
         x = "i", y = "Posterior mean of causal effect") + 
    geom_hline(yintercept = post_means_nocluster[j], color = "darkblue", lwd = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") 
    
}


pdf(file = "CE_AntiHER2.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 7.5)
# The height of the plot in inches

ggplot_fun(8, causal_cl1, causal_cl2, post_means_nocluster, int_var)
dev.off()
lapply(1:length(causal_cl1), function (j) ggplot_fun(j, causal_cl1, causal_cl2, post_means_nocluster, int_var))

############

post_means_ordered = rbind(matrix(unlist(causal_cl1), 
                                  ncol =  length(causal_cl1), 
                                  nrow = length(causal_cl1[[1]]), byrow = F),
                           matrix(unlist(causal_cl2), 
                                  ncol =  length(causal_cl2), 
                                  nrow = length(causal_cl2[[1]]), byrow = F ))
str(post_means_ordered)
mat_df = as.data.frame(post_means_ordered)
dim(mat_df)
names(mat_df) = int_var

# Aggiungi un'indicizzazione delle righe
mat_df$rw = 1:nrow(mat_df)
mat_df = reshape2::melt(mat_df, id.vars = "rw")
head(mat_df)

pl <- ggplot(mat_df, aes(x = variable, y = rw, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colors = c("brown2", "white", "dodgerblue"),
                       values = scales::rescale(c(-1, 0, 1))) +
  theme_minimal() +
  labs(title = "Posterior means of subject-specific causal effects") +
  xlab("variables") +
  ylab("individuals")
# 0 bianco, rosso neg, verde positivo 

print(pl)
dev.off()


## boxplots: non si vede la differenza per niente 
ind_cl1 = sample(1:length(causal_cl1[[1]]), round(length(causal_cl1[[1]])*0.05))
ind_cl2 = sample(1:length(causal_cl1[[2]]), round(length(causal_cl2[[1]])*0.05))

df = data.frame(values = as.vector(out_causal[[1]][c(ind_cl1, ind_cl2), ]))
df$id = as.factor(rep(1:(length(ind_cl1) + length(ind_cl2)), each = 90000))
head(df)
mean_values = aggregate(values ~ id, data = df, FUN = mean)
dim(mean_values)
dim(df)
df$cl = c(rep("A", length(ind_cl1)*90000),rep("B", length(ind_cl2)*90000) )
# with medians 
p = ggplot(df, aes(x=id, y=values)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-0.1, 0.1))
p

# with means
p = ggplot(df, aes(x = id, y = values, fill = cl)) +
  geom_boxplot(fatten = NULL, outlier.shape = NA) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               size = 1, linetype = "solid") +
  scale_fill_manual(values = c("A" = "gold1", "B" = "seagreen2")) +
  coord_cartesian(ylim = c(-0.1, 0.1))


########################
p.pos <- ggplot(data.frame(pos = positive.probs, clusters = as.factor(vi)), aes(x=clusters, y=pos)) + 
  geom_boxplot() + ylab("Posterior prob of positive causal effects") 
p.pos

p.neg <- ggplot(data.frame(neg = negative.probs, clusters = as.factor(vi)), aes(x=clusters, y=neg)) + 
  geom_boxplot() + ylab("Posterior prob of negative causal effects") 
p.neg

p.null <- ggplot(data.frame(null = null.probs, clusters = as.factor(vi)), aes(x=clusters, y=null)) + 
  geom_boxplot() + ylab("Posterior prob of null causal effects")
p.null
###############
cols = c("CTRCD", "age", "heart_rate", "heart_rhythm", "LVEF", "AC", "antiHER2",    
         "HTA","DL","DM","smoker","exsmoker","ACprev","antiHER2prev",
         "RTprev","CIprev","ICMprev","ARRprev","VALVprev","cxvalv", "BMI")
colnames(X) = cols
p.age = ggplot(data.frame(age = as.factor(X$age), causal_effects = post_means), aes(x=age, y=causal_effects)) + 
  geom_boxplot() + ylab("Causal effect") 
p.age # strange behaviour here


p.heart_rythm = ggplot(data.frame(heart_rythm = as.factor(X$heart_rhythm),
                                  causal_effects = post_means), aes(x=heart_rythm, y=causal_effects)) + 
  geom_boxplot() + ylab("Causal effect") 
p.heart_rythm

p.lvef = ggplot(data.frame(lvef= as.factor(X$LVEF),
                           causal_effects = post_means), aes(x=lvef, y=causal_effects)) + 
  geom_boxplot() + ylab("Causal effect") 
p.lvef

p.antiHER2 = ggplot(data.frame(antiHER2 = as.factor(X$antiHER2),
                               causal_effects = post_means), aes(x=antiHER2, y=causal_effects)) + 
  geom_boxplot() + ylab("Causal effect") 
p.antiHER2

p.HTA = ggplot(data.frame(HTA = as.factor(X$HTA),
                          causal_effects = post_means), aes(x=HTA, y=causal_effects)) + 
  geom_boxplot() + ylab("Causal effect") 
p.HTA
