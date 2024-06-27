library(mcclust.ext)
library(ggplot2)
library(scales)
library(devtools)
library(ggradar)
library(dplyr)

load("data/out_breastcancer1.RData")
X = read.csv("data/data.csv")

X = X[,-1]
X = X[, c(3,2,1,4:21)]
colnames(X)

A = matrix(unlist(lapply(1:ncol(X), function(j) as.factor(X[,j]))), nrow = nrow(X), ncol=ncol(X), byrow = F)
head(A);head(X)
colnames(A) = colnames(X)
X[] <- lapply( X, factor)
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
prop_all[1:2, ] %>%
  ggradar(group.line.width = 1,
          group.point.size = 1.5, group.colours = c("#1b9e77", "grey")) +
  theme(legend.position = "bottom", legend.title = element_text(size = 17))

dev.off()
pdf('data/spider_cl2.pdf', pointsize=10, width=9, height=9)
prop_all[c(1,3), ] %>%
  ggradar(group.line.width = 1,
          group.point.size = 1.5,  group.colours = c("#d95f02", "grey")) +
  theme(legend.position = "bottom", legend.title = element_text(size = 17))
dev.off()


