library(ggplot2)
library(MASS)
library(dplyr)
library(umap)
library(viridis)
library(gridExtra)
library(reticulate)

use_python("/N/soft/rhel8/python/gnu/3.10.5/bin/python")
py_config()
py_available()

# source("../../../Algorithms/DR metrics.R")
source("Algorithms/DR metrics.R")

# data = read.delim("../../../../../CyTOF data/Exp1_NK_uns_CD4_inf.txt")
data = read.delim("Data/Exp1_NK_uns_CD4_inf.txt")
data = data %>%
  select(-c(1,2)) %>%
  sample_n(5000) %>%
  mutate_all(function(x) log(1+x))

###

pca = prcomp(data, center = TRUE, scale. = TRUE, rank. = 30)
Z = pca$x
Z_dist = as.matrix(dist(Z))

###

n = dim(Z)[1]
p = dim(Z)[2]
r = 5

Y = Z[,c(1,2)]
Y_dist = as.matrix(dist(Y))

###

neighbors = seq(from=10, to=300, by=10)
b = 20
k = 20

trusts = trusts_noise = vector(length = length(neighbors)*b)

print("Starting loop!")

count = 1
for (i in 1:length(neighbors)) {
  for (j in 1:b) {
    trans_umap = umap(Z, method = "umap-learn", n_neighbors = neighbors[i], n_components = 2, preserve.seed=FALSE)
    X_umap = trans_umap$layout
    X_dist = as.matrix(dist(X_umap))
    
    trusts[count] = trustworthiness_full2(Y_dist, X_dist, k)
    trusts_noise[count] = trustworthiness_full2(Z_dist, X_dist, k)
    
    print(paste0(count, "/", length(neighbors)*b, " loops complete!"))
    
    count = count + 1
  }
}

###

df = data.frame(n_neighbors = rep(neighbors, each = b),
                trust = trusts,
                trust_noise = trusts_noise)

# save(Y, df, p1, p2, file = trust vs shep (umap CyTOF).Rda")
save(Y, df, file = "~/DR noise/Output/trust vs shep (umap CyTOF 2).Rda")
