library(dplyr)
library(ggplot2)
library(umap)
library(reticulate)

use_python("/N/soft/rhel8/python/gnu/3.10.5/bin/python")
py_config()
py_available()

source("Algorithms/DR metrics.R")

load("Data/10x data.Rda")

###
Z = unname(pca$x)

n = dim(Z)[1]
p = dim(Z)[2]
r = 4

Y = unname(pca$x)[,1:r]

###

neighbors = seq(from=50, to=200, by=10)
b = 20
k = 50

trusts = sheps = vector(length = length(neighbors)*b)

print("Starting loop!")

count = 1
best_X = NA
best_n_neighbors = NA
best_trust = Inf

for (i in 1:length(neighbors)) {
  for (j in 1:b) {
    trans_umap = umap(Z, method = "umap-learn", n_neighbors = neighbors[i], n_components = 2)
    X_umap = trans_umap$layout
    
    sample_indices = sample(1:n, 1000)
    
    trust = trustworthiness_full_approx(Y, X_umap, k, sample_indices)
    shep = dist_cor_full_approx(Y, X_umap, sample_indices)
    
    if (trust < best_trust) {
      best_X = X_umap
      best_trust = trust
      best_n_neighbors = neighbors[i]
    }
    
    trusts[count] = trust
    sheps[count] = shep
    
    print(paste0(count, "/", length(neighbors)*b, " loops complete!"))
    
    count = count + 1
  }
}

df = data.frame(n_neighbors = rep(neighbors, each = b),
                trust = trusts,
                shep = sheps)

save(Y, df, best_X, best_n_neighbors, best_trust, file = "~/DR noise/Output/10x example.Rda")