library(dplyr)
library(Rtsne)
library(EFAtools)

# source("../../../Algorithms/DR metrics.R")
source("Algorithms/DR metrics.R")

# data = read.delim("../../../../../scRNA-seq data/molecules.txt")
data = read.delim("Data/molecules.txt")

data = as.data.frame(t(data))
data = data %>%
  dplyr::select(where(function(x) mean(x == 0) < 0.2)) %>%
  mutate_all(function(x) log(1+x))

pca = prcomp(data, center = TRUE, scale. = TRUE, rank. = 30)
Z = pca$x

###

n = dim(Z)[1]
p = dim(Z)[2]
r = 10

Y = Z[,(1:r)]

###

b = 20
k = 50

best_trust = -Inf
best_trust_noise = -Inf

print("Starting loop!")

for (j in 1:b) {
  X_tsne = Rtsne(Z, perplexity = 60)$Y
  X_tsne_noise = Rtsne(Z, perplexity = 50)$Y
  
  trust = trustworthiness_full_approx(Y, X_tsne, k)
  trust_noise = trustworthiness_full_approx(Z, X_tsne_noise, k)
  
  if (trust > best_trust) {
    best_X = X_tsne
    best_trust = trust
  }
  
  if (trust_noise > best_trust_noise) {
    best_X_noise = X_tsne_noise
    best_trust_noise = trust_noise
  }
}

###

# save(best_X, best_X_noise, best_trust, best_trust_noise, file = "best rep (10 dim scRNA-seq).Rda")
save(best_X, best_X_noise, best_trust, best_trust_noise, file = "~/DR noise/Output/best rep (10 dim scRNA-seq).Rda")