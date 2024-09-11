library(Rtsne)
library(dplyr)
library(rgl)

source("Algorithms/DR metrics.R")
# source("../../../Algorithms/DR metrics.R")

n = 500
r = 3
p = 10

theta = seq(from = 0, to = 2*pi, length.out = n/2 + 1)[-1]

x1 = 5*cos(theta)
y1 = 5*sin(theta)
z1 = rep(0, n/2)

x2 = 5*cos(theta) + 5
y2 = rep(0, n/2)
z2 = 5*sin(theta)

Y = matrix(c(x1, x2, y1, y2, z1, z2), nrow = n, ncol = 3)

# plot3d(x=Y[,1], y=Y[,2], z=Y[,3], col=rep(c("#00BFC4","#F8766D"), each=250), xlab="", ylab="", zlab="", axes=FALSE)
# rgl.snapshot("links.png")

noise = matrix(rnorm(n*p, sd = 1), nrow = n, ncol = p)

perplexity = seq(from = 70, to = 150, by = 5)
perplexity_noise = seq(from = 10, to = 60, by = 5)
b = 40
k = 10

####### sd = 0.5 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 0.5*noise

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity_noise)*b)

for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i])$Y
    
    trusts[(i-1)*b + j] = trustworthiness_full(Y, X_tsne, k)
    sheps[(i-1)*b + j] = dist_cor_full(Y, X_tsne)
  }
}

for (i in 1:length(perplexity_noise)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity_noise[i])$Y
    
    trusts_noise[(i-1)*b + j] = trustworthiness_full(Z, X_tsne, k)
    sheps_noise[(i-1)*b + j] = dist_cor_full(Z, X_tsne)
  }
}

df_sd1 = data.frame(perplexity = rep(perplexity, each = b),
                trust = trusts,
                shep = sheps)
df_sd1_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

####### sd = 1 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + noise

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity_noise)*b)

best_trust = -Inf
best_trust_noise = -Inf

for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i])$Y
    
    trusts[(i-1)*b + j] = trustworthiness_full(Y, X_tsne, k)
    sheps[(i-1)*b + j] = dist_cor_full(Y, X_tsne)
    
    if (trusts[(i-1)*b + j] > best_trust) {
      best_X = X_tsne
      best_trust = trusts[(i-1)*b + j]
    }
  }
}

for (i in 1:length(perplexity_noise)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity_noise[i])$Y
    
    trusts_noise[(i-1)*b + j] = trustworthiness_full(Z, X_tsne, k)
    sheps_noise[(i-1)*b + j] = dist_cor_full(Z, X_tsne)
    
    if (trusts_noise[(i-1)*b + j] > best_trust_noise) {
      best_X_noise = X_tsne
      best_trust_noise = trusts_noise[(i-1)*b + j]
    }
  }
}

df_sd2 = data.frame(perplexity = rep(perplexity, each = b),
                    trust = trusts,
                    shep = sheps)
df_sd2_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

####### sd = 1.5 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 1.5*noise

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity_noise)*b)

for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i])$Y
    
    trusts[(i-1)*b + j] = trustworthiness_full(Y, X_tsne, k)
    sheps[(i-1)*b + j] = dist_cor_full(Y, X_tsne)
  }
}

for (i in 1:length(perplexity_noise)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity_noise[i])$Y
    
    trusts_noise[(i-1)*b + j] = trustworthiness_full(Z, X_tsne, k)
    sheps_noise[(i-1)*b + j] = dist_cor_full(Z, X_tsne)
  }
}

df_sd3 = data.frame(perplexity = rep(perplexity, each = b),
                    trust = trusts,
                    shep = sheps)
df_sd3_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

####### sd = 2 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 2*noise

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity_noise)*b)

for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i])$Y
    
    trusts[(i-1)*b + j] = trustworthiness_full(Y, X_tsne, k)
    sheps[(i-1)*b + j] = dist_cor_full(Y, X_tsne)
  }
}

for (i in 1:length(perplexity_noise)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity_noise[i])$Y
    
    trusts_noise[(i-1)*b + j] = trustworthiness_full(Z, X_tsne, k)
    sheps_noise[(i-1)*b + j] = dist_cor_full(Z, X_tsne)
  }
}

df_sd4 = data.frame(perplexity = rep(perplexity, each = b),
                    trust = trusts,
                    shep = sheps)
df_sd4_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

####### sd = 2.5 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 2.5*noise

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity_noise)*b)

for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i])$Y
    
    trusts[(i-1)*b + j] = trustworthiness_full(Y, X_tsne, k)
    sheps[(i-1)*b + j] = dist_cor_full(Y, X_tsne)
  }
}

for (i in 1:length(perplexity_noise)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity_noise[i])$Y
    
    trusts_noise[(i-1)*b + j] = trustworthiness_full(Z, X_tsne, k)
    sheps_noise[(i-1)*b + j] = dist_cor_full(Z, X_tsne)
  }
}

df_sd5 = data.frame(perplexity = rep(perplexity, each = b),
                    trust = trusts,
                    shep = sheps)
df_sd5_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

####### sd = 3 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 3*noise

trusts = sheps = vector(length = length(perplexity)*b)
trusts_noise = sheps_noise = vector(length = length(perplexity_noise)*b)

for (i in 1:length(perplexity)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity[i])$Y
    
    trusts[(i-1)*b + j] = trustworthiness_full(Y, X_tsne, k)
    sheps[(i-1)*b + j] = dist_cor_full(Y, X_tsne)
  }
}

for (i in 1:length(perplexity_noise)) {
  for (j in 1:b) {
    X_tsne = Rtsne(Z, perplexity = perplexity_noise[i])$Y
    
    trusts_noise[(i-1)*b + j] = trustworthiness_full(Z, X_tsne, k)
    sheps_noise[(i-1)*b + j] = dist_cor_full(Z, X_tsne)
  }
}

df_sd6 = data.frame(perplexity = rep(perplexity, each = b),
                    trust = trusts,
                    shep = sheps)
df_sd6_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

save(df_sd1, df_sd1_noise,
     df_sd2, df_sd2_noise,
     df_sd3, df_sd3_noise,
     df_sd4, df_sd4_noise,
     df_sd5, df_sd5_noise,
     df_sd6, df_sd6_noise,
     best_trust, best_X,
     best_trust_noise, best_X_noise,
     file = "~/DR noise/Output/Links vary sd.Rda")