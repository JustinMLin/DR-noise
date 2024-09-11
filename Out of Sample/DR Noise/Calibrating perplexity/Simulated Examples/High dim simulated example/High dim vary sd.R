library(Rtsne)
library(dplyr)
library(MASS)

source("Algorithms/DR metrics.R")
# source("../../../Algorithms/DR metrics.R")

n = 350
r = 7
p = 60

Y = rbind(mvrnorm(n=50, mu=c(10,0,0,0,0,0,0), Sigma=diag(abs(rnorm(n=7, mean=0, sd=7)))),
          mvrnorm(n=50, mu=c(0,10,0,0,0,0,0), Sigma=diag(abs(rnorm(n=7, mean=0, sd=7)))),
          mvrnorm(n=50, mu=c(0,0,10,0,0,0,0), Sigma=diag(abs(rnorm(n=7, mean=0, sd=7)))),
          mvrnorm(n=50, mu=c(0,0,0,10,0,0,0), Sigma=diag(abs(rnorm(n=7, mean=0, sd=7)))),
          mvrnorm(n=50, mu=c(0,0,0,0,10,0,0), Sigma=diag(abs(rnorm(n=7, mean=0, sd=7)))),
          mvrnorm(n=50, mu=c(0,0,0,0,0,10,0), Sigma=diag(abs(rnorm(n=7, mean=0, sd=7)))),
          mvrnorm(n=50, mu=c(0,0,0,0,0,0,10), Sigma=diag(abs(rnorm(n=7, mean=0, sd=7)))))
noise = matrix(rnorm(n*p, sd = 1), nrow = n, ncol = p)

perplexity = seq(from = 40, to = 100, by = 5)
perplexity_noise = seq(from = 20, to = 80, by = 5)
b = 40
k = 20

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

df_sd1 = data.frame(perplexity = rep(perplexity, each = b),
                trust = trusts,
                shep = sheps)
df_sd1_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
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

df_sd2 = data.frame(perplexity = rep(perplexity, each = b),
                    trust = trusts,
                    shep = sheps)
df_sd2_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

####### sd = 3 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 3*noise

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

df_sd3 = data.frame(perplexity = rep(perplexity, each = b),
                    trust = trusts,
                    shep = sheps)
df_sd3_noise = data.frame(perplexity = rep(perplexity_noise, each = b),
                          trust = trusts_noise,
                          shep = sheps_noise)

####### sd = 3.5 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 3.5*noise

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

####### sd = 4 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 4*noise

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

####### sd = 4.5 #######
Z = cbind(Y, matrix(0, nrow = n, ncol = p-r)) + 4.5*noise

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
     file = "~/DR noise/Output/High dim vary sd.Rda")