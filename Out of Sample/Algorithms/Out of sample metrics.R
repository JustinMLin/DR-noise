library(dbscan)

################################################################################
## Helper Functions

get_nns = function(pt, data, k) {
  kNN(data, k, matrix(pt, nrow = 1))$id
}

nn_rank = function(pt, j, data) {
  n = length(data[,1])
  rank = match(j, kNN(data, n-1, matrix(pt, nrow = 1))$id)
  
  if (is.na(rank)) n else rank
}

d = function(x, y) {
  norm(x - y, type = "2")
}

################################################################################
## Out of Sample Metrics

trustworthiness = function(Z, w, X, y, k) {
  n = length(Z[,1])
  
  high_dim_neighbors = get_nns(w, Z, k)
  low_dim_neighbors = get_nns(y, X, k)
  
  U = setdiff(low_dim_neighbors, high_dim_neighbors)
  
  if (length(U) == 0) {
    1
  }
  else {
    1 - 2/(k*(2*n - 3*k + 1)) * sum(sapply(U, function(j) nn_rank(w, j, Z)) - k)
  }
}

continuity = function(Z, w, X, y, k) {
  n = length(Z[,1])
  
  high_dim_neighbors = get_nns(w, Z, k)
  low_dim_neighbors = get_nns(y, X, k)
  
  V = setdiff(high_dim_neighbors, low_dim_neighbors)
  
  if (length(V) == 0) {
    1
  }
  else {
    1 - 2/(k*(2*n - 3*k + 1)) * sum(sapply(V, function(j) nn_rank(y, j, X)) - k)
  }
}

precision = function(Z, w, X, y, k) {
  high_dim_neighbors = get_nns(w, Z, k)
  
  dw = vector(length = k)
  dy = vector(length = k)
  
  for (i in 1:k) {
    dw[i] = d(w, Z[high_dim_neighbors[i],])
    dy[i] = d(y, X[high_dim_neighbors[i],])
  }
  
  1 - 0.5*(norm(dw/norm(dw, type = "2") - dy/norm(dy, type = "2"), type = "2"))
}

stress = function(Z, w, X, y) {
  n = length(Z[,1])
    
  w_mat = matrix(w, nrow = n, ncol = length(w), byrow = TRUE)
  w_dists = sqrt(rowSums((Z - w_mat)^2))
  
  y_mat = matrix(y, nrow = n, ncol = length(y), byrow = TRUE)
  y_dists = sqrt(rowSums((X - y_mat)^2))
  
  sum((w_dists - y_dists)^2) / sum((w_dists)^2)
}

local_stress = function(Z, w, X, y, k) {
  neighbors = get_nns(w, Z, k)
  
  w_mat = matrix(w, nrow = k, ncol = length(w), byrow = TRUE)
  w_dists = sqrt(rowSums((Z[neighbors,] - w_mat)^2))
  
  y_mat = matrix(y, nrow = k, ncol = length(y), byrow = TRUE)
  y_dists = sqrt(rowSums((X[neighbors,] - y_mat)^2))
  
  sum((w_dists - y_dists)^2) / sum((w_dists)^2)
}

dist_cor = function(Z, w, X, y) {
  n = length(Z[,1])
  
  w_mat = matrix(w, nrow = n, ncol = length(w), byrow = TRUE)
  w_dists = sqrt(rowSums((Z - w_mat)^2))
  
  y_mat = matrix(y, nrow = n, ncol = length(y), byrow = TRUE)
  y_dists = sqrt(rowSums((X - y_mat)^2))
  
  cor(w_dists, y_dists, method = "spearman")
}