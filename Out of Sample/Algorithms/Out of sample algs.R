library(RSpectra)
library(ucminf)
library(minqa)

###############################################################################
## Min Cost

min_cost = function(X, delta, range = 10, iter = 1e+7) {
  n = length(X[,1])
  q = length(X[1,])
  
  min_func = function(y) {
    y_mat = matrix(data = y, nrow = n, ncol = q, byrow = TRUE)
    
    sq_diff_coord = (y_mat - X)^2
    
    y_dists_sq = rowSums(sq_diff_coord)
    
    sum(abs(delta-y_dists_sq))
  }
  
  minqa::bobyqa(par = rep(0,q), fn = min_func)$par
}

min_cost_sq = function(X, delta, range = 10, iter = 1e+7) {
  n = length(X[,1])
  q = length(X[1,])
  
  min_func = function(y) {
    y_mat = matrix(data = y, nrow = n, ncol = q, byrow = TRUE)
    
    sq_diff_coord = (y_mat - X)^2
    
    y_dists = sqrt(rowSums(sq_diff_coord))
    
    sum((sqrt(delta)-y_dists)^2)
  }
  
  minqa::bobyqa(par = rep(0,q), fn = min_func)$par
}

###############################################################################
## Trosset

tau = function(M) {
  n = length(M[1,]) - 1
  e = rep(1, n+1)
  w = c(rep(1,n), 0)
  
  -1/2 * (diag(n+1) - (e %*% t(w))/n) %*% M %*% (diag(n+1) - (w %*% t(e))/n)
}

Trosset = function(X, A2, range = 10, iter = 1e+7) {
  n = length(X[,1])
  q = length(X[1,])
  
  B = tau(A2)
  b = B[1:n,n+1]
  beta = B[n+1,n+1]
  
  min_func = function(y) {
    2 * norm(b - X %*% y, type="F")^2 + (beta - norm(y, type="2")^2)^2
  }
  
  grad = function(y) {
    -4 * t(X) %*% b + 4 * t(X) %*% X %*% y - 4 * beta * y + 4 * norm(y, type="2")^2 * y 
  }
  
  pca = prcomp(Z, rank. = 2, scale. = FALSE, center = FALSE)
  start = predict(pca, matrix(w, nrow = 1))
  
  optim(par = start, fn = min_func, gr = grad, method = "BFGS")$par
}

Trosset_approx = function(X, A2) {
  n = length(X[,1])
  B = tau(A2)
  b = B[1:n,n+1]
  
  c(solve(t(X) %*% X) %*% t(X) %*% b)
}

###############################################################################
## Bengio

K_mds_lp = function(p) {
  K = function(a, b, data) {
    n = length(data[,1])
    
    exp_with_a = sum(colSums(abs(a-t(data))^p)^(2/p))/n
    exp_with_b = sum(colSums(abs(b-t(data))^p)^(2/p))/n
    
    exp = 0
    for (i in 1:n) {
      exp = exp + sum(colSums(abs(data[i,]-t(data))^p)^(2/p))
    }
    exp = exp/n^2
    
    -0.5 * (sum(abs(a-b)^p)^(2/p) - exp_with_a - exp_with_b + exp)
  }
  
  return(K)
}

colMax = function(data) apply(data, 2, max, na.rm = TRUE)

K_mds_inf = function(a, b, data) {
  n = length(data[,1])
  
  exp_with_a = sum(colMax(abs(a-t(data)))^2)/n
  exp_with_b = sum(colMax(abs(b-t(data)))^2)/n
  
  exp = 0
  for (i in 1:n) {
    exp = exp + sum(colMax(abs(data[i,]-t(data)))^2)
  }
  exp = exp/n^2
  
  -0.5 * (max(abs(a-b))^2 - exp_with_a - exp_with_b + exp)
}

Bengio = function(Z, D2, w, q, K) {
  n = length(Z[,1])
  
  e = rep(1,n)
  C = diag(n) - e %*% t(e)/n
  
  M = -0.5 * C %*% D2 %*% C
  
  eigen_system = eigs(M, q)
  lambda =  eigen_system$values
  if (any(lambda <= 0)) {
    stop("Not enough positive eigenvalues")
  }
  
  v = eigen_system$vectors
  
  Kw = vector(length =  n)
  for (i in 1:n) {
    Kw[i] = K(w, Z[i,], Z)
  }
  
  diag(1/sqrt(lambda)) %*% t(v) %*% Kw
}


Bengio_multiple = function(Z, D2, W, q, K) {
  n = length(Z[,1])
  m = length(W[,1])
  
  e = rep(1,n)
  C = diag(n) - e %*% t(e)/n

  M = -0.5 * C %*% D2 %*% C
  
  eigen_system = eigs(M, q)
  lambda =  eigen_system$values
  if (any(lambda <= 0)) {
    stop("Not enough positive eigenvalues")
  }
  
  v = eigen_system$vectors
  
  KW = matrix(nrow = m, ncol = n)
  for (i in 1:m) {
    for (j in 1:n) {
      KW[i,j] = K(W[i,], Z[j,], Z)
    }
  }
  
  KW %*% v %*% diag(1/sqrt(lambda))
}