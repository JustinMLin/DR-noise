dist_sq = function(x, y) {
  sum((x-y)^2)
}

spherical_dist_sq = function(r,p1,p2) {
  theta1 = p1[1]
  phi1 = p1[2]
  theta2 = p2[1]
  phi2 = p2[2]
  sigma = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(theta1-theta2))
  
  (r*sigma)^2
}

compute_dist_sq_mat = function(data) {
  dist_mat = matrix(nrow = length(data[,1]), ncol = length(data[,1]))
  
  for (i in 1:length(data[,1])) {
    for (j in 1:i) {
      dist_mat[i,j] = dist_mat[j,i] = dist_sq(data[i,], data[j,])
    }
  }
  
  return(dist_mat)
}

sq_norm = function(x) {
  sum(x^2)
}

rel_error = function(pt, data, a2) {
  n = length(data[,1])
  ret = vector(length = n)
  
  for (i in 1:n) {
    ret[i] = abs(dist_sq(pt, data[i,]) - a2[i])/a2[i]
  }
  
  ret
}

abs_error = function(pt, data, a2) {
  n = length(data[,1])
  ret = vector(length = n)
  
  for (i in 1:n) {
    ret[i] = abs(dist_sq(pt, data[i,]) - a2[i])
  }
  
  ret
}

################################################################################
# Min Cost

# proj_data = coordinates of embedding of original data (n x d matrix)
# a2 = squared dissimilarities of new object from original n object (n x 1 vector)
# range: optimization starting point uniform over [-range,range]^d


min_cost = function(proj_data, a2, range = 5, iter = 10000000) {
  dim = length(proj_data[1,])
  
  min_func = function(y) {
    y_mat = matrix(data = y, nrow = length(proj_data[,1]), ncol = length(y), byrow = TRUE)
    
    sq_diff_coord = (y_mat - proj_data)^2
    
    y_dists = sqrt(rowSums(sq_diff_coord))
    
    sum((y_dists - sqrt(a2))^2)
  }
  
  best = optim(par = rep(0,dim), fn = min_func)
  rep({current = optim(par = runif(dim,-range,range), fn = min_func)
  if (current$val < best$val) {
    best = current
  }},
  iter)
  
  best
}

################################################################################
# Sammon Min Cost

# proj_data = coordinates of embedding of original data (n x d matrix)
# a2 = squared dissimilarities of new object from original n object (n x 1 vector)
# range: optimization starting point uniform over [-range,range]^d


min_cost_Sammon = function(proj_data, a2, range = 5) {
  dim = length(proj_data[1,])
  
  min_func_Sammon = function(y) {
    y_mat = matrix(data = y, nrow = length(proj_data[,1]), ncol = length(y), byrow = TRUE)
    
    sq_diff_coord = (y_mat - proj_data)^2
    
    y_dists = sqrt(rowSums(sq_diff_coord))
    
    sum((y_dists - sqrt(a2))^2/sqrt(a2))
  }
  
  best = optim(par = rep(0,dim), fn = min_func_Sammon)
  rep({current = optim(par = runif(dim,-range,range), fn = min_func_Sammon)
  if (current$val < best$val) {
    best = current
  }},
  10000000)
  
  best
}

################################################################################
# Trosset

# proj_data = coordinates of embedding of original data (n x d matrix)
# D2 = squared dissimilarities for original n objects (n x n matrix)
# a2 = squared dissimilarities of new object from original n object (n x 1 vector)
# range: optimization starting point uniform over [-range,range]^d


Trosset = function(proj_data, D2, a2, range = 5) {
  n = length(a2)
  dim = length(proj_data[1,])
  
  A2 = cbind(rbind(D2, a2), c(a2,0))
  
  e = rep(1,n+1)
  w = c(rep(1,n),0)
  
  B = -0.5 * (diag(n+1) - e %*% t(w) / n) %*% A2 %*% (diag(n+1) - w %*% t(e) / n)
  
  min_func = function(y) {
    Y = rbind(proj_data, y)
    norm(B - Y %*% t(Y), type = "F")^2
  }
  
  best = optim(par = rep(0,dim), fn = min_func)
  rep({current = optim(par = runif(dim,-range,range), fn = min_func)
  if (current$val < best$val) {
    best = current
  }},
  10000000)
  
  best
}

################################################################################
# Trosset convex

# proj_data = coordinates of embedding of original data (n x d matrix)
# D2 = squared dissimilarities for original n objects (n x n matrix)
# a2 = squared dissimilarities of new object from original n object (n x 1 vector)


Trosset_convex = function(proj_data, D2, a2) {
  n = length(a2)
  dim = length(proj_data[1,])
  
  A2 = cbind(rbind(D2, a2), c(a2,0))
  
  e = rep(1,n+1)
  w = c(rep(1,n),0)
  
  B = -0.5 * (diag(n+1) - e %*% t(w) / n) %*% A2 %*% (diag(n+1) - w %*% t(e) / n)
  
  b = B[1:n,n+1]
  
  min_func = function(y) {
    sum((b - proj_data %*% y)^2)
  }
  
  optim(par = rep(0,dim), fn = min_func)
}


################################################################################
# Bengio

K_mds = function(a,b,data) {
  n = length(data[,1])
  
  exp_with_a = dist_sq(a,t(data))/n
  exp_with_b = dist_sq(b,t(data))/n
  
  exp = 0
  for (i in 1:n) {
    exp = exp + dist_sq(data[i,],t(data))
  }
  exp = exp/n^2
  
  -0.5 * (dist_sq(a,b) - exp_with_a - exp_with_b + exp)
}

Bengio = function(data, x, proj_dim, K) {
  n = length(data[,1])
  d = proj_dim
  
  e = rep(1,n)
  C = diag(n) - e %*% t(e)/n
  D = compute_dist_sq_mat(data)
  
  M = -0.5 * C %*% D %*% C
  
  eigen_system = eigen(M)
  lambda =  eigen(M)$values[1:d]
  if (any(lambda <= 0)) {
    stop("Not enough positive eigenvalues")
  }
  
  v = eigen(M)$vectors[,1:d]
  
  Kx = vector(length =  n)
  for (i in 1:n) {
    Kx[i] = K(x, data[i,], data)
  }
  
  e = vector(length = d)
  for (i in 1:d) {
    e[i] = 1/sqrt(lambda[i]) * t(v[,i]) %*% Kx
  }
  
  e
}

Bengio_data = function(data, proj_dim, K) {
  n = length(data[,1])
  d = proj_dim
  
  ret_mat = matrix(nrow = n, ncol = d)
  
  M = matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:i) {
      M[i,j] = M[j,i] = K(data[i,], data[j,], data)
    }
  }
  
  eigen_system = eigen(M)
  lambda =  eigen(M)$values[1:d]
  if (any(lambda <= 0)) {
    stop("Not enough positive eigenvalues")
  }
  
  v = eigen(M)$vectors[,1:d]
  
  for (i in 1:n) {
    Kx = vector(length = n)
    for (j in 1:n) {
      Kx[j] = K(data[i,], data[j,], data)
    }
    
    e = vector(length = d)
    for (j in 1:d) {
      e[j] = 1/sqrt(lambda[j]) * t(v[,j]) %*% Kx
    }
    
    ret_mat[i,] = e
  }
  
  ret_mat
}
