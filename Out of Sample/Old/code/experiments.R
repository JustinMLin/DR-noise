source("out of sample algs.R")
library(ggplot2)
library(gridExtra)
library(rgl)

temp = runif(50,-10,10)
data = matrix(c(temp,temp^2,rep(0,50)), ncol = 3)
data[,1] = data[,1] + rnorm(50,0,2)
data[,2] = data[,2] + rnorm(50,0,2)

n = length(data[,1])

Bengio_data = matrix(nrow = n, ncol = 2)
for (i in 1:n) {
  Bengio_data[i,] = Bengio(data, data[i,], proj_dim = 2, K = K_mds)
}

x = c(5,25,runif(1,-10,10) + 1)

a2 = vector(length = n)
for (i in 1:n) {
  a2[i] = dist_sq(x, data[i,])
}
D2 = compute_dist_sq_mat(Bengio_data)

Trosset_pt = Trosset(Bengio_data, D2, a2, range = 5)$par
min_cost_pt = min_cost(Bengio_data, a2, range = 5)$par
Bengio_pt = Bengio(data, x, 2, K_mds)

d = data.frame(x = c(Bengio_data[,1], Trosset_pt[1], min_cost_pt[1], Bengio_pt[1]), 
           y = c(Bengio_data[,2], Trosset_pt[2], min_cost_pt[2], Bengio_pt[2]),
           type = c(rep("original", n), "Trosset", "Min Cost", "Bengio"))
p = ggplot(d, aes(x = x, y = y, col = type)) + geom_point()
p



###############################################################################

n = 60
proj_data = matrix(rnorm(2*n, -10, 10), ncol = 2)
D2 = compute_dist_sq_mat(proj_data)

rand_point = runif(2, -20, 20)
rand_point_dist = vector(length = n)
for (i in 1:n) {
  rand_point_dist[i] = dist_sq(rand_point, proj_data[i,])
}

a2 = rand_point_dist + runif(n, -10, 10)

Trosset_pt = Trosset(proj_data, D2, a2, range = 20)
min_cost_pt = min_cost(proj_data, a2, range = 20)

d = data.frame(x = c(proj_data[,1], Trosset_pt$par[1], min_cost_pt$par[1]), 
           y = c(proj_data[,2], Trosset_pt$par[2], min_cost_pt$par[2]),
           type = c(rep("Original", n), "Trosset", "Min Cost"))
ggplot(data = d, aes(x = x, y = y, col = type)) + geom_point()


Trosset_dists = vector(length = n)
min_cost_dists = vector(length= n)
for (i in 1:n) {
  Trosset_dists[i] = dist_sq(Trosset_pt$par, proj_data[i,])
  min_cost_dists[i] = dist_sq(min_cost_pt$par, proj_data[i,])
}

Trosset_data = data.frame(x = sqrt(a2), 
                          y = sqrt(Trosset_dists), 
                          subtracted = sqrt(Trosset_dists) - sqrt(a2))
Trosset_data$group = "Trosset"
min_cost_data = data.frame(x = sqrt(a2), 
                           y = sqrt(min_cost_dists), 
                           subtracted = sqrt(min_cost_dists) - sqrt(a2))
min_cost_data$group = "Min Cost"

combine = rbind(Trosset_data,min_cost_data)

ggplot(data = combine, aes(x = x, y = subtracted, group = group, col = group)) +
  geom_point() + 
  geom_smooth(se = FALSE)

###############################################################################

experiment = function(n) {
  proj_data = matrix(runif(2*n, -10, 10), ncol = 2)
  D2 = compute_dist_sq_mat(proj_data)
  
  rand_point = c(runif(2, -20, 20), runif(1, 0, 20))
  a2 = vector(length = n)
  for (i in 1:n) {
    a2[i] = dist_sq(c(proj_data[i,],0), rand_point)
  }
  
  new_D2 = rbind(cbind(D2,a2),c(a2,0))
  e = rep(1,n+1)
  P = diag(n+1) - e %*% t(e) / (n+1)
  B = -0.5 * P %*% new_D2 %*% P
  d = svd(B)$d
  emb = round((d[1] + d[2]) / sum(d), digits = 5)
  
  Trosset_pt = Trosset(proj_data, D2, a2, range = 20)
  Trosset_error = rel_error(Trosset_pt$par, proj_data, a2)
  
  min_cost_pt = min_cost(proj_data, a2, range = 20)
  min_cost_error = rel_error(min_cost_pt$par, proj_data, a2)
  
  data.frame(Ddata_sq = a2, Trosset_error, min_cost_error, emb = rep(emb,n))
}

get_data = function(n, runs) {
  data = experiment(n)
  
  for (i in 2:runs) {
    data = rbind(data, experiment(n))
  }
  data
}

data = get_data(60,10)

p = ggplot(data, aes(x = Ddata_sq, y = Trosset_error, col = "Trosset")) + 
  geom_point() +
  geom_smooth(se = FALSE) + 
  geom_point(data = data, aes(x = Ddata_sq, y = min_cost_error, col = "Min Cost")) +
  geom_smooth(data = data, aes(x = Ddata_sq, y = min_cost_error, col = "Min Cost"), se = FALSE)
p + facet_grid(rows = vars(emb), scale = "free_y")


###############################################################################

sim = function(n = 20, b = 100, base_pt = c(runif(2,-10,10),0)) {
  res_sq = matrix(nrow = b, ncol = 2)
  
  data = matrix(runif(2*n,-10,10), nrow = n)
  data3D = cbind(data, rep(0,n))
  Bengio_base_pt = Bengio(data3D, base_pt, 2, K_mds)
  
  for (j in 1:b) {
    error = rnorm(3, 0, 10)
    new_pt = base_pt + error
    
    Bengio_pt = Bengio(data3D, new_pt, 2, K_mds)
    
    a2 = vector(length = n)
    for (i in 1:n) {
      a2[i] = dist_sq(new_pt, data3D[i,])
    }
    
    min_cost_pt = min_cost(data, a2, range = 15, iter = +1e06)$par
    
    res_sq[j,1] = dist_sq(min_cost_pt, base_pt[-3])
    res_sq[j,2] = dist_sq(Bengio_pt, Bengio_base_pt[-3])
  }
  res_sq
}

for(i in seq(from = -8, to = 8, by = 8)) {
  for(j in seq(from = -8, to = 8, by = 8)) {
    name = paste("p",i+8,j+8,sep = "")
    
    res = sim(base_pt = c(i,j,0))
    data = data.frame(x = c(res[,1],res[,2]), 
                      type = c(rep("Min Cost", 100), rep("Bengio", 100)))
    assign(name, ggplot(data, aes(x = x, fill = type)) + 
             geom_histogram(position = "identity", alpha = 0.3, bins = 50))
  }
}
grid.arrange(p016,p816,p1616,p08,p88,p168,p00,p80,p160, ncol = 3)

###############################################################################

sim = function(n = 20, b = 100, base_pt, sd) {
  res_sq = matrix(nrow = b, ncol = 2)
  
  r = 5*sqrt(runif(n))
  theta = runif(n,0,2*pi)
  
  data = matrix(c(r*cos(theta), r*sin(theta)), ncol = 2)
  data3D = cbind(data, rep(0,n))
  Bengio_base_pt = Bengio(data3D, base_pt, 2, K_mds)
  
  for (j in 1:b) {
    error = rnorm(3, 0, sd)
    new_pt = base_pt + error
    
    Bengio_pt = Bengio(data3D, new_pt, 2, K_mds)
    
    a2 = vector(length = n)
    for (i in 1:n) {
      a2[i] = dist_sq(new_pt, data3D[i,])
    }
    
    min_cost_pt = min_cost(data, a2, range = 15, iter = +1e06)$par
    
    res_sq[j,1] = dist_sq(min_cost_pt, base_pt[-3])
    res_sq[j,2] = dist_sq(Bengio_pt, Bengio_base_pt[-3])
  }
  res_sq
}

build_box = function(sds, base_pt, data_size = 20) {
  n = length(sds)
  mat = matrix(nrow = n*100, ncol = 3)
    
  for(i in 1:n) {
    mat[((i-1)*100+1):(i*100),1] = rep(sds[i],100)
    mat[((i-1)*100+1):(i*100),2:3] = sim(n = data_size, base_pt = base_pt, sd = sds[i])
  }
  table = data.frame(sd = c(mat[,1],mat[,1]),
                     r2 = c(mat[,2],mat[,3]),
                     type = c(rep("Min Cost", n*100), rep("Bengio", n*100)))
  
  ggplot(table, aes(x = factor(sd), y = r2, color = type)) + geom_boxplot()
}

p1 = build_box_diff(c(0.5,1,5,10), c(5,0,0)) + ylim(-20,20) 
p2 = build_box_diff(c(0.5,1,2,5), c(10,0,0)) + ylim(-20,20)
p3 = build_box_diff(c(0.5,1,2,5), c(15,0,0)) + ylim(-20,20)
p4 = build_box_diff(c(0.5,1,2,5), c(25,0,0)) + ylim(-20,20)
grid.arrange(p1,p2,p3,p4,ncol=1)

###############################################################################

sim = function(n = 20, b = 100, new_pt, sd) {
  res_sq = matrix(nrow = b, ncol = 2)
  
  theta = runif(n,0,2*pi)
  phi = runif(n, 0, pi)
  r = 5*(runif(n,0,1))^(1/3)
  
  x = r*sin(phi)*cos(theta)
  y = r*sin(phi)*sin(theta)
  z = r*cos(phi)
  
  data = matrix(c(x,y,z), ncol = 3)
  
  data2D = Bengio_data(data, 2, K_mds)
  
  a2_base = vector(length = n)
  for (i in 1:n) {
    a2_base[i] = dist_sq(new_pt, data[i,])
  }
  min_cost_base_pt = min_cost(data2D, a2_base, range = 15, iter = 100000)$par
  
  Bengio_base_pt = Bengio(data, new_pt, 2, K_mds)
  
  for (j in 1:b) {
    error = rnorm(3, 0, sd)
    error_pt = new_pt + error
    
    a2 = vector(length = n)
    for (i in 1:n) {
      a2[i] = dist_sq(error_pt, data[i,])
    }
    min_cost_pt = min_cost(data2D, a2, range = 15, iter = 100000)$par
    
    Bengio_pt = Bengio(data, error_pt, 2, K_mds)
    
    res_sq[j,1] = dist_sq(min_cost_pt, min_cost_base_pt[-3])
    res_sq[j,2] = dist_sq(Bengio_pt, Bengio_base_pt[-3])
  }
  res_sq
}

build_box = function(sds, base_pt, data_size = 20) {
  n = length(sds)
  mat = matrix(nrow = n*100, ncol = 3)
  
  for(i in 1:n) {
    mat[((i-1)*100+1):(i*100),1] = rep(sds[i],100)
    mat[((i-1)*100+1):(i*100),2:3] = sim(n = data_size, new_pt = base_pt, sd = sds[i])
  }
  table = data.frame(sd = c(mat[,1],mat[,1]),
                     r2 = c(mat[,2],mat[,3]),
                     type = c(rep("Min Cost", n*100), rep("Bengio", n*100)))
  
  ggplot(table, aes(x = factor(sd), y = r2, color = type)) + geom_boxplot()
}

build_box_diff = function(sds, base_pt, data_size = 20) {
  n = length(sds)
  mat = matrix(nrow = n*100, ncol = 3)
  
  for(i in 1:n) {
    mat[((i-1)*100+1):(i*100),1] = rep(sds[i],100)
    mat[((i-1)*100+1):(i*100),2:3] = sim(n = data_size, new_pt = base_pt, sd = sds[i])
  }
  table = data.frame(sd = mat[,1],
                     r2_diff = mat[,2] - mat[,3])
  
  ggplot(table, aes(x = factor(sd), y = r2_diff)) + geom_boxplot()
}

p2 = build_box(c(1,2,5,10,15), c(5,0,0))
p3 = build_box(c(1,2,5,10,15), c(10,0,0))
p4 = build_box(c(1,2,5,10,15), c(15,0,0))
grid.arrange(p2,p3,p4,ncol=1)

###############################################################################

sim = function(n = 20, b = 100, new_pt, sd) {
  res_sq = matrix(nrow = b, ncol = 2)
  
  theta = runif(n,0,2*pi)
  phi = runif(n, 0, pi)
  r = 5*(runif(n,0,1))^(1/3)
  
  x = r*sin(phi)*cos(theta)
  y = r*sin(phi)*sin(theta)
  z = r*cos(phi)
  
  data = matrix(c(x,y,z), ncol = 3)
  
  data2D = Bengio_data(data, 2, K_mds)
  
  a2_base = vector(length = n)
  for (i in 1:n) {
    a2_base[i] = dist_sq(new_pt, data[i,])
  }
  min_cost_base_pt = min_cost(data2D, a2_base, range = 15, iter = 100000)$par
  
  Bengio_base_pt = Bengio(data, new_pt, 2, K_mds)
  
  for (j in 1:b) {
    error = rnorm(3, 0, sd)
    error_pt = new_pt + error
    
    a2 = vector(length = n)
    for (i in 1:n) {
      a2[i] = dist_sq(error_pt, data[i,])
    }
    min_cost_pt = min_cost(data2D, a2, range = 15, iter = +1e06)$par
    
    Bengio_pt = Bengio(data, error_pt, 2, K_mds)
    
    res_sq[j,1] = dist_sq(min_cost_pt, min_cost_base_pt[-3])
    res_sq[j,2] = dist_sq(Bengio_pt, Bengio_base_pt[-3])
  }
  res_sq
}

build_box = function(sds, base_pt, data_size = 20) {
  n = length(sds)
  mat = matrix(nrow = n*100, ncol = 3)
  
  for(i in 1:n) {
    mat[((i-1)*100+1):(i*100),1] = rep(sds[i],100)
    mat[((i-1)*100+1):(i*100),2:3] = sim(n = data_size, new_pt = base_pt, sd = sds[i])
  }
  table = data.frame(sd = c(mat[,1],mat[,1]),
                     r2 = c(mat[,2],mat[,3]),
                     type = c(rep("Min Cost", n*100), rep("Bengio", n*100)))
  
  ggplot(table, aes(x = factor(sd), y = r2, color = type)) + geom_boxplot()
}

p1 = build_box(c(1,2,5), c(2,0,0))
p2 = build_box(c(1,2,5), c(5,0,0))
p3 = build_box(c(1,2,5), c(10,0,0))
p4 = build_box(c(1,2,5), c(15,0,0))
grid.arrange(p1,p2,p3,p4,ncol=1)

###############################################################################

get_rel_errors = function(data, proj_data, base_pt, error_sd) {
  error_pt = base_pt + rnorm(3, 0, error_sd)
  Bengio_pt = Bengio(data, error_pt, 2, K_mds)
  
  n = length(data[,1])
  a2 = vector(length = n)
  for (i in 1:n) {
    a2[i] = dist_sq(error_pt,data[i,])
  }
  min_cost_pt = min_cost(proj_data, a2, range = 15, iter = +1e06)$par
  
  base_pt_dists = vector(length = n)
  Bengio_dists = vector(length = n)
  min_cost_dists = vector(length = n)
  
  for (i in 1:n) {
    base_pt_dists[i] = dist_sq(base_pt, data[i,])
    Bengio_dists[i] = dist_sq(Bengio_pt, proj_data[i,])
    min_cost_dists[i] = dist_sq(min_cost_pt, proj_data[i,])
  }
  
  cbind(abs(base_pt_dists - min_cost_dists)/base_pt_dists, 
        abs(base_pt_dists - Bengio_dists)/base_pt_dists)
}

sim = function(n = 20, b = 100, new_pt, sd) {
  theta = runif(n,0,2*pi)
  phi = runif(n, 0, pi)
  r = 5*(runif(n,0,1))^(1/3)
  
  x = r*sin(phi)*cos(theta)
  y = r*sin(phi)*sin(theta)
  z = r*cos(phi)
  
  data = matrix(c(x,y,z), ncol = 3)
  proj_data = Bengio_data(data, 2, K_mds)
  
  mean_errors = matrix(nrow = b, ncol = 2)
  for (i in 1:b) {
    errors = get_rel_errors(data, proj_data, new_pt, sd)
    mean_errors[i,] = colMeans(errors)
  }
  
  mean_errors
}

build_box = function(sds, base_pt, data_size = 20) {
  n = length(sds)
  mat = matrix(nrow = n*100, ncol = 3)
  
  for(i in 1:n) {
    mat[((i-1)*100+1):(i*100),1] = rep(sds[i],100)
    mat[((i-1)*100+1):(i*100),2:3] = sim(n = data_size, new_pt = base_pt, sd = sds[i])
  }
  table = data.frame(sd = c(mat[,1],mat[,1]),
                     mean_rel_error = c(mat[,2],mat[,3]),
                     type = c(rep("Min Cost", n*100), rep("Bengio", n*100)))
  
  ggplot(table, aes(x = factor(sd), y = mean_rel_error, color = type)) + geom_boxplot()
}

p1 = build_box(c(1,2,5,10), c(2,0,0))
p2 = build_box(c(1,2,5,10), c(5,0,0))
p3 = build_box(c(1,2,5,10), c(10,0,0))
grid.arrange(p1,p2,p3,ncol=1)

###############################################################################

get_rel_errors = function(spherical_data, proj_data, base_pt, error_sd) {
  n = length(spherical_data[,1])
  
  error_theta = rnorm(1, 0, error_sd)
  error_phi = rnorm(1, 0, error_sd/4)
  error_pt = base_pt + c(error_theta, error_phi)
  error_pt_cart = c(5*sin(error_pt[2]*cos(error_pt[1])), 
                    5*sin(error_pt[2]*sin(error_pt[1])),
                    5*cos(error_pt[2]))
  
  x = 5*sin(spherical_data[,2])*cos(spherical_data[,1])
  y = 5*sin(spherical_data[,2])*sin(spherical_data[,1])
  z = 5*cos(spherical_data[,2])
  
  cartesian_data = cbind(x,y,z)
  
  Bengio_error_pt = Bengio(cartesian_data, error_pt_cart, 2, K_mds)
  
  a2_error_pt = vector(length = n)
  for (j in 1:n) {
    a2_error_pt[j] = spherical_dist_sq(5, error_pt, spherical_data[j,])^2
  }
  min_cost_error_pt = min_cost(proj_data, a2_error_pt, range = 5, iter = 1e+06)$par
  
  base_pt_dists = vector(length = n)
  Bengio_dists = vector(length = n)
  min_cost_dists = vector(length = n)
  
  for (i in 1:n) {
    base_pt_dists[i] = spherical_dist_sq(5, base_pt, spherical_data[i,])
    Bengio_dists[i] = dist_sq(Bengio_error_pt, proj_data[i,])
    min_cost_dists[i] = dist_sq(min_cost_error_pt, proj_data[i,])
  }
  
  cbind(abs(base_pt_dists - min_cost_dists)/base_pt_dists, 
        abs(base_pt_dists - Bengio_dists)/base_pt_dists)
}

sim = function(n = 20, b = 100, new_pt, sd) {
  theta = runif(n,0,2*pi)
  phi = runif(n, 0, pi)
  r = 5*(runif(n,0,1))^(1/3)
  
  spherical_data = cbind(theta,phi)
  
  x = r*sin(phi)*cos(theta)
  y = r*sin(phi)*sin(theta)
  z = r*cos(phi)
  
  cartesian_data = matrix(c(x,y,z), ncol = 3)
  proj_data = Bengio_data(cartesian_data, 2, K_mds)
  
  mean_errors = matrix(nrow = b, ncol = 2)
  for (i in 1:b) {
    errors = get_rel_errors(spherical_data, proj_data, new_pt, sd)
    mean_errors[i,] = colMeans(errors)
  }
  
  mean_errors
}

build_box = function(sds, base_pt, data_size = 20) {
  n = length(sds)
  mat = matrix(nrow = n*100, ncol = 3)
  
  for(i in 1:n) {
    mat[((i-1)*100+1):(i*100),1] = rep(sds[i],100)
    mat[((i-1)*100+1):(i*100),2:3] = sim(n = data_size, new_pt = base_pt, sd = sds[i])
  }
  table = data.frame(sd = c(mat[,1],mat[,1]),
                     mean_rel_error = c(mat[,2],mat[,3]),
                     type = c(rep("Min Cost", n*100), rep("Bengio", n*100)))
  
  ggplot(table, aes(x = factor(sd), y = mean_rel_error, color = type)) + geom_boxplot()
}

p1 = build_box(sds = c(pi/16,pi/8,pi/4), base_pt = c(0,pi/8))
p2 = build_box(sds = c(pi/16,pi/8,pi/4), base_pt = c(0,pi/4))
p3 = build_box(sds = c(pi/16,pi/8,pi/4), base_pt = c(0,pi/2))
grid.arrange(p1,p2,p3,ncol = 1)

################################################################################

sim = function(n, new_pt, error_sd) {
  theta = runif(n,0,2*pi)
  phi = runif(n, 0, pi)
  r = 5*(runif(n,0,1))^(1/3)
  
  spherical_data = cbind(theta,phi)
  
  x = r*sin(phi)*cos(theta)
  y = r*sin(phi)*sin(theta)
  z = r*cos(phi)
  
  cartesian_data = matrix(c(x,y,z), ncol = 3)
  proj_data = Bengio_data(cartesian_data, 2, K_mds)
  
  new_pt_cart = c(5*sin(new_pt[2]*cos(new_pt[1])), 
                  5*sin(new_pt[2]*sin(new_pt[1])),
                  5*cos(new_pt[2]))
  
  Bengio_base = Bengio(cartesian_data, new_pt_cart, 2, K_mds)
  
  error_theta = rnorm(1, 0, error_sd)
  error_phi = rnorm(1, 0, error_sd/4)
  error_pt = new_pt + c(error_theta, error_phi)
  error_pt_cart = c(5*sin(error_pt[2]*cos(error_pt[1])), 
                    5*sin(error_pt[2]*sin(error_pt[1])),
                    5*cos(error_pt[2]))
  
  Bengio_pt = Bengio(cartesian_data, error_pt_cart, 2, K_mds)
  
  a2 = vector(length = n)
  for (i in 1:n) {
    a2[i] = spherical_dist_sq(5, error_pt, spherical_data[1,])^2
  }
  min_cost_pt = min_cost(proj_data, a2, range = 5, iter = 1e+07)$par
  
  data = data.frame(x = c(proj_data[,1], Bengio_base[1], min_cost_pt[1], Bengio_pt[1]), 
                    y = c(proj_data[,2], Bengio_base[2], min_cost_pt[2], Bengio_pt[2]),
                    type = c(rep("Data", n), "Bengio Base", "Min Cost", "Bengio"))
  ggplot(data, aes(x = x, y = y, color = factor(type))) + geom_point()
}

sim(60, new_pt = c(0, pi/4), error_sd = pi/4)

