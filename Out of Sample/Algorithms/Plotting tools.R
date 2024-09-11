library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)

create_plots2 = function(X1, risks1, X2, risks2, indices) {
  d1 = data.frame(x = X1[indices,1],
                  y = X1[indices,2],
                  trust = risks1[,1],
                  cont = risks1[,2],
                  prec = risks1[,3],
                  stress = risks1[,4],
                  dist_cor = risks1[,5])
  
  d2 = data.frame(x = X2[indices,1],
                  y = X2[indices,2],
                  trust = risks2[,1],
                  cont = risks2[,2],
                  prec = risks2[,3],
                  stress = risks2[,4],
                  dist_cor = risks2[,5])
  
  p1_trust = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = trust), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$trust, d2$trust)), 1)) + 
    ggtitle("Trustworthiness")
  
  p2_trust = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = trust), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$trust, d2$trust)), 1)) + 
    ggtitle("Trustworthiness")
  
  p1_cont = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = cont), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$cont, d2$cont)), 1)) + 
    ggtitle("Continuity")
  
  p2_cont = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = cont), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$cont, d2$cont)), 1)) + 
    ggtitle("Continuity")
  
  p1_prec = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = prec), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$prec, d2$prec)), 1)) + 
    ggtitle("Projected Precision Score")
  
  p2_prec = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = prec), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$prec, d2$prec)), 1)) + 
    ggtitle("Projected Precision Score")
  
  p1_stress = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = stress), size = 2) + 
    scale_color_viridis(limits = c(min(d1$stress), max(d1$stress))) + 
    ggtitle("Normalized Stress")
  
  p2_stress = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = stress), size = 2) + 
    scale_color_viridis(limits = c(min(d2$stress), max(d2$stress))) + 
    ggtitle("Normalized Stress")
  
  p1_cor = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = dist_cor), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$dist_cor, d2$dist_cor)), 1)) + 
    ggtitle("Correlation of Distances")
  
  p2_cor = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = dist_cor), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$dist_cor, d2$dist_cor)), 1)) + 
    ggtitle("Correlation of Distances")
  
  grid.arrange(p1_trust, p2_trust, 
               p1_cont, p2_cont, 
               p1_prec, p2_prec,
               p1_stress, p2_stress,
               p1_cor, p2_cor,
               ncol = 2)
}

create_plots3 = function(X1, risks1, X2, risks2, X3, risks3, indices) {
  d1 = data.frame(x = X1[indices,1],
                  y = X1[indices,2],
                  trust = risks1[,1],
                  cont = risks1[,2],
                  prec = risks1[,3],
                  stress = risks1[,4],
                  dist_cor = risks1[,5])
  
  d2 = data.frame(x = X2[indices,1],
                  y = X2[indices,2],
                  trust = risks2[,1],
                  cont = risks2[,2],
                  prec = risks2[,3],
                  stress = risks2[,4],
                  dist_cor = risks2[,5])
  
  d3 = data.frame(x = X3[indices,1],
                  y = X3[indices,2],
                  trust = risks3[,1],
                  cont = risks3[,2],
                  prec = risks3[,3],
                  stress = risks3[,4],
                  dist_cor = risks3[,5])
  
  p1_trust = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = trust), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$trust, d2$trust, d3$trust)), 1)) + 
    ggtitle("Trustworthiness")
  
  p2_trust = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = trust), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$trust, d2$trust, d3$trust)), 1)) + 
    ggtitle("Trustworthiness")
  
  p3_trust = d3 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = trust), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$trust, d2$trust, d3$trust)), 1)) + 
    ggtitle("Trustworthiness")
  
  p1_cont = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = cont), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$cont, d2$cont, d3$cont)), 1)) + 
    ggtitle("Continuity")
  
  p2_cont = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = cont), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$cont, d2$cont, d3$cont)), 1)) + 
    ggtitle("Continuity")
  
  p3_cont = d3 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = cont), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$cont, d2$cont, d3$cont)), 1)) + 
    ggtitle("Continuity")
  
  p1_prec = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = prec), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$prec, d2$prec, d3$prec)), 1)) + 
    ggtitle("Projected Precision Score")
  
  p2_prec = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = prec), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$prec, d2$prec, d3$prec)), 1)) + 
    ggtitle("Projected Precision Score")
  
  p3_prec = d3 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = prec), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$prec, d2$prec, d3$prec)), 1)) + 
    ggtitle("Projected Precision Score")
  
  p1_stress = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = stress), size = 2) + 
    scale_color_viridis(limits = c(min(d1$stress), max(d1$stress))) + 
    ggtitle("Normalized Stress")
  
  p2_stress = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = stress), size = 2) + 
    scale_color_viridis(limits = c(min(d2$stress), max(d2$stress))) + 
    ggtitle("Normalized Stress")
  
  p3_stress = d3 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = stress), size = 2) + 
    scale_color_viridis(limits = c(min(d3$stress), max(d3$stress))) + 
    ggtitle("Normalized Stress")
  
  p1_cor = d1 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = dist_cor), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$dist_cor, d2$dist_cor, d3$dist_cor)), 1)) + 
    ggtitle("Correlation of Distances")
  
  p2_cor = d2 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = dist_cor), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$dist_cor, d2$dist_cor, d3$dist_cor)), 1)) + 
    ggtitle("Correlation of Distances")
  
  p3_cor = d3 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = dist_cor), size = 2) + 
    scale_color_viridis(limits = c(min(c(d1$dist_cor, d2$dist_cor, d3$dist_cor)), 1)) + 
    ggtitle("Correlation of Distances")
  
  grid.arrange(p1_trust, p2_trust, p3_trust,
               p1_cont, p2_cont, p3_cont,
               p1_prec, p2_prec, p3_prec,
               p1_stress, p2_stress, p3_stress,
               p1_cor, p2_cor, p3_cor,
               ncol = 3)
}

create_plots = function(...) {
  if (length(list(...)) == 5) {
    create_plots2(...)
  }
  else if (length(list(...)) == 7) {
    create_plots3(...)
  }
  else {
    print("Incorrect number of arguments!")
  }
}

###############################################################################################

create_hist2 = function(risks1, risks2, types, num_per_type) {
  d = data.frame(trust = c(risks1[,1], risks2[,1]),
                 cont = c(risks1[,2], risks2[,2]),
                 prec = c(risks1[,3], risks2[,3]),
                 dist_cor = c(risks1[,5], risks2[,5]),
                 type = c(rep(types[1], num_per_type), rep(types[2], num_per_type)))
  
  p1 = ggplot(d, aes(x = trust)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3)
  
  p2 = ggplot(d, aes(x = cont)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3)
  
  p3 = ggplot(d, aes(x = prec)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3)
  
  p4 = ggplot(d, aes(x = dist_cor)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3)
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

create_hist3 = function(risks1, risks2, risks3, types, num_per_type) {
  d = data.frame(trust = c(risks1[,1], risks2[,1], risks3[,1]),
                 cont = c(risks1[,2], risks2[,2], risks3[,2]),
                 prec = c(risks1[,3], risks2[,3], risks3[,3]),
                 dist_cor = c(risks1[,5], risks2[,5], risks3[,5]),
                 type = c(rep(types[1], num_per_type), rep(types[2], num_per_type), rep(types[3], num_per_type)))
  
  p1 = ggplot(d, aes(x = trust)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3) +
    geom_density(data = subset(d, type == types[3]), color = "green", fill = "green", alpha = 0.3)
  
  p2 = ggplot(d, aes(x = cont)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3) +
    geom_density(data = subset(d, type == types[3]), color = "green", fill = "green", alpha = 0.3)
  
  p3 = ggplot(d, aes(x = prec)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3) +
    geom_density(data = subset(d, type == types[3]), color = "green", fill = "green", alpha = 0.3)
  
  p4 = ggplot(d, aes(x = dist_cor)) +
    geom_density(data = subset(d, type == types[1]), color = "red", fill = "red", alpha = 0.3) +
    geom_density(data = subset(d, type == types[2]), color = "blue", fill = "blue", alpha = 0.3) +
    geom_density(data = subset(d, type == types[3]), color = "green", fill = "green", alpha = 0.3)
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

create_hist = function(...) {
  if (length(list(...)) == 4) {
    create_hist2(...)
  }
  else if (length(list(...)) == 5) {
    create_hist3(...)
  }
  else {
    print("Incorrect number of arguments!")
  }
}