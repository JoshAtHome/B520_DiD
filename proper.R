rm(list = ls())
library(dplyr)
library(fixest)
library(did)
options(scipen = 999)

T = 7
N = 20
shock <- T

sim_data = function(heterogeneity, violation){
  dat = expand.grid(t = 1:T,i = 1:N)
  normal_random <- rnorm(T, mean = 0.5, sd = 0.1) 
  trend <- cumsum(normal_random)
  ##############################################################################
  dat$effect <- 1
  if (heterogeneity!=0){dat$effect <- rep(rnorm(N, mean = 1, sd = heterogeneity*0.1), each = T)} # heterogeneity across units
  if (heterogeneity!=0){dat$effect <- rep(rgamma(N, shape = heterogeneity, rate = 1), each = T)}
  dat <- mutate(dat, group = ifelse(i > N/2,"treat","control"), G = shock*(group == "treat"),
                treat = 1L*(group == "treat"), control = 1L*(group == "control"),  current = 1L*(G == t),
                pre = 1L*(t < T), post = 1L*(t >= T),
                tre = trend[t],
                random_change = rnorm(N*T, mean = 0, sd = 0.1),
                external_shock = ifelse(i > N/6,1,0),
                y = tre + treat * 0.2 + 
                  control * random_change + 
                  post * treat * effect + 
                  post * external_shock * violation + rnorm(N*T, mean = 0, sd = 0.001) ,
                #################################################################
                X1 = tre + post * treat * violation + rnorm(N*T, mean = 0, sd = 0.2))
  print(cor(dat$external_shock, dat$y))
  dat
}

##############################################################################
v <- 0
dat1 <- sim_data(1, v)
dat2 <- sim_data(2, v)
dat3 <- sim_data(3, v)
dat4 <- sim_data(4, v)

att1 <- aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat1, panel=TRUE), type = 'simple')$overall.att
att2 <- aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat2, panel=TRUE), type = 'simple')$overall.att

abs(att1 - mean(dat1$effect[dat1$current == 1]))
abs(att2 - mean(dat2$effect[dat2$current == 1]))

aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat1, panel=TRUE))
aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat2, panel=TRUE))





##############################################################################
a <- att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat1, panel=TRUE)
ggdid(a)
a

get_descriptive = function(data, treat){
  treated_outcomes <-data$y[data$treat == treat & data$t==T]
  mean_outcome <- mean(treated_outcomes)
  lb <- t.test(treated_outcomes)$conf.int[1]
  ub <- t.test(treated_outcomes)$conf.int[2]
  sd <- sd(treated_outcomes)
  all <- c(mean_outcome, lb, ub, sd)
  return(all)
}

dat <- sim_data(0,0)
get_descriptive(dat,1)


att_estimate <- function(heterogeneity, violation){
  dat <- sim_data(heterogeneity, violation)
  descr0 <- get_descriptive(dat, 0)
  descr1 <- get_descriptive(dat, 1)
  att <- aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat, panel=TRUE))
  att_cond <- aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~external_shock, data = dat, panel=TRUE))
  est <- att$overall.att
  est_cond <- att_cond$overall.att
  se <- att$overall.se
  act <- mean(dat$effect[dat$treat == 1])
  return(list(est=est, est_cond=est_cond, act=act, se=se, descr0=descr0, descr1=descr1))
}

repl <- function(vec){
  imp <- vec  
  missing <- is.na(vec)
  imp[missing] <- mean(unlist(vec), na.rm = TRUE)
  return(imp)
}

sim <- function(amount, heterogeneity, violation) {
  effect_est <- numeric(amount)
  effect_est_cond <- numeric(amount)                                            # all estimated effects
  effect_act <- numeric(amount)                                                 # all actual effects
  standard_errors <- numeric(amount)  
  
  descr0_mat <- matrix(NA, nrow = amount, ncol = 4)
  descr1_mat <- matrix(NA, nrow = amount, ncol = 4)
  cat("Computing simulation....\n")
  
  for (i in 1:amount) {
    res <- att_estimate(heterogeneity, violation)                               # estimate under standard
    effect_est[i] <- res$est
    effect_est_cond[i] <- res$est_cond
    effect_act[i] <- res$act
    standard_errors[i] <- res$se
    descr0_mat[i, ] <- res$descr0
    descr1_mat[i, ] <- res$descr1}
  
  bias <- abs(mean(repl(effect_est) - repl(effect_act)))                        # calculate difference between estimate and true effect
  bias_cond <- abs(mean(repl(effect_est_cond) - repl(effect_act)))              # calculate difference between estimate and true effect
  se <- mean(standard_errors)
  descr0_means <- colMeans(descr0_mat, na.rm = TRUE)
  descr1_means <- colMeans(descr1_mat, na.rm = TRUE)
  
  return(list(bias=bias, bias_cond=bias_cond, se=se, descr0_means=descr0_means, descr1_means=descr1_means)) }

sim(1,10,1)$bias

bias_hom <- function(n,viol) {
  result <- sim(n,1,viol)
  bias_table <- data.frame(Treatment_effect = c('homogeneous'),
                           Bias = round(result$bias, 5),
                           Bias_cond = round(result$bias_cond, 5))
  bias_table
}

bias_hom(10,1)


##############################################################################
##############################################################################

library(reshape2)
library(ggplot2)

bias_table <- function(n) {
  biases <- matrix(data=NA, nrow = 5, ncol = 5)
  biases_cond <- matrix(data=NA, nrow = 5, ncol = 5)
  standard_errors <- matrix(data=NA, nrow = 5, ncol = 5)
  for (i in 0:4){
    for (j in 0:4){
      result <- sim(n,i,j)
      if (i==0 && j==0){
        desc0 <- result$descr0_means
        desc1 <- result$descr1_means
        desc <- rbind(desc0, desc1)
      }
      biases[i+1,j+1] <- result$bias
      biases_cond[i+1,j+1] <- result$bias_cond
      standard_errors[i+1,j+1] <- result$se
    }
  }
  mean_se <- colMeans(standard_errors, na.rm = TRUE)
  mean_bias <- rowMeans(biases, na.rm = TRUE)
  return(list(biases=biases,biases_cond=biases_cond, mean_bias=mean_bias, mean_se=mean_se, desc=desc))
}

desc_df <- function(desc) {
  desc_df <- data.frame(group =  c('Untreated', 'Treated'), outcome = desc[, 1], lb = desc[, 2], ub = desc[, 3], se = desc[, 4])
  rownames(desc_df) <- NULL
  return(desc_df)
}

plot_se <- function(bia, st_e) {
  tab <- data.frame(matrix(ncol = length(st_e), nrow = 2))
  colnames(tab) <- c('homogeneity ', 'heterogeneity (1) ', 'heterogeneity (2) ','heterogeneity (3) ','heterogeneity (4) ')
  rownames(tab) <- c('Mean Bias', 'Mean Standard Error')
  tab[1, ] <- bia
  tab[2, ] <- st_e
  return(tab)
}

plot_heatmap <- function(mat, title) {
  mat_df <- as.data.frame(mat)
  colnames(mat_df) <- c("V0", "V1", "V2", 'V3', 'V4')
  mat_df$row <- 0:4
  mat_df <- melt(mat_df, id.vars = "row", variable.name = "column", value.name = "value")
  map <- ggplot(mat_df, aes(x = column, y = row, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "green", high = "red", name = "Bias") +
    labs(x = "Violation", y = "Heterogeneity", title = 'Unconditional Biases') +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, vjust=-2, size=20),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        aspect.ratio = 1)
  print(map)
}

####################
data <- bias_table(5)
####################

biases <- data$biases
biases_cond <- data$biases_cond
mean_se <- data$mean_se
mean_bias <- data$mean_bias
descriptive_hom_noviol <- data$desc


# Data
desc_df(descriptive_hom_noviol) # shows descriptive statistics of homogeneous effects with no violation
# Results
plot_se(mean_bias, mean_se) # shows that heterogeneity does not increase bias, but does increase standard error
plot_heatmap(biases, 'Conditional Biases') # shows that heterogneity does not increase bias


