rm(list = ls())
library(dplyr)
library(fixest)
library(did)
options(scipen = 999)

wd <- '/Users/josh/Documents/Laptop/SS24/AdvancedTopics/codes2'
setwd(wd)


T = 7
N = 20
shock <- T

sim_data = function(heterogeneity, violation){
  dat = expand.grid(t = 1:T,i = 1:N)
  normal_random <- rnorm(T, mean = 0.5, sd = 0.1) 
  trend <- cumsum(normal_random)
  ##############################################################################
  dat$effect <- 1
  if (heterogeneity!=0){dat$effect <- rep(rnorm(N, mean = heterogeneity, sd = 1), each = T)} # heterogeneity across units
 # if (heterogeneity!=0){dat$effect <- rep(rnorm(N, mean = 1, sd = heterogeneity*0.1), each = T)} # heterogeneity across units
  
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
#  print(cor(dat$external_shock, dat$y))
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

##############################################################################
a <- att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat1, panel=TRUE)
ggdid(a)
a

get_descriptive = function(data, treat){
  treated_outcomes <-data$y[data$treat == treat & data$t==T]
  M <- mean(treated_outcomes)
  SD <- sd(treated_outcomes)
  lower <- t.test(treated_outcomes)$conf.int[1]
  upper <- t.test(treated_outcomes)$conf.int[2]
  all <- c(M, SD, lower, upper)
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

sim(1,10,2)

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
  desc0 <- matrix(data=NA, nrow = n^2, ncol = 4)
  desc1 <- matrix(data=NA, nrow = n^2, ncol = 4)
  ### create desc0 and desc1 for all combinations. Het 0 vio 0, Het 0 vio 1, Het 0 vio 2
  # all 25. 
  # if het i==0 and j==0 desc0
  counter <- 1
  for (i in 0:4){
    for (j in 0:4){
      result <- sim(n,i,j)
      desc0[counter,] <- result$descr0_means
      desc1[counter,] <- result$descr1_means
      biases[i+1,j+1] <- result$bias
      biases_cond[i+1,j+1] <- result$bias_cond
      standard_errors[i+1,j+1] <- result$se
      counter <- counter + 1
    }
  }
  mean_se <- colMeans(standard_errors, na.rm = TRUE)
  mean_bias <- rowMeans(biases, na.rm = TRUE)
  desc0 <- colMeans(desc0, na.rm = TRUE)
  desc1 <- colMeans(desc1, na.rm = TRUE)
  desc <- rbind(desc0,desc1)
  return(list(biases=biases,biases_cond=biases_cond, mean_bias=mean_bias, mean_se=mean_se, desc=desc))
}

desc_df <- function(desc) {
  desc_df <- data.frame(group =  c('Untreated', 'Treated'), M = desc[, 1], SE = desc[, 2], lower = desc[, 3], upper = desc[, 4])
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
    theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        aspect.ratio = 1) + 
     theme(panel.background = element_blank())
  return(map)
}



############################################################
data <- bias_table(5)
############################################################
biases <- data$biases
biases_cond <- data$biases_cond
mean_se <- data$mean_se
mean_bias <- data$mean_bias
descriptive_hom_noviol <- data$desc
############################################################


library(xtable)
# Data
tab1 <-desc_df(descriptive_hom_noviol) # shows descriptive statistics of homogeneous effects with no violation
tab1

# Results
tab2 <- plot_se(mean_bias, mean_se) # shows that heterogeneity does not increase bias, but does increase standard error
tab2

figure1 <- plot_heatmap(biases, 'Conditional Biases') # shows that heterogneity does not increase bias
figure1

ggsave("figure1.png", plot = figure1, width = 4, height = 4, units = "in", bg='#ffffff')

print(xtable(tab1), type = "latex")
print(xtable(tab2), type = "latex")


