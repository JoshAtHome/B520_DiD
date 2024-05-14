rm(list = ls())
library(dplyr)
library(fixest)
library(did)
options(scipen = 999)
options(digits=4)

wd <- 'INPUT'
setwd(wd)


T = 7
N = 20
shock <- T

sim_data = function(heterogeneity, violation){
  dat = expand.grid(t = 1:T,i = 1:N) # create grid
  normal_random <- rnorm(T, mean = 0.5, sd = 0.1) 
  trend <- cumsum(normal_random) # create trend
  ##############################################################################
  dat$effect <- 1 # homogeneous effect
 # if (heterogeneity!=0){dat$effect <- rep(rnorm(N, mean = heterogeneity, sd = 1), each = T)} # heterogeneity across units
 # if (heterogeneity!=0){dat$effect <- rep(rgamma(N, shape = heterogeneity, rate = 1), each = T)}
  if (heterogeneity!=0){dat$effect <- rep(rnorm(N, mean = 1, sd = heterogeneity*1), each = T)} # heterogeneity across units
  dat <- mutate(dat, group = ifelse(i > N/2,"treat","control"), G = shock*(group == "treat"),
                treat = 1L*(group == "treat"), control = 1L*(group == "control"),  current = 1L*(G == t),
                pre = 1L*(t < T), post = 1L*(t >= T),
                tre = trend[t],
                random_change = rnorm(N*T, mean = 0, sd = 0.1), # include a tiny random change
                external_shock = ifelse(i > N/6,1,0), # shock for violations of parallel trends
                y = tre + treat * 0.2 + 
                  control * random_change + 
                  post * treat * effect + 
                  post * external_shock * violation + rnorm(N*T, mean = 0, sd = 0.001) , # outcome formula
                #################################################################
                X1 = tre + post * treat * violation + rnorm(N*T, mean = 0, sd = 0.2)) # covariate for testing
#  print(cor(dat$external_shock, dat$y))
  dat
}


##############################################################################


get_descriptive = function(data, treat){ # gets descriptive statistics for each dataset 
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
  dat <- sim_data(heterogeneity, violation) # simulate dataset
  descr0 <- get_descriptive(dat, 0) # descriptive for control
  descr1 <- get_descriptive(dat, 1) # descriptive for treatment
  att <- aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~1, data = dat, panel=TRUE)) # compute standard DiD
  att_cond <- aggte(att_gt(yname = "y", tname = "t", idname = "i", gname = "G", xformla = ~external_shock, data = dat, panel=TRUE)) # compute conditional DiD
  est <- att$overall.att # return estimate
  est_cond <- att_cond$overall.att # return conditional estimate
  se <- att$overall.se # return standard errors
  act <- mean(dat$effect[dat$treat == 1]) # return actual effect, because we know it and need it 
  return(list(est=est, est_cond=est_cond, act=act, se=se, descr0=descr0, descr1=descr1))
}


repl <- function(vec){
  imp <- vec  
  missing <- is.na(vec)
  imp[missing] <- mean(unlist(vec), na.rm = TRUE)
  return(imp)
} # helping function for na values in vectors

sim <- function(amount, heterogeneity, violation) {
  effect_est <- numeric(amount) # all estimated effects
  effect_est_cond <- numeric(amount)  # all conditional estimated effects                                      
  effect_act <- numeric(amount)    # all actual effects                                            
  standard_errors <- numeric(amount)   # all standard errors
  
  descr0_mat <- matrix(NA, nrow = amount, ncol = 4) # all descriptive control
  descr1_mat <- matrix(NA, nrow = amount, ncol = 4) # all descriptive treatment
  cat("Computing simulation....\n")
  
  for (i in 1:amount) {
    res <- att_estimate(heterogeneity, violation)  # compute simulation, get estimates and values                             
    effect_est[i] <- res$est
    effect_est_cond[i] <- res$est_cond
    effect_act[i] <- res$act
    standard_errors[i] <- res$se
    descr0_mat[i, ] <- res$descr0
    descr1_mat[i, ] <- res$descr1} # .... insert them into the matrices
  
  bias <- abs(mean(repl(effect_est) - repl(effect_act)))   # calculate difference between estimate and true effect
  bias_cond <- abs(mean(repl(effect_est_cond) - repl(effect_act)))   # calculate difference between conditional estimate and true effect
  se <- mean(standard_errors) # mean standard error
  descr0_means <- colMeans(descr0_mat, na.rm = TRUE) # compute means of descriptive matrices control
  descr1_means <- colMeans(descr1_mat, na.rm = TRUE) # compute means of descriptive matrices treatment
  
  return(list(bias=bias, bias_cond=bias_cond, se=se, descr0_means=descr0_means, descr1_means=descr1_means)) }

sim(1,10,2)

bias_hom <- function(n,viol) {
  result <- sim(n,1,viol)
  bias_table <- data.frame(Treatment_effect = c('homogeneous'),
                           Bias = round(result$bias, 5),
                           Bias_cond = round(result$bias_cond, 5))
  bias_table
} # testing function

bias_hom(10,1)


##############################################################################
##############################################################################

library(reshape2)
library(ggplot2)

bias_table <- function(n) {
  biases <- matrix(data=NA, nrow = 5, ncol = 5) # matrix biases
  biases_cond <- matrix(data=NA, nrow = 5, ncol = 5) # matrix biases
  standard_errors <- matrix(data=NA, nrow = 5, ncol = 5) # matrix standard errors
  desc0 <- matrix(data=NA, nrow = n^2, ncol = 4) # matrix descriptive control
  desc1 <- matrix(data=NA, nrow = n^2, ncol = 4) # matrix descriptive treatment
  counter <- 1
  for (i in 0:4){ # iterate through heterogeneities
    for (j in 0:4){ # iterate through violations
      result <- sim(n,i,j) # compute simulation
      desc0[counter,] <- result$descr0_means # import descriptive control
      desc1[counter,] <- result$descr1_means # import descriptive treatment
      biases[i+1,j+1] <- result$bias # import biases
      biases_cond[i+1,j+1] <- result$bias_cond # import conditional biases
      standard_errors[i+1,j+1] <- result$se # import standard errors
      counter <- counter + 1
    }
  }
  mean_se <- colMeans(standard_errors, na.rm = TRUE) # means standard errors
  mean_bias <- rowMeans(biases, na.rm = TRUE) # mean biases
  desc0 <- colMeans(desc0, na.rm = TRUE) 
  desc1 <- colMeans(desc1, na.rm = TRUE)
  desc <- rbind(desc0,desc1)
  return(list(biases=biases,biases_cond=biases_cond, mean_bias=mean_bias, mean_se=mean_se, desc=desc))
} 

desc_df <- function(desc) {
  desc_df <- data.frame(group =  c('Untreated', 'Treated'), M = desc[, 1], SE = desc[, 2], lower = desc[, 3], upper = desc[, 4])
  rownames(desc_df) <- NULL
  return(desc_df)
} # put descriptive into nice df

plot_se <- function(bia, st_e) {
  tab <- data.frame(matrix(ncol = length(st_e), nrow = 2))
  colnames(tab) <- c('homogeneity ', 'heterogeneity (1) ', 'heterogeneity (2) ','heterogeneity (3) ','heterogeneity (4) ')
  rownames(tab) <- c('Mean Bias', 'Mean Standard Error')
  tab[1, ] <- bia
  tab[2, ] <- st_e
  return(tab)
} # compare standard errors of biases

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
} # plot the biases returned from bias_table

tab_bias_cond <- function(bia) {
  tab <- data.frame(matrix(ncol = 5, nrow = 5))
  rownames(tab) <- c('homogeneity ', 'heterogeneity (1) ', 'heterogeneity (2) ','heterogeneity (3) ','heterogeneity (4)')
  colnames(tab) <- c('V0', 'V1', 'V2', 'V3', 'V4')
  tab[] <- bia
  return(tab)
} # put conditional biases into df

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

figure1 <- plot_heatmap(biases, 'Conditional Biases') # shows that heterogeneity does not increase bias
figure1

ggsave("figure1.png", plot = figure1, width = 4, height = 4, units = "in", bg='#ffffff')

tab3 <- tab_bias_cond(biases_cond) # biases of conditional did for comparison 
tab3

print(xtable(tab1), type = "latex")
print(xtable(tab2), type = "latex")
print(xtable(tab3), type = "latex")



