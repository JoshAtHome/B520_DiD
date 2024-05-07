
library(dplyr)
library(ggplot2)
library(did)

T = 7
N = 200
shock <- T-1

################################################################################
# @param heterogeneity decides if homogenous (0) or heterogenous (1) treatment effect 
# @param violation decides if a violation (1) or no violation (0) in parallel trends occurs
sim_data = function(heterogeneity, violation){
  dat = expand.grid(t = 1:T,i = 1:N) 
  normal_random <- rnorm(T, mean = 0.5, sd = 0.5)
  trend <- cumsum(normal_random)
  ################################################# create treatment effects
  dat$effect <- rnorm(N*T, mean = 3, sd = 3)
  if (heterogeneity!=1){dat$effect <- rnorm(N*T, mean = heterogeneity, sd = 3)}
  if (heterogeneity==0){dat$effect <- 3}
  ################################################# create violations
  dat$violation <- rnorm(N*T, mean = violation, sd = 3)
  if (violation==1){dat$violation <- rnorm(N*T, mean = 0, sd = 3)}
  if (violation==0){dat$violation <- 0}
  
  dat <- mutate(dat, group = ifelse(i > N/2,"treat","control"), G = shock*(group == "treat"),
                treat = 1L*(group == "treat"), control = 1L*(group == "control"),  current = 1L*(G == t), 
                pre = 1L*(t < T-1), post = 1L*(t >= T-1),
                pre_tre = pre*treat, pre_con = pre*control,
                post_tre = post*treat, post_con = post*control,
                tre = trend[t], X1 = rep(rnorm(N, mean = 3, sd = 1), each = T),
                y = tre +              # the basic trend
                    treat * 2 * X1 -   # the basic differences between treat / control
                    control * X1  +    # the basic differences between treat / control
                    violation * post_tre +          # the violation of parallel trends
                    current * effect)  # the treatment effect (homogeneous/heterogeneous)
  
  dat
}

################################################################################

show_graph = function(dat, label = "", show.means = TRUE) {
  gdat = dat %>% group_by(group, t, post, treat) %>% summarize(y = mean(y))
  gg = ggplot(gdat, aes(y = y, x = t, color = group)) + geom_line() +  geom_vline(xintercept = T-1) +
    theme_bw()
  gg
}
# data <- sim_data(0,0)
# show_graph(data)

################################################################################

get_descriptive = function(data, treat){
  treated_outcomes <-data$y[data$treat == treat]
  mean_outcome <- mean(treated_outcomes)
  lb <- t.test(treated_outcomes)$conf.int[1]
  ub <- t.test(treated_outcomes)$conf.int[2]
  sd <- sd(treated_outcomes)
  all <- c(mean_outcome, lb, ub, sd)
  return(all)
}

# get_descriptive(data, 0)

################################################################################

att_estimate <- function(heterogeneity, violation, x_formula){
  # @param heterogeneity 
  # @param violation
  # @param x_formula, either ~1 (standard) or ~X1 (conditional) or any other covariate for that matter
  dat <- sim_data(heterogeneity, violation)                                     # simulates the data accordingly
  descr0 <- get_descriptive(dat, 0)   
  descr1 <- get_descriptive(dat, 1)                                              
  att <- att_gt(yname = "y", tname = "t", idname = "i", gname = "G", 
                xformla = x_formula, data = dat)                                # estimates the att-values
  est <- att$att[5]                                                         # estimated treatment effect
  act <- mean(dat$effect[dat$treat == 1])                                   # actual treatment effect                                       # var-cov-matrix
  return(list(est=est, act=act, descr0=descr0, descr1=descr1))
}

################################################################################

sim <- function(amount, conditional, heterogeneity, violation) {
  # @param amount of times to simulate 
  # @param conditional (1) or standard (0)
  # @param heterogeneity 
  # @param violation
  effect_est <- numeric(amount)                                                 # all estimated effects
  effect_act <- numeric(amount)                                                 # all actual effects
  descr0_mat <- matrix(NA, nrow = amount, ncol = 4)
  descr1_mat <- matrix(NA, nrow = amount, ncol = 4)
  cat("Computing simulation....\n")
  if (conditional==0){for (i in 1:amount) {
      res <- att_estimate(heterogeneity, violation, ~1)                         # estimate under standard
      effect_est[i] <- res$est
      effect_act[i] <- res$act
      descr0_mat[i, ] <- res$descr0
      descr1_mat[i, ] <- res$descr1}}
  else if (conditional==1){for (i in 1:amount) {
      res <- att_estimate(heterogeneity, violation, ~X1)                        # estimate under conditional
      effect_est[i] <- res$est
      effect_act[i] <- res$act
      descr0_mat[i, ] <- res$descr0
      descr1_mat[i, ] <- res$descr1}}
  bias <- abs(mean(unlist(effect_est), na.rm = TRUE) -                          # calculate difference between estimate and true effect 
                mean(unlist(effect_act), na.rm = TRUE))
  descr0_means <- colMeans(descr0_mat, na.rm = TRUE)
  descr1_means <- colMeans(descr1_mat, na.rm = TRUE)
  
  return(list(bias=bias, descr0_means=descr0_means, descr1_means=descr1_means))
}

# sim(10,0,0,0)

################################################################################

create_bias_table <- function(names) {
  biases <- numeric(length(names))
  for (i in seq_along(names)) {biases[i] <- get(names[i])$bias}
  bias_table <- data.frame(DiD_Model = c('Standard', 'Standard', 'Standard', 'Conditional', 'Conditional', 'Conditional'),
                           Treatment_effect = c('homogeneous','heterogeneous (mean=3)','heterogeneous (mean=5)'),
                           Bias = biases)
  bias_table
}

################################################################################

create_desc_table_app <- function(names) {
  desc0 <- matrix(NA, nrow = length(names), ncol = 4)
  for (i in seq_along(names)) {desc0[i,] <- get(names[i])$descr0_means}
  desc0_tab <- data.frame(group = 'control', model = names, mean = desc0[, 1], lower = desc0[, 2], upper = desc0[, 3], sd = desc0[, 4])
  desc1 <- matrix(NA, nrow = length(names), ncol = 4)
  for (i in seq_along(names)) {desc1[i,] <- get(names[i])$descr1_means}
  desc1_tab <- data.frame(group = 'treatment', model = names, mean = desc1[, 1], lower = desc1[, 2], upper = desc1[, 3], sd = desc1[, 4])
  
  desc_app <- rbind(desc0_tab, desc1_tab)
  return(desc_app)
}

################################################################################






