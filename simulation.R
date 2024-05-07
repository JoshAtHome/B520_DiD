library(knitr)

n<-1000

# sim(n, conditional, heterogeneity, violation) --- REAL ATT = 3
# @param conditional 0,1
# @param heterogeneity 0: effect = 3, x: effect = random, mean x
# @param violation 0: no violation, 1: mean violation of 0, x: mean violation of x

#################### no violation ##############################################
namev0 <- c('sim000','sim010','sim050','sim100','sim110','sim150')
sim000 <- sim(n,0,0,0) # standard / homogeneity
sim010 <- sim(n,0,1,0) # standard / mean heterogeneity = true effect
sim050 <- sim(n,0,5,0) # standard / mean heterogeneity = 5
sim100 <- sim(n,1,0,0) # conditional / homogeneity
sim110 <- sim(n,1,1,0) # conditional / mean heterogeneity = true effect
sim150 <- sim(n,1,5,0) # conditional / mean heterogeneity = 5

#################### violation mean 0 ##########################################
namev1 <- c('sim001','sim011','sim051','sim101','sim111','sim151')
sim001 <- sim(n,0,0,1) # standard / homogeneity
sim011 <- sim(n,0,1,1) # standard / mean heterogeneity = true effect
sim051 <- sim(n,0,5,1) # standard / mean heterogeneity = 5
sim101 <- sim(n,1,0,1) # conditional / homogeneity
sim111 <- sim(n,1,1,1) # conditional / mean heterogeneity = true effect
sim151 <- sim(n,1,5,1) # conditional / mean heterogeneity = 5

#################### violation mean 3 ##########################################
namev3 <- c('sim003','sim013','sim053','sim103','sim113','sim153')
sim003 <- sim(n,0,0,3) # standard / homogeneity
sim013 <- sim(n,0,1,3) # standard / mean heterogeneity = true effect
sim053 <- sim(n,0,5,3) # standard / mean heterogeneity = 5
sim103 <- sim(n,1,0,3) # conditional / homogeneity
sim113 <- sim(n,1,1,3) # conditional / mean heterogeneity = true effect
sim153 <- sim(n,1,5,3) # conditional / mean heterogeneity = 5


##  create descritptive statistics of outcome
create_desc_table <- function() {
  desc <- matrix(NA, nrow = 2, ncol = 4)
  desc[1,] <- sim000$descr0_means
  desc[2,] <- sim000$descr1_means
  desc_table <- data.frame(group = c('control', 'treatment'), mean = desc[, 1], lower = desc[, 2],
                           upper = desc[, 3], sd = desc[, 4])
  desc_table
}

