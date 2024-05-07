rm(list = ls())
set.seed(182)
library(knitr)

################################################################################
################################################################################
wd <- 'REPLACE'
################################################################################
################################################################################
setwd(wd)

source("FINAL.R")  # T=7 * N=200 = 1400 observations per dataset
source("simulation.R") # n=1000

################################################################################
################################# Paper #########################################

desc01 <- create_desc_table() # standard model, homogeneity, no violation
tab0 <- create_bias_table(namev0) # no violation
tab1 <- create_bias_table(namev1) # violation mean 0
tab3 <- create_bias_table(namev3) # violation mean 3

desc01
tab0
tab1
tab3

# > tab0
# DiD_Model       Treatment_effect        Bias
# 1    Standard            homogeneous 0.000000000
# 2    Standard heterogeneous (mean=3) 0.002973617
# 3    Standard heterogeneous (mean=5) 0.007019343
# 4 Conditional            homogeneous 0.000000000
# 5 Conditional heterogeneous (mean=3) 0.027452656
# 6 Conditional heterogeneous (mean=5) 0.016141803
# > tab1
# DiD_Model       Treatment_effect       Bias
# 1    Standard            homogeneous 0.01389261
# 2    Standard heterogeneous (mean=3) 0.04798729
# 3    Standard heterogeneous (mean=5) 0.03343478
# 4 Conditional            homogeneous 0.02601735
# 5 Conditional heterogeneous (mean=3) 0.09713398
# 6 Conditional heterogeneous (mean=5) 0.05671153
# > tab3
# DiD_Model       Treatment_effect     Bias
# 1    Standard            homogeneous 3.014494
# 2    Standard heterogeneous (mean=3) 2.936238
# 3    Standard heterogeneous (mean=5) 2.983325
# 4 Conditional            homogeneous 2.974003
# 5 Conditional heterogeneous (mean=3) 3.019969
# 6 Conditional heterogeneous (mean=5) 3.000552

# > tab0
# DiD_Model       Treatment_effect        Bias
# 1    Standard            homogeneous 0.000000000
# 2    Standard heterogeneous (mean=3) 0.068357290
# 3    Standard heterogeneous (mean=5) 0.027624570
# 4 Conditional            homogeneous 0.000000000
# 5 Conditional heterogeneous (mean=3) 0.017695966
# 6 Conditional heterogeneous (mean=5) 0.001270277
# > tab1
# DiD_Model       Treatment_effect        Bias
# 1    Standard            homogeneous 0.016599659
# 2    Standard heterogeneous (mean=3) 0.040555909
# 3    Standard heterogeneous (mean=5) 0.001538203
# 4 Conditional            homogeneous 0.030233084
# 5 Conditional heterogeneous (mean=3) 0.087096667
# 6 Conditional heterogeneous (mean=5) 0.005147153
# > tab3
# DiD_Model       Treatment_effect     Bias
# 1    Standard            homogeneous 3.023842
# 2    Standard heterogeneous (mean=3) 2.953147
# 3    Standard heterogeneous (mean=5) 3.055548
# 4 Conditional            homogeneous 3.020872
# 5 Conditional heterogeneous (mean=3) 2.928136
# 6 Conditional heterogeneous (mean=5) 2.979227


################################################################################
############################### Appendix #######################################

desc0 <- create_desc_table_app(namev0) # no violation
desc1 <- create_desc_table_app(namev1) # violation mean 0
desc3 <- create_desc_table_app(namev3) # violation mean 3


# tab0 %>% #  'No violation', ifelse(vio == 1, 'Violation (mean=0)', ifelse(vio == 3, 'Violation (mean=3)', 0)))
#   kable(caption = "Biases for No Violation of the Parallel Trends Assumption",
#         format = "latex",
#         col.names = c("DiD_Model", "Treatment_effect", "Bias"),
#         align = "c")

