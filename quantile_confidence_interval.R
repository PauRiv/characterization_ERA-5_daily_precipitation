#####################################################################
#####                                                           #####
#####             Confidence intervals of quantiles             #####
#####                                                           #####
#####   (based on the Extended Generalized Pareto distribution  #####
#####               EGPD   Tencaliec et al. 2019)               #####
#####                                                           #####
#####             Author: Pauline Rivoire 19.02.2020            #####
#####                                                           #####
#####################################################################



# This file contains the function Bootstrap_quantiles()


# Note that the functions fitBBGP and dEGP.BB are EGPD functions available on request to the authors Tencaliec et al. 2019
#                                  (Flexible semiparametric generalized Pareto modeling of the entire range of rainfall amount)
# fitBBGP(X, m)         fits to the precipitation sample X an EGPD with Berstein polynomials of degree m
# qEGP.BB(p, theta)     gives the quantiles of non-exceedance proba p for given fitted EGPD parameters theta


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bootstrap_quantiles
# Description
# 
# Bootstrapped values of the EGPD parameters and quantiles
# 
# Arguments
# season_precip    Time series of seasonal positive precipitation
# B                size of the bootstrap  (default B = 300)
# m                polynomial order of the Bernstein approximation  (default m = 30)
#
# Value
# Bootstrap_quantiles gives B boostrapped value of xi, sigma, w_1, ... w_m (the EGPD parameters) 
#                                               and of the quantiles for probabilities 0.1, 0.3, 0.5, 0.75, 0.9, 0.95
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bootstrap_quantiles <- function(season_precip, B=300, m=30){
  n_wetdays <- length(season_precip)
  
  indices_bootstrap <- sample.int(n_wetdays, size = (floor(2* n_wetdays / 3) + 1))
  ind_boot_number_1 <- sort(sample(indices_bootstrap, floor(length(indices_bootstrap) / 2)))
  ind_boot_number_2 <- setdiff(indices_bootstrap, ind_boot_number_1)
  
  boot_prec_mat_1 <- matrix(data = sample(season_precip[ind_boot_number_1],
                                          size = length(ind_boot_number_1)*B, replace = T),
                            nrow = length(ind_boot_number_1), ncol = B)
  
  boot_prec_mat_2 <- matrix(data = sample(season_precip[ind_boot_number_2],
                                          size = length(ind_boot_number_1)*B, replace = T),
                            nrow = length(ind_boot_number_1), ncol = B)
  
  boot_para <- cbind(apply(X = boot_prec_mat_1, FUN = fitBBGP, MARGIN = 2, m = m),
                     apply(X = boot_prec_mat_2, FUN = fitBBGP, MARGIN = 2, m = m))
  
  q_boot <- apply(X=boot_para, MARGIN = 2, FUN = qEGP.BB, x=c(0.1,0.3,0.5,0.75,0.9,0.95))
  row.names(RL_boot)=c("q0.1", "q0.3", "q0.5", "q0.75", "q0.9", "q0.95")
  
  return(rbind(boot_para, q_boot))
}

  
  
