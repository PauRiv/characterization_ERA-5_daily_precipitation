#####################################################################
#####                                                           #####
#####              Kullback-Leibler divergence test             #####
#####                                                           #####
#####   (based on the Extended Generalized Pareto distribution  #####
#####               EGPD   Tencaliec et al. 2019)               #####
#####                                                           #####
#####             Author: Pauline Rivoire 09.06.2020            #####
#####                                                           #####
#####################################################################

# This file contains the functions Kdiv.EGP() and test.Kdiv()


# Note that the functions fitBBGP and dEGP.BB are EGPD functions available on request to the authors Tencaliec et al. 2019
#                                  (Flexible semiparametric generalized Pareto modeling of the entire range of rainfall amount)
# fitBBGP(X, m)         fits to the precipitation sample X an EGPD with Berstein polynomials of degree m
# dEGP.BB(X, theta)     gives the density of precip X for a given fitted EGPD parameters theta

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kdiv.EGP
# Description
# 
# Kullback divergence associated with the fitting of a EGPdistribution
# 
# Arguments
# X    a first time serie of positive precipitation
# Y    a second time serie of positive precipitation
# m    polynomial order of the Bernstein approximation  (default m = 30)
#
# Value
# Kdiv.EGP gives the Kullback divergence between X and Y, with densities of the type EGP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Kdiv.EGP <- function(X, Y, m=30){
  theta.f.bar <- fitBBGP(X, m)
  theta.g.bar <- fitBBGP(Y, m)
  
  if(!is.na(sum(theta.f.bar)) & !is.na(sum(theta.g.bar))){
    # apply the fitted densities to the two time series
    f.hat.X <- dEGP.BB(X, theta.f.bar)
    g.hat.X <- dEGP.BB(X, theta.g.bar)
    f.hat.Y <- dEGP.BB(Y, theta.f.bar)
    g.hat.Y <- dEGP.BB(Y, theta.g.bar)
    
    # compute the Kullback divergence, with the symmetrical part
    K <- sum(log(f.hat.X/g.hat.X)/length(X)) + sum(log(g.hat.Y/f.hat.Y)/length(Y))
  } else {
    K <- NA
  }#end if else NA
  
  return(K)
}#end function Kdiv.EGP




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test.Kdiv
# Description
# 
# Perform a two-sample test based on Kullback divergence
# 
# Arguments
# X       a first time serie of positive precipitation.
# Y       a second time serie of positive precipitation.
# m       polynomial order of the Bernstein approximation (default m = 30).
# n_K_h0  number of intermadiate values of K to compute to build the K distribution under the
#         null-hypothesis h0: f=g (default n_K_h0 = 500).
#
# Value
# K.stat        value of the test statistic.
# K.distrib.h0  vector (size n_K_h0) of values of K.stat used as distribution of K under h0.
# mean.K.h0     mean value of K.distrib.h0.
# p.value       the p-value of the bilateral test, P[K=K.stat|h0].
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test.Kdiv <- function(X, Y, m = 30, n_K_h0 = 500){
  
  K.stat <- Kdiv.EGP(X, Y, m = m)
  
  Z <- c(X,Y)
  K.distrib.h0 <- vector("numeric", length = n_K_h0)
  # K.distrib.h0 is the vector of Kullback divergences computed from several X_suffle
  #and Y_shuffle, which are randomly chosen in Z
  
  for (ii in 1:n_K_h0) {
    pi_x <- sample(1:length(Z), size = floor(length(Z)/2))
    pi_y <- setdiff(1:length(Z),pi_x)
    
    X_shuffle <- Z[pi_x]
    Y_shuffle <- Z[pi_y]
    K.distrib.h0[ii] <- Kdiv.EGP(X_shuffle,Y_shuffle)
  }#end for ii
  
  mean.K.h0 <- mean(K.distrib.h0, na.rm = T)
  if (K.stat>mean.K.h0){
    p.value <- (1 - ecdf(K.distrib.h0)(K.stat))
  } else {
    p.value <- "K.stat<mean.K.h0"
  }
  
  return(list(K.stat = K.stat, mean.K.h0 = mean.K.h0, p.value = p.value, K.distrib.h0=K.distrib.h0))
  
}#end function test.Kdiv