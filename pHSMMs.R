# https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study2736765655


elephant = read.csv("/Users/jan-ole/Desktop/periodic inference ISEC/ivory_elephant1.csv")

# Preprocessing -----------------------------------------------------------


library(tidyverse)
elephant$time = as.character(elephant$timestamp) %>% 
  str_sub(start = 12, end = 13)
elephant$time = as.integer(elephant$time)

library(moveHMM)
library(optimParallel)
library(LaMa)

elephant = moveHMM::prepData(elephant, coordNames = c("location.long", "location.lat"))
head(elephant)

elephant$step = elephant$step*1000

boxplot(elephant$step~elephant$time, bor = "white")

# mllk_hsmm = function(theta.star, X, N, agsizes){
#   mu_dwell = exp(theta.star[1:N])
#   mu = exp(theta.star[N + 1:N])
#   sigma = exp(theta.star[2*N + 1:N])
#   kappa = exp(theta.star[3*N + 1:N])
#   size = exp(theta.star[4*N + 1:N])
#   if(N>2){ # only needed if N>2
#   omega = matrix(0,N,N)
#   omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(theta.star[(1+2*degree)*N+2*N+1:(N*(N-2))])),N,N-1)))
#   omega = t(omega)/apply(omega,2,sum)
#   } else{ omega = matrix(c(0,1,1,0),2,2) }
#   dm = vector("list")
#   for(j in 1:2){ dm[[j]] = dnbinom((1:agsizes[j])-1, mu = mu_dwell[j], size = size[j]) }
#   Gamma = LaMa::tpm_hsmm(omega, dm)
#   delta = LaMa::stationary(Gamma)
#   allprobs = matrix(1, nrow(X), N)
#   ind = which(!is.na(X$step) & !is.na(X$angle))
#   for (j in 1:N){
#     allprobs[ind,j] = dgamma(X$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])*
#       CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
#   }
#   -LaMa::forward_s(delta, Gamma, allprobs, agsizes)
# }
# 
# theta.star = c(log(c(
#   10, 7, 100, 1500, 100, 800, 0.5, 2, 100, 100
# )))
#   
# N = 2
# mod1 = nlm(mllk_hsmm, theta.star, X = elephant, N = 2, agsizes = c(14,14), print.level = 2, iterlim = 1e3)  
# 
# agsizes = c(14,14)
# mu_dwell_hat1 = exp(mod1$estimate[1:N])
# mu_hat1 = exp(mod1$estimate[N + 1:N])
# sigma_hat1 = exp(mod1$estimate[2*N + 1:N])
# kappa_hat1 = exp(mod1$estimate[3*N + 1:N])
# size_hat1 = exp(mod1$estimate[4*N + 1:N])
# if(N>2){ # only needed if N>2
#   omega_hat1 = matrix(0,N,N)
#   omega_hat1[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(mod1$estimate[5*N+1:(N*(N-2))])),N,N-1)))
#   omega_hat1 = t(omega_hat1)/apply(omega_hat1,2,sum)
# } else{ omega_hat1 = matrix(c(0,1,1,0),2,2) }
# dm_hat1 = vector("list")
# for(j in 1:2){ dm_hat1[[j]] = dnbinom((1:agsizes[j])-1, mu = mu_dwell_hat1[j], size = size_hat1[j]) }
# Gamma_hat1 = LaMa::tpm_hsmm(omega_hat1, dm_hat1)
# delta_hat1 = LaMa::stationary(Gamma_hat1)
# 
# color = c("orange", "deepskyblue")
# par(mfrow = c(1,2))
# plot(dm_hat1[[1]], type = "h", lwd = 2, col = color[1], bty = "n", 
#      xlab = "dwell time (hours)", ylab = "probabilities", xaxt = "n")
# axis(1, at = seq(0,14, by = 2), labels = seq(0,28,by=4))
# 
# plot(dm_hat1[[2]], type = "h", lwd = 2, col = color[2], bty = "n", 
#      xlab = "dwell time (hours)", ylab = "probabilities", xaxt = "n")
# axis(1, at = seq(0,14, by = 2), labels = seq(0,28,by=4))


## inhomogeneous model
L = 12
elephant$tod = elephant$time/2 +1
Z = LaMa::trigBasisExp(tod = 0:11, L=12, degree = 1)

mllk_phsmm = function(theta.star, X, Z, N, K = 1, agsizes){
  beta = matrix(theta.star[1:((1+2*K)*N)], nrow = N, ncol = 1+2*K)
  Mu_dwell = exp(cbind(1,Z)%*%t(beta))
  mu = exp(theta.star[(1+2*K)*N + 1:N])
  sigma = exp(theta.star[(1+2*K)*N + N + 1:N])
  kappa = exp(theta.star[(1+2*K)*N + 2*N + 1:N])
  size = exp(theta.star[(1+2*K)*N + 3*N + 1:N])
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N), exp(theta.star[(1+2*K)*N+4*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  dm = vector("list")
  for(j in 1:2){ 
    dm[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu = Mu_dwell[,j], size = size[j])
  }
  Gamma = LaMa::tpm_phsmm(omega, dm)
  delta = LaMa::stationary_p(Gamma, X$tod[1])
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for (j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  allprobs = t(apply(allprobs, 1, rep, times = agsizes))
  -LaMa::forward_p(delta, Gamma, allprobs, tod = X$tod)
}

K = 1
theta.star = c(log(c(5,4)), rep(0,2*K*N), log(c(
  300, 1100, 300, 800, 0.2, 1, 100, 100
)))


cl = makeCluster(7); setDefaultCluster(cl=cl)
mod2 = optimParallel(theta.star, mllk_phsmm, X = elephant, N = 2, Z = Z, K = 1, agsizes = c(20,20),
                     control = list(trace = 4, maxit = 1e4, ndeps=rep(1e-6,length(theta.star))))
stopCluster(cl)

mod2.1 = nlm(mllk_phsmm, mod2$par, X = elephant, Z = Z, N = 2, agsizes = c(20,20),
           print.level = 2, iterlim = 1e3, stepmax = 10)

K = 1; N = 2
beta_hat2 = matrix(mod2.1$estimate[1:((1+2*K)*N)], nrow = N, ncol = 1+2*K)
Mu_dwell_hat2 = exp(cbind(1,Z)%*%t(beta_hat2))
mu_hat2 = exp(mod2.1$estimate[(1+2*K)*N + 1:N])
sigma_hat2 = exp(mod2.1$estimate[(1+2*K)*N + N + 1:N])
kappa_hat2 = exp(mod2.1$estimate[(1+2*K)*N + 2*N + 1:N])
size_hat2 = exp(mod2.1$estimate[(1+2*K)*N + 3*N + 1:N])
omega_hat2 = matrix(c(0,1,1,0),2,2)
dm_hat2 = vector("list")
for(j in 1:2){ 
  dm_hat2[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu = Mu_dwell_hat2[,j], size = size_hat2[j])
}
Gamma_hat2 = LaMa::tpm_phsmm(omega_hat2, dm_hat2)
delta_hat2 = LaMa::stationary_p(Gamma_hat2, elephant$tod[1])

plot(Mu_dwell_hat2[,1], pch = 19, col = color[1], bty = "n")
plot(Mu_dwell_hat2[,2], pch = 19, col = color[2], bty = "n")

source("./dwell_functions.R")

pmf1 = ddwell_hsmm(1:14, 1, dm_hat2, Gamma_hat2, agsizes)
pmf2 = ddwell_hsmm(1:14, 2, dm_hat2, Gamma_hat2, agsizes)

plot(pmf1)
plot(pmf2)

# more flexible modelling of size parameter


mllk_phsmm2 = function(theta.star, X, Z, N, K = 1, agsizes){
  beta = array(theta.star[1:((1+2*K)*N*2)], dim = c(N, 1+2*K, 2))
  Mu_dwell = exp(cbind(1,Z)%*%t(beta[,,1]))
  Size = exp(cbind(1,Z)%*%t(beta[,,2]))
  mu = exp(theta.star[(1+2*K)*N*2 + 1:N])
  sigma = exp(theta.star[(1+2*K)*N*2 + N + 1:N])
  kappa = exp(theta.star[(1+2*K)*N*2 + 2*N + 1:N])
  if(N>2){ # only needed if N>2
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N), exp(theta.star[(1+2*K)*N+4*N+1:(N*(N-2))])),N,N-1)))
    omega = t(omega)/apply(omega,2,sum)
  } else{ omega = matrix(c(0,1,1,0),2,2) }
  dm = vector("list")
  for(j in 1:2){ 
    dm[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu = Mu_dwell[,j], size = Size[,j])
  }
  Gamma = LaMa::tpm_phsmm(omega, dm)
  delta = LaMa::stationary_p(Gamma, X$tod[1])
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for (j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  allprobs = t(apply(allprobs, 1, rep, times = agsizes))
  -LaMa::forward_p(delta, Gamma, allprobs, tod = X$tod)
}


K = 1; N = 2
theta.star = c(log(c(5,4)), 1.5, -1.5, -0.2, -0.6, log(c(10000,10000)), rep(0,2*K*N),
               log(c(
  300, 1100, 300, 800, 0.2, 1
)))

cl = makeCluster(7); setDefaultCluster(cl=cl)
mod3.1 = optimParallel(mod3$estimate, mllk_phsmm2, X = elephant, N = 2, Z = Z, K = 1, agsizes = c(35,35),
                     control = list(trace = 4, maxit = 1e4, ndeps=rep(1e-8,length(theta.star))))
stopCluster(cl)

mod3 = nlm(mllk_phsmm2, theta.star, X = X, Z = Z, N = 2, K = 1, agsizes = c(35,35),
             print.level = 2, iterlim = 1e3, stepmax = 5)

K = 1; N = 2; agsizes = c(35,35)
beta_hat3 = array(mod3$estimate[1:((1+2*K)*N*2)], dim = c(N, 1+2*K, 2))
Mu_dwell_hat3 = exp(cbind(1,Z)%*%t(beta_hat3[,,1]))
Size_hat3 = exp(cbind(1,Z)%*%t(beta_hat3[,,2]))
mu_hat3 = exp(mod3$estimate[(1+2*K)*N*2 + 1:N])
sigma_hat3 = exp(mod3$estimate[(1+2*K)*N*2 + N + 1:N])
kappa_hat3 = exp(mod3$estimate[(1+2*K)*N*2 + 2*N + 1:N])
omega_hat3 = matrix(c(0,1,1,0),2,2)
dm_hat3 = vector("list")
for(j in 1:2){ 
  dm_hat3[[j]] = sapply(1:agsizes[j]-1, dnbinom, mu = Mu_dwell_hat3[,j], size = Size_hat3[,j])
}
Gamma_hat3 = LaMa::tpm_phsmm(omega_hat3, dm_hat3)
Delta_hat3 = LaMa::stationary_p(Gamma_hat3)

pmf1 = ddwell_hsmm(1:20, 1, dm_hat3, Gamma_hat3, agsizes)
pmf2 = ddwell_hsmm(1:20, 2, dm_hat3, Gamma_hat3, agsizes)

plot(pmf1, type = "h")
plot(pmf2, type = "h")
