
# functions
source("./functions/sim_study_functions.R")

# libraries
library(parallel)
library(LaMa)


# Part 2: aggregate sizes -------------------------------------------------

# parameter for sim model
beta = matrix(c(log(c(8,7,6)), -0.2, 0.2, -0.6, 0.3, -0.2, -0.4), nrow = 3)
omega = matrix(c(0, 0.7, 0.3, 
                 0.2, 0, 0.8,
                 0.5, 0.5, 0), nrow = 3, byrow = TRUE)
stateparams = list(
  mu = c(20, 200, 800),
  sigma = c(20, 150, 500),
  kappa = c(0.2, 1, 2.5)
)

# trigonometric basis expansion
Z = cbind(1, LaMa::trigBasisExp(1:24, 24))
# calculating all dwell time means
dM = exp(Z%*%t(beta))

# get maximum mean for each state and calculate aggreagte sizes
maxlambda = apply(dM, 2, max)
agsizes = ceiling(qpois(0.995, maxlambda)+2)


# Simulation --------------------------------------------------------------

# one really long data set

set.seed(123)
data = sim_phsmm(1e6, beta, omega, stateparams)

## plotting state-dependent distributions
mu = stateparams$mu
sigma = stateparams$sigma
kappa = stateparams$kappa
delta = c(sum(data$C==1), sum(data$C==2), sum(data$C==3))/nrow(data)

pdf("./figures/simulation/state_dependent.pdf", width=8, height = 4)
par(mfrow = c(1,2), xpd = F)
hist(data$step, prob = T, breaks = 100, bor = "white", main = "", 
     ylim = c(0, 0.0025), xlim = c(0,2500), ylab = "density", xlab = "turning angle")
curve(delta[1]*dgamma(x, shape=mu[1]^2/sigma[1]^2, scale=sigma[1]^2/mu[1]),add=T, lwd=2, col = color[1], n = 500)
curve(delta[2]*dgamma(x, shape=mu[2]^2/sigma[2]^2, scale=sigma[2]^2/mu[2]),add=T, lwd=2, col = color[2], n = 500)
curve(delta[3]*dgamma(x, shape=mu[3]^2/sigma[3]^2, scale=sigma[3]^2/mu[3]),add=T, lwd=2, col = color[3], n = 500)
curve(delta[1]*dgamma(x, shape=mu[1]^2/sigma[1]^2, scale=sigma[1]^2/mu[1])+
        delta[2]*dgamma(x, shape=mu[2]^2/sigma[2]^2, scale=sigma[2]^2/mu[2])+
        delta[3]*dgamma(x, shape=mu[3]^2/sigma[3]^2, scale=sigma[3]^2/mu[3]), add=T, lwd=2, lty=2, n=500)
hist(data$angle, prob = T, breaks = 50, bor = "white", main = "", ylab = "density", xlab = "turning angle")
curve(delta[1]*CircStats::dvm(x, 0, kappa[1]), add=T, lwd=2, col = color[1], n = 500)
curve(delta[2]*CircStats::dvm(x, 0, kappa[2]), add=T, lwd=2, col = color[2], n = 500)
curve(delta[3]*CircStats::dvm(x, 0, kappa[3]), add=T, lwd=2, col = color[3], n = 500)
curve(delta[1]*CircStats::dvm(x, 0, kappa[1])+
        delta[2]*CircStats::dvm(x, 0, kappa[2])+
        delta[3]*CircStats::dvm(x, 0, kappa[3]), add=T, lwd=2, lty=2, n=500)
dev.off()


# Model fitting -----------------------------------------------------------

library(optimParallel)

## Model 1: true periodic HSMM

L = 24; K = 1
Z = LaMa::trigBasisExp(1:L-1, L=L, degree=K)
# indexshift is because in LaMa gamma^(t) = Pr(S_t | S_t-1) while in this paper gamma^(t) = Pr(S_t+1 | S_t)
ind = which(!is.na(data$step) & !is.na(data$angle))

thetainit = c(log(c(stateparams$mu,
                    stateparams$sigma,
                    stateparams$kappa)),
              as.numeric(beta), 
              rep(0,N*(N-2))
)

# cl = makeCluster(7); setDefaultCluster(cl=cl)
# modpHSMM_optim = optimParallel(par=thetainit, fn=mllkpHSMM, X=data, Z=Z, ind=ind, agsizes=agsizes,
#                          parallel=list(forward=TRUE),
#                          control = list(trace=4, maxit=1e4, ndeps=rep(1e-5,length(thetainit))))
# stopCluster(cl)
# 
# modpHSMM = nlm(mllkpHSMM, modpHSMM_optim$par, X=data, Z=Z, ind=ind, agsizes=agsizes,
#                iterlim = 500, print.level = 2)
# saveRDS(modpHSMM, "./simulation_study/models/modpHSMM.rds")

modpHSMM = readRDS("./simulation_study/models/modpHSMM.rds")


## Model 2: periodic HMM

Z = LaMa::trigBasisExp(1:L-1, L=L, degree=K)

thetainit = c(log(c(stateparams$mu,
                    stateparams$sigma,
                    stateparams$kappa)), 
              rep(-2,N*(N-1)), rep(0, 2*N*(N-1)*K))

# cl = makeCluster(7); setDefaultCluster(cl=cl)
# modpHMM_optim = optimParallel(par=thetainit, fn=mllkHMM, X=data, Z=Z, ind=ind, K=1,
#                                parallel=list(forward=T),
#                                control = list(trace=4, maxit=1e4, ndeps=rep(1e-5,length(thetainit))))
# stopCluster(cl)
# 
# modpHMM = nlm(mllkHMM, modpHMM_optim$par, X=data, Z=Z, ind=ind, K=1,
#                iterlim = 500, print.level = 2)
# saveRDS(modpHMM, "./simulation_study/models/modpHMM.rds")

modpHMM = readRDS("./simulation_study/models/modpHMM.rds")

## Model 3: homogeneous HSMM

lambda = c(8,7,6) # exponentiated intercepts as starting value

thetainit = c(log(c(stateparams$mu,
                    stateparams$sigma,
                    stateparams$kappa,
                    lambda)), 
              rep(0,N*(N-2))
)

# cl = makeCluster(7); setDefaultCluster(cl=cl)
# modHSMM_optim = optimParallel(par=thetainit, fn=mllkHSMM, X=data, ind=ind, agsizes=agsizes,
#                                parallel=list(forward=T),
#                                control = list(trace=4, maxit=1e4, ndeps=rep(1e-5,length(thetainit))))
# stopCluster(cl)
# modHSMM = nlm(mllkHSMM, modHSMM_optim$par, X=data, ind=ind, agsizes=agsizes,
#               iterlim = 500, print.level = 2, stepmax = 10)
# saveRDS(modHSMM, "./simulation_study/models/modHSMM.rds")

modHSMM = readRDS("./simulation_study/models/modHSMM.rds")


# Transform parameters ----------------------------------------------------
# pHSMM
N = 3; K = 1; L = 24
modpHSMM$mu = exp(modpHSMM$estimate[1:N])
modpHSMM$sigma = exp(modpHSMM$estimate[N+1:N])
modpHSMM$kappa = exp(modpHSMM$estimate[2*N + 1:N])
modpHSMM$beta = matrix(modpHSMM$estimate[3*N + 1:(N*(1+2*K))], nrow = N, ncol = 1+2*K)
Z = LaMa::trigBasisExp(1:L-1, L=L, degree=K)
modpHSMM$Lambda = exp(cbind(1,Z)%*%t(modpHSMM$beta))
if(N>2){ # only needed if N>2
  omega = matrix(0,N,N)
  omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(modpHSMM$estimate[3*N+N*(1+2*K)+1:(N*(N-2))])),N,N-1)))
  omega = t(omega)/apply(omega,2,sum)
}else{ omega = matrix(c(0,1,1,0),2,2) }
modpHSMM$omega = omega
dm = list()
for(j in 1:N){ dm[[j]] = sapply(1:agsizes[j]-1, dpois, lambda = modpHSMM$Lambda[,j]) }
modpHSMM$dm = dm
modpHSMM$Gamma = LaMa::tpm_phsmm(omega, dm)
modpHSMM$delta = LaMa::stationary_p(modpHSMM$Gamma, t = data$tod[1])

ind = which(!is.na(data$step) & !is.na(data$angle))
modpHSMM$allprobs = matrix(NA, nrow = 1e6, ncol = N)
for(j in 1:N){
  modpHSMM$allprobs[ind,j] = dgamma(data$step[ind], shape=modpHSMM$mu[j]^2/modpHSMM$sigma[j]^2, scale=modpHSMM$sigma[j]^2/modpHSMM$mu[j])*
    CircStats::dvm(data$angle[ind], mu = 0, kappa = modpHSMM$kappa[j])
}
allprobs_large = t(apply(modpHSMM$allprobs, 1, rep, times = agsizes))
states = LaMa::viterbi_p(modpHSMM$delta, modpHSMM$Gamma, allprobs_large, data$tod)
modpHSMM$rawstates = states
states = rep(NA, 1e6)
for(j in 1:N){
  states[which(modpHSMM$rawstates %in% c(c(0,cumsum(agsizes)[-N])[j]+1:agsizes[j]))] = j
}
modpHSMM$states = states
modpHSMM$dec_accuracy = sum(modpHSMM$states==data$C)/1e6


## pHMM
modpHMM$mu = exp(modpHMM$estimate[1:N])
modpHMM$sigma = exp(modpHMM$estimate[N+1:N])
modpHMM$kappa = exp(modpHMM$estimate[2*N + 1:N])
Z = LaMa::trigBasisExp(1:L, L=L, degree=K)
modpHMM$beta = matrix(modpHMM$estimate[3*N + 1:(N*(N-1)*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
modpHMM$Gamma = LaMa::tpm_p(tod=1:L, L=L, beta=modpHMM$beta, degree=K, Z=Z)
modpHMM$delta = LaMa::stationary_p(modpHMM$Gamma, t = data$tod[1])

ind = which(!is.na(data$step) & !is.na(data$angle))
allprobs = matrix(NA, nrow = 1e6, ncol = N)
for(j in 1:N){
  allprobs[ind,j] = dgamma(data$step[ind], shape=modpHMM$mu[j]^2/modpHMM$sigma[j]^2, scale=modpHMM$sigma[j]^2/modpHMM$mu[j])*
    CircStats::dvm(data$angle[ind], mu = 0, kappa = modpHMM$kappa[j])
}
modpHMM$allprobs = allprobs
modpHMM$states = LaMa::viterbi_p(modpHMM$delta, modpHMM$Gamma, modpHMM$allprobs, data$tod)
modpHMM$dec_accuracy = sum(modpHMM$states==data$C)/1e6

## HSMM
modHSMM$mu = exp(modHSMM$estimate[1:N])
modHSMM$sigma = exp(modHSMM$estimate[N+1:N])
modHSMM$kappa = exp(modHSMM$estimate[2*N + 1:N])
modHSMM$lambda = exp(modHSMM$estimate[3*N + 1:N])
if(N>2){ # only needed if N>2
  omega = matrix(0,N,N)
  omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N),exp(modHSMM$estimate[4*N+1:(N*(N-2))])),N,N-1)))
  omega = t(omega)/apply(omega,2,sum)
}else{ omega = matrix(c(0,1,1,0),2,2) }
modHSMM$omega = omega
for(j in 1:N){ dm[[j]] = dpois(1:agsizes[j]-1, modHSMM$lambda[j]) }
modHSMM$dm = dm
modHSMM$Gamma = LaMa::tpm_hsmm(omega, dm)
modHSMM$delta = LaMa::stationary(modHSMM$Gamma)

ind = which(!is.na(data$step) & !is.na(data$angle))
modHSMM$allprobs = matrix(NA, nrow = 1e6, ncol = N)
for(j in 1:N){
  modHSMM$allprobs[ind,j] = dgamma(data$step[ind], shape=modHSMM$mu[j]^2/modHSMM$sigma[j]^2, scale=modHSMM$sigma[j]^2/modHSMM$mu[j])*
    CircStats::dvm(data$angle[ind], mu = 0, kappa = modHSMM$kappa[j])
}
allprobs_large = t(apply(modHSMM$allprobs, 1, rep, times = agsizes))
modHSMM$rawstates = LaMa::viterbi(modHSMM$delta, modHSMM$Gamma, allprobs_large)
states = rep(NA, 1e6)
for(j in 1:N){ states[which(modHSMM$rawstates %in% c(c(0,cumsum(agsizes)[-N])[j]+1:agsizes[j]))] = j }
modHSMM$states = states
modHSMM$dec_accuracy = sum(modHSMM$states==data$C)/1e6



# Decoding accuracies -----------------------------------------------------

modnames = c("periodic HSMM", "periodic HMM", "homogeneous HSMM")

accuracies = round(c(modpHSMM$dec_accuracy, modpHMM$dec_accuracy, modHSMM$dec_accuracy), 4)
names(accuracies) = modnames
accuracies*100

# Compare dwell-time distributions ----------------------------------------

# function to conveniently compute empirical dwell-time distribution of decoded states


# functions needed for overall dwell-time distribution
source("./functions/dwell_functions.R")

## pHSMM
# overall dwell-time distribution
pmf = list()
for(j in 1:N){
  pmf[[j]] = ddwell_hsmm(1:agsizes[j], j, modpHSMM$dm, modpHSMM$Gamma, agsizes)
}
modpHSMM$pmf = pmf
# empirical dwell-time distributions of vitebi decoded states
modpHSMM$emppmf = empirical_pmf(modpHSMM$states, agsizes)

## pHMM
pmf = list()
for(j in 1:N){
  pmf[[j]] = ddwell(1:agsizes[j], j, modpHMM$Gamma)
}
modpHMM$pmf = pmf
# empirical dwell-time distributions of vitebi decoded states
modpHMM$emppmf = empirical_pmf(modpHMM$states, agsizes)

## HSMM
modHSMM$pmf = modHSMM$dm
# empirical dwell-time distributions of vitebi decoded states
modHSMM$emppmf = empirical_pmf(modHSMM$states, agsizes)

saveRDS(modpHSMM, "./simulation_study/models/full/modpHSMM.rds")
saveRDS(modpHMM, "./simulation_study/models/full/modpHMM.rds")
saveRDS(modHSMM, "./simulation_study/models/full/modHSMM.rds")

modpHSMM = readRDS("./simulation_study/models/full/modpHSMM.rds")
modpHMM = readRDS("./simulation_study/models/full/modpHMM.rds")
modHSMM = readRDS("./simulation_study/models/full/modHSMM.rds")


# comparison
mods = list(modpHSMM, modpHMM, modHSMM)

# plotting model-implied overall dwell-time distributions against viterbi decoded

# pdf("./figures/simulation/dwell-time_distribution.pdf", width = 7, height = 5)
# 
# color = c("orange", "deepskyblue", "seagreen2")
# statenames = c("resting", "foraging", "travelling")
# modnames = c("periodic HSMM", "periodic HMM", "homogeneous HSMM")
# 
# par(mfrow = c(3,3), mar = c(2.5,3,2,0.5)+0.1)
# for(mod in 1:3){
#   if(mod == 3){
#     par(mar = c(5,4,3,0.5)+0.1)
#   }
#   for(j in 1:N){
#     if(mod==1){
#       main = statenames[j]
#     } else{main = ""}
#     if(mod == 3){
#       xlab = "dwell time (hours)"
#     } else{ xlab = ""}
#     if(j == 1){
#       ylab = "probabilities"
#       par(mar = c(4,4,2,0.5)+0.1)
#     } else{ylab = ""
#       mar = c(2.5,3,2,0.5)+0.1
#     }
#     
#     plot(1:agsizes[j], mods[[mod]]$pmf[[j]], type = "h", lwd = 2, col = color[j],
#          ylim = c(0, 0.2),
#          xlab = xlab, ylab = ylab, bty = "n", main = main)
#     points(1:agsizes[j], mods[[mod]]$emppmf[[j]], pch = 19, col = "#00000030")
#     if(j == 1){
#       legend("topright", legend = modnames[mod], bty = "n")
#     }
#   }
#   legend("topright", pch = c(19, NA), lwd = c(NA,2), col = c("#00000030", "gray"),
#          legend = c("viterbi-decoded","model-implied"), bty = "n")
# }
# 
# dev.off()



# plotting model-implied overall dwell-time distributions against viterbi decoded

pdf("./figures/simulation/dwell-time_distribution.pdf", width = 7, height = 5)

N=3

m = matrix(c(1,1,1,2:10), nrow = 4, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(0.3, 1, 1, 1))
par(mar = c(0.5,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", xlim = c(0,10), ylim = c(1,10))

legend(x = 1.7, y = 9.5, inset = c(0.3,0), pch = c(NA, NA, 19), lwd = c(2, 1, NA), 
       lty = c(1, 1, NA), col = c("gray", "#00000030", "#00000030"),
       legend = c("model-implied", "true empirical", "Viterbi-decoded"), bty = "n", horiz = TRUE,
       text.width = c(2, 1.8, 0))


color = c("orange", "deepskyblue", "seagreen2")
statenames = c("state 1", "state 2", "state 3")

emppmf_true = empirical_dwell_times(data$C, agsizes)

par(mar = c(5.5,4,2,1))

for(mod in 1:3){
  if(mod == 3){
    par(mar = c(5,4,3,0.5)+0.1)
  }
  for(j in 1:N){
    if(mod==1){
      main = statenames[j]
    } else{main = ""}
    if(mod == 3){
      xlab = "dwell time (hours)"
    } else{ xlab = ""}
    if(j == 1){
      ylab = "probabilities"
      par(mar = c(4,4,2,0.5)+0.1)
    } else{ylab = ""
    mar = c(2.5,3,2,0.5)+0.1
    }
    
    plot(1:agsizes[j], mods[[mod]]$pmf[[j]], type = "h", lwd = 2, col = color[j],
         ylim = c(0, 0.2), xlim = c(0,20),
         xlab = xlab, ylab = ylab, bty = "n", main = main)
    points(1:agsizes[j], mods[[mod]]$emppmf[[j]], pch = 19, col = "#00000030")
    lines(1:agsizes[j], emppmf_true[[j]], lwd = 1, lty = 1, col = "#00000030")
    
    if(j == 3){
      legend("topright", legend = modnames[mod], bty = "n")
    }
  }
}

dev.off()

# Pseudo resiudals --------------------------------------------------------

# functions to calculate cdf transform
pres_step = function(step, mu, sigma, stateprobs){
  N = ncol(stateprobs)
  n = nrow(stateprobs)
  pres = rep(0, n)
  shape = mu^2/sigma^2; scale = sigma^2/mu
  for(j in 1:N){
    pres = pres + stateprobs[,j] * pgamma(step, shape=shape[j], scale=scale[j])
  }
  qnorm(pres)
}
pres_angle = function(angle, kappa, stateprobs){
  N = ncol(stateprobs)
  n = nrow(stateprobs)
  pres = rep(0, n)
  for(j in 1:N){
    pres = pres + stateprobs[,j] * CircStats::pvm(angle, mu=0, kappa = kappa[j])
  }
  qnorm(pres)
}

# pHSMM
stateprobs = matrix(NA, 1e6, 3)
cumagsizes = c(0, cumsum(agsizes[-N]))
# calculating state probs split because of memory
for(i in 1:100){
  cat("\n",i)
  allprobs_large = t(apply(modpHSMM$allprobs[1e4*(i-1)+1:1e4,], 1, rep, times = agsizes))
  stateprobs_large = LaMa::stateprobs_p(modpHSMM$delta, modpHSMM$Gamma, allprobs_large, tod = data$tod[1e4*(i-1)+1:1e4])
  for(j in 1:N){
    stateprobs[1e4*(i-1)+1:1e4, j] = rowSums(stateprobs_large[,cumagsizes[j]+1:agsizes[j]])
  }
}
modpHSMM$stateprobs = stateprobs
modpHSMM$pseudores_step = pres_step(data$step, modpHSMM$mu, modpHSMM$sigma, modpHSMM$stateprobs)

# pHMM
modpHMM$stateprobs = LaMa::stateprobs_p(modpHMM$delta, modpHMM$Gamma, modpHMM$allprobs, tod = data$tod)
modpHMM$pseudores_step = pres_step(data$step, modpHMM$mu, modpHMM$sigma, modpHMM$stateprobs)

# HSMM
stateprobs = matrix(NA, 1e6, 3)
cumagsizes = c(0, cumsum(agsizes[-N]))
# calculating state probs split because of memory
for(i in 1:100){
  cat("\n",i)
  allprobs_large = t(apply(modHSMM$allprobs[1e4*(i-1)+1:1e4,], 1, rep, times = agsizes))
  stateprobs_large = LaMa::stateprobs(modHSMM$delta, modHSMM$Gamma, allprobs_large)
  for(j in 1:N){
    stateprobs[1e4*(i-1)+1:1e4, j] = rowSums(stateprobs_large[,cumagsizes[j]+1:agsizes[j]])
  }
}
modHSMM$stateprobs = stateprobs
modHSMM$pseudores_step = pres_step(data$step, modHSMM$mu, modHSMM$sigma, modHSMM$stateprobs)


## plotting
mods = list(modpHSMM, modpHMM, modHSMM)
par(mfrow = c(2,3))
for(m in 1:3){
  hist(mods[[mod]]$pseudores_step, prob = TRUE, xlab = "pseudo residual (step length)",
       ylab = "density", main = modnames[m], bor = "white", ylim = c(0,0.5), xlim = c(-3,3))
  curve(dnorm(x), add = TRUE, col = "orange", lwd = 2, lty = 2)
}
for(m in 1:3){
  acf(mods[[mod]]$pseudores_step, bty = "n", xlim = c(0,5), main = "")
}



