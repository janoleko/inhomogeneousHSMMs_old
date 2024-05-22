# https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study2736765655

setwd("/Users/jan-ole/Desktop/periodic inference ISEC/")
elephant = read.csv("/Users/jan-ole/Desktop/periodic inference ISEC/data/ivory_elephant1.csv")

# Preprocessing -----------------------------------------------------------


library(tidyverse)
elephant$time = as.character(elephant$timestamp) %>% 
  str_sub(start = 12, end = 13)
elephant$time = as.integer(elephant$time)

library(moveHMM)

elephant = moveHMM::prepData(elephant, coordNames = c("location.long", "location.lat"))
head(elephant)

elephant$step = elephant$step*1000

boxplot(elephant$step~elephant$time)


# EDA ---------------------------------------------------------------------

pdf("./figures/track.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
plot(elephant$x[1:5000], elephant$y[1:5000], col = "#00000040", type = "l", bty = "n", xlab = "longitude", ylab = "latitude")
plot(elephant$step[1:200], type = "h", lwd = 0.5, bty = "n", xlab = "time", ylab = "step length")
abline(h=500, col = "orange", lty = 2)
dev.off()


# Homogeneous model -------------------------------------------------------

mllk = function(theta.star, X, N){
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  mu = exp(theta.star[N*(N-1)+1:N])
  sigma = exp(theta.star[N*(N-1)+N+1:N])
  kappa = exp(theta.star[N*(N-1)+2*N+1:N])
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step)&!is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape = mu[j]^2/sigma[j]^2, scale = sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], 0, kappa[j])
  }
  -LaMa::forward(delta, Gamma, allprobs)
}

hist(elephant$step)

theta.star = c(rep(-1,2), log(10), log(2000), log(10), log(1000), log(1), log(1))
N = 2

X = elephant %>% dplyr::select(step, angle)

mod1 = nlm(mllk, theta.star, N = 2, X = X, iterlim = 1000, print.level = 2)

library(LaMa)
N = 2
Gamma_hat = tpm(mod1$estimate[1:2])
delta_hat = stationary(Gamma_hat)
mu_hat = exp(mod1$estimate[N*(N-1)+1:N])
sigma_hat = exp(mod1$estimate[N*(N-1)+N+1:N])
kappa_hat = exp(mod1$estimate[N*(N-1)+2*N+1:N])

color = c("orange", "deepskyblue")

pdf("./figures/marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,1))
hist(elephant$step, breaks = 100, prob = T, bor = "white", xlim = c(0,5000), main = "", xlab = "step length", ylab = "density")
curve(delta_hat[1]*dgamma(x, shape = mu_hat[1]^2/sigma_hat[1]^2, scale = sigma_hat[1]^2/mu_hat[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_hat[2]*dgamma(x, shape = mu_hat[2]^2/sigma_hat[2]^2, scale = sigma_hat[2]^2/mu_hat[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_hat[1]*dgamma(x, shape = mu_hat[1]^2/sigma_hat[1]^2, scale = sigma_hat[1]^2/mu_hat[1])+
        delta_hat[2]*dgamma(x, shape = mu_hat[2]^2/sigma_hat[2]^2, scale = sigma_hat[2]^2/mu_hat[2]), add = T, lwd = 2, lty = 2, n = 500)
# dev.off()

hist(elephant$angle, breaks = 30, prob = T, bor = "white", main = "", xlab = "turning angle", ylab = "density")
curve(delta_hat[1]*CircStats::dvm(x, 0, kappa_hat[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_hat[2]*CircStats::dvm(x, 0, kappa_hat[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_hat[1]*CircStats::dvm(x, 0, kappa_hat[1])+
        delta_hat[2]*CircStats::dvm(x, 0, kappa_hat[2]), add = T, lwd = 2, lty = 2, n = 500)

dev.off()

## decoding states
allprobs = matrix(1, nrow(elephant), 2)
ind = which(!is.na(elephant$step)&!is.na(elephant$angle))
for(j in 1:N){
  allprobs[ind,j] = dgamma(elephant$step[ind], shape = mu_hat[j]^2/sigma_hat[j]^2, scale = sigma_hat[j]^2/mu_hat[j])*
    CircStats::dvm(elephant$angle[ind], 0, kappa_hat[j])
}

states = viterbi(delta_hat, Gamma_hat, allprobs)

pdf("./figures/decoded.pdf", width = 8, height = 4)
par(mfrow = c(1,2))

plot(elephant$x[1:500], elephant$y[1:500], type = "l", bty = "n", xlab = "longitude", ylab = "latitude", col = "white")
segments(x0 = elephant$x[1:499], y0 = elephant$y[1:499], 
         x1 = elephant$x[2:500], y1 = elephant$y[2:500], col = color[states], lwd = 1)

plot(elephant$step[201:400], type = "h", lwd = 0.5, col = color[states[201:400]], bty = "n", xlab = "time", ylab = "step length")

dev.off()

## dwell-time distribution

pdf("./figures/dwell_distr_hom.pdf", height = 3, width = 8)
par(mfrow = c(1,2))
plot(1:12, dgeom(0:11, prob = Gamma_hat[1,2]), type = "h", lwd = 2, bty = "n", xlim = c(0,12),
     col = color[1], ylab = "probabilities", xlab = "dwell time (hours)", xaxt = "n", main = "inactive state")
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
plot(1:12, dgeom(0:11, prob = Gamma_hat[2,1]), type = "h", lwd = 2, bty = "n", xlim = c(0,12),
     col = color[2], ylab = "probabilities", xlab = "dwell time (hours)", xaxt = "n", main = "active state")
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
dev.off()

# Inhomogeneous model -----------------------------------------------------

mllk2 = function(theta.star, X, N, K){
  beta = matrix(theta.star[1:((N*(N-1))*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
  Gamma = LaMa::tpm_p(tod = 1:12*2, L = 24, beta = beta, degree = K)
  delta = LaMa::stationary_p(Gamma, t = X$time[1])
  mu = exp(theta.star[(N*(N-1))*(1+2*K)+1:N])
  sigma = exp(theta.star[(N*(N-1))*(1+2*K)+N+1:N])
  kappa = exp(theta.star[(N*(N-1))*(1+2*K)+2*N+1:N])
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step)&!is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape = mu[j]^2/sigma[j]^2, scale = sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], 0, kappa[j])
  }
  -LaMa::forward_p(delta, Gamma, allprobs, X$time/2+1)
}


theta.star = c(rep(-1,2), rep(0,4), log(10), log(2000), log(10), log(1000), log(1), log(1))
N = 2; K = 1

X = elephant %>% dplyr::select(step, angle, time)

mod2 = nlm(mllk2, theta.star, X = X, N = 2, K = 1, iterlim = 1000, print.level = 2, hessian = T)

K = 1
beta_hat2 = matrix(mod2$estimate[1:((N*(N-1))*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
Gamma_hat2 = tpm_p(tod = 1:12*2, L = 24, beta = beta_hat2, degree = K)
Delta_hat2 = stationary_p(Gamma_hat2)
mu_hat2 = exp(mod2$estimate[(N*(N-1))*(1+2*K)+1:N])
sigma_hat2 = exp(mod2$estimate[(N*(N-1))*(1+2*K)+N+1:N])
kappa_hat2 = exp(mod2$estimate[(N*(N-1))*(1+2*K)+2*N+1:N])

# build continuous Delta
get_delta_cont = function(theta.star, t, state = c(1:2), N=2, L=12, K=1){
  coef = matrix(theta.star[1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  timeseq = (t + 1:L) %% L
  Gamma = array(NA, dim = c(N,N,L))
  for (k in 1:L){
    eta = iHSMM::pv(coef, time = timeseq[k], degree = K, L = L)
    G = diag(N)
    G[!G] = exp(eta)
    G = G/rowSums(G)
    Gamma[,,k] = G }
  GammaT=Gamma[,,1]
  for(k in 2:(L-1)){
    GammaT = GammaT%*%Gamma[,,k] }
  delta = solve(t(diag(N)-GammaT+1), rep(1,N))
  return(delta[state])
}

n = 200; N = 2
Delta_cont = matrix(nrow = n, ncol = N)
timeseq = seq(0,12,length = n)

for(i in 1:n){
  Delta_cont[i,] = get_delta_cont(mod2$estimate, t = timeseq[i])
}

# confidence intervals
L = 12
I = solve(mod2$hessian)
thetas = mvtnorm::rmvnorm(1000, mod2$estimate, sigma = I)
Deltas = array(dim = c(n, N, 1000))

for(t in 1:n){
  print(t)
  for(i in 1:1000){
    Deltas[t,,i] = get_delta_cont(thetas[i,], t = timeseq[t])
  }
}
DeltaCI = array(dim = c(n,N,2))
for(t in 1:n){
  for(j in 1:N){
    DeltaCI[t,j,] = quantile(Deltas[t,j,], probs = c(0.025, 0.975))
  }
}

Z = trigBasisExp(timeseq, L = 12, degree = 1)
Gamma_cont = tpm_p(L = 12, beta = beta_hat2, degree = 1, Z = Z)

rho = matrix(nrow = n, ncol = N)
for(i in 1:n){
  rho[i,] = stationary(Gamma_cont[,,i])
}

pdf("./figures/stationary_p1.pdf", width = 7, height = 4)
par(mfrow = c(1,1))
plot(timeseq, Delta_cont[,2], type = "l", bty = "n", lwd = 2, ylim = c(0,1), 
     ylab = "Pr(active)", col = color[2], xlab = "time of day", xaxt = "n")
points(Delta_hat2[,2], pch = 19, col = color[2])
lines(timeseq, rho[,2], lwd = 2, col = "gold2")
polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,2,1], rev(DeltaCI[,2,2])), col = alpha(color[2], 0.2), border = FALSE)

axis(1, at = seq(0,12,by = 2), labels = seq(0,24,by=4))
legend(x = 3, y = 1, bty = "n", lwd = 2, col = c(color[2], "gold2"), legend = c("stationary distribution", "biased approximation"))

dev.off()

source("./functions/dwell_functions.R")
source("./functions/auxiliary_functions.R")

dpmf1 = ddwell(1:12, 1, Gamma_hat2)
dpmf2 = ddwell(1:12, 2, Gamma_hat2)

par(mfrow = c(1,2))
plot(1:12, dpmf1, type = "h", lwd = 2, bty = "n", col = color[1], 
     xaxt = "n", xlab = "dwell time", ylab = "probabilities", ylim = c(0, 0.35))
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
plot(1:12, dpmf2, type = "h", lwd = 2, bty = "n", col = color[2], 
     xaxt = "n", xlab = "dwell time", ylab = "probabilities", ylim = c(0, 0.25))
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))


# Gamma plot

pdf("./figures/tpm_p.pdf", width = 7, height = 4)
par(mfrow = c(2,2), mar = c(3,4,2,1))
for(i in 1:N){
  for(j in 1:N){
    plot(Gamma_hat2[i,j,], type = "l", lwd = 2, ylim = c(0,1), bty = "n", ylab = paste0("Pr(",i,"->",j,")"), xaxt = "n", xlab = "")
    points(Gamma_hat2[i,j,], pch = 19)
    axis(1, at = seq(0,12, by = 2), labels = seq(0,24,by = 4))
  }
}
dev.off()

## decoding states
N=2
allprobs2 = matrix(1, nrow(X), N)
ind = which(!is.na(elephant$step)&!is.na(elephant$angle))
for(j in 1:N){
  allprobs2[ind,j] = dgamma(elephant$step[ind], shape = mu_hat2[j]^2/sigma_hat2[j]^2, scale = sigma_hat2[j]^2/mu_hat2[j])*
    CircStats::dvm(elephant$angle[ind], 0, kappa_hat2[j])
}

states2 = viterbi_p(Delta_hat2[elephant$time[1]/2+1,], Gamma_hat2, allprobs2, elephant$time/2+1)

dwellfreq = list()
for(j in 1:2){ 
  dwelltimes = get_dwell_times(states2, j)
  dwellfreq[[j]] = prop.table(table(factor(dwelltimes, levels = 1:max(dwelltimes))))
}

dpmf1 = ddwell(1:12, 1, Gamma_hat2)
dpmf2 = ddwell(1:12, 2, Gamma_hat2)

# dwelltimes homogeneous model
dwellfreq_h = list()
for(j in 1:2){ 
  dwelltimes_h = get_dwell_times(states, j)
  dwellfreq_h[[j]] = prop.table(table(factor(dwelltimes_h, levels = 1:max(dwelltimes_h))))
}

pdf("./figures/dwell_times_modelchecking.pdf", width = 8, height = 4)
par(mfrow = c(1,2))

plot(1:12, dgeom(1:12-1, prob = 1-Gamma_hat[1,1]), type = "h", lwd = 2, bty = "n", col = color[1], 
     xaxt = "n", xlab = "dwell time", ylab = "probabilities", ylim = c(0, 0.35), xlim = c(0,12), main = "no periodic variation")
points(1:12, dwellfreq_h[[1]][1:12], col = "#00000040", pch = 16)
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
# legend("top", bty = "n", legend = c("no periodic variation"))


plot(1:12, dpmf1, type = "h", lwd = 2, bty = "n", col = color[1],
     xaxt = "n", xlab = "dwell time", ylab = "probabilities", ylim = c(0, 0.35), xlim = c(0,12), main = "periodic variation")
points(1:10, dwellfreq[[1]], col = "#00000040", pch = 16)
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
# legend("top", bty = "n", legend = c("periodic variation"))

legend("topright", pch=c(NA,16), lwd = c(2, NA), col = c(color[1], "#00000040"), 
       legend = c("model-implied", "empirical"), bty = "n")

dev.off()

# plot(1:12, dpmf2, type = "h", lwd = 2, bty = "n", col = color[2], 
#      xaxt = "n", xlab = "dwell time", ylab = "probabilities", ylim = c(0, 0.25))
# points(1:12, dwellfreq[[2]][1:12], col = "#00000040", pch = 20)
# axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))


# time-varying dwell-time distribution

dpmft = array(dim = c(L, 20, 2))
for(j in 1:2){
  for(t in 1:L){
    dpmft[t,,j] = ddwell_t(1:20, t, j, Gamma_hat2)
  }
}

pdf("./figures/time_varying.pdf", width = 9, height = 5.5)
par(mfrow = c(4,3), mar = c(4,2,1,1))
j = 2
for(t in 1:12){
  plot(1:12, dpmft[t,,j][1:12], type = "h", lwd = 2, bty = "n", col = color[j],
       xaxt = "n", xlab = "dwell time", ylab = "probabilities", ylim = c(0, 0.35), xlim = c(0,12))
  axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
  legend("topright", bty = "n", legend = c(paste0(t*2-2,":00")))
}
dev.off()

# mean time-varying dwell times

Md = matrix(nrow = L, ncol = 2)
for(j in 1:2){
  for(t in 1:12){
    Md[t,j] = sum(1:50 * ddwell_t(1:50, t, j, Gamma_hat2))
  }
}

Mds = array(dim = c(L,2,1000))
K = 1; N = 2
for(j in 1:2){
  for(t in 1:12){
    for(i in 1:1000){
      beta = matrix(thetas[i,1:6], nrow = 2, ncol = 3)
      Gamma = tpm_p(tod = 1:12*2, L = 24, beta = beta, degree = K)
      Mds[t,j,i] = sum(1:20 * ddwell_t(1:20, t, j, Gamma))
    }
  }
}

MdCI = array(dim = c(L,N,2))
for(t in 1:L){
  for(j in 1:N){
    MdCI[t,j,] = quantile(Mds[t,j,], probs = c(0.025, 0.975))
  }
}

pdf("./figures/mean_dwell_times.pdf", width=7, height = 4)
par(mfrow = c(1,1), mar = c(5,4,4,2)+0.1)
plot(0:11, Md[,1], bty = "n", type = "l", lwd = 2, col = color[1], ylim = c(0,10),
     xaxt = "n", yaxt = "n", xlim = c(0,12), xlab = "time of day", ylab = "mean dwell time (hours)")
points(0:11, Md[,1], pch = 19, col = color[1])
segments(x0 = 0:11, x1 = 0:11, y0 = MdCI[,1,1], y1 = MdCI[,1,2], col = alpha(color[1], 0.8), lwd = 2)
segments(x0 = 0:11-0.1, x1 = 0:11+0.1, y0 = MdCI[,1,1], y1 = MdCI[,1,1], col = alpha(color[1], 0.8), lwd = 2)
segments(x0 = 0:11-0.1, x1 = 0:11+0.1, y0 = MdCI[,1,2], y1 = MdCI[,1,2], col = alpha(color[1], 0.8), lwd = 2)

# polygon(c(0:11, rev(0:11)), c(MdCI[,1,1], rev(MdCI[,1,2])), col = alpha(color[1], 0.2), border = FALSE)

lines(0:11, Md[,2], lwd = 2, col = color[2])
points(0:11, Md[,2], pch = 19, col = color[2])
# polygon(c(0:11, rev(0:11)), c(MdCI[,2,1], rev(MdCI[,2,2])), col = alpha(color[2], 0.2), border = FALSE)
segments(x0 = 0:11, x1 = 0:11, y0 = MdCI[,2,1], y1 = MdCI[,2,2], col = alpha(color[2], 0.8), lwd = 2)
segments(x0 = 0:11-0.1, x1 = 0:11+0.1, y0 = MdCI[,2,1], y1 = MdCI[,2,1], col = alpha(color[2], 0.8), lwd = 2)
segments(x0 = 0:11-0.1, x1 = 0:11+0.1, y0 = MdCI[,2,2], y1 = MdCI[,2,2], col = alpha(color[2], 0.8), lwd = 2)

axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
axis(2, at = seq(0,10, by = 2), labels = seq(0,20, by = 4))

legend("top", lwd = 2, col = color, legend = c("inactive state", "active state"), bty = "n")

dev.off()



## K = 2

theta.star = c(rep(-1,2), rep(0,8), log(300), log(1000), log(300), log(800), log(0.2), log(1))
N = 2; K = 2

X = elephant %>% dplyr::select(step, angle, time)

mod3 = nlm(mllk2, theta.star, X = X, N = 2, K = K, iterlim = 1000, print.level = 2, hessian = T)

beta_hat3 = matrix(mod3$estimate[1:((N*(N-1))*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
Gamma_hat3 = tpm_p(tod = 1:12*2, L = 24, beta = beta_hat3, degree = K)
Delta_hat3 = stationary_p(Gamma_hat3)
mu_hat3 = exp(mod3$estimate[(N*(N-1))*(1+2*K)+1:N])
sigma_hat3 = exp(mod3$estimate[(N*(N-1))*(1+2*K)+N+1:N])
kappa_hat3 = exp(mod3$estimate[(N*(N-1))*(1+2*K)+2*N+1:N])

n = 200; N = 2
Delta_cont3 = matrix(nrow = n, ncol = N)
timeseq = seq(0,12,length = n)

get_delta_cont2 = function(theta.star, t, state = c(1:2), N=2, L=12, K=1){
  coef = matrix(theta.star[1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  timeseq = (t + 1:L) %% L
  Z = trigBasisExp(timeseq, L = 12, degree = K)
  Gamma = tpm_p(L = 12, beta = coef, degree = K, Z = Z)
  GammaT=Gamma[,,1]
  for(k in 2:(L-1)){
    GammaT = GammaT%*%Gamma[,,k] }
  delta = solve(t(diag(N)-GammaT+1), rep(1,N))
  return(delta[state])
}

for(i in 1:n){ Delta_cont3[i,] = get_delta_cont2(mod3$estimate, t = timeseq[i], K=2) }

# confidence intervals
L = 12
I = solve(mod3$hessian)
thetas = mvtnorm::rmvnorm(2000, mod3$estimate, sigma = I)
Deltas = array(dim = c(n, N, 2000))

for(t in 1:n){
  print(t)
  for(i in 1:2000){
    Deltas[t,,i] = get_delta_cont2(thetas[i,], t = timeseq[t], K = 2)
  }
}
DeltaCI = array(dim = c(n,N,2))
for(t in 1:n){
  for(j in 1:N){
    DeltaCI[t,j,] = quantile(Deltas[t,j,], probs = c(0.025, 0.975))
  }
}


# Gamma plot

Z = trigBasisExp(timeseq, L = 12, degree = 2)
Gamma_cont = tpm_p(L = 12, beta = beta_hat3, degree = 2, Z = Z)

pdf("./figures/tpm_p.pdf", width = 7, height = 4)
par(mfrow = c(2,2), mar = c(3,4,2,1))
for(i in 1:N){
  for(j in 1:N){
    plot(timeseq, Gamma_cont[i,j,], type = "l", lwd = 2, ylim = c(0,1), bty = "n", ylab = paste0("Pr(",i,"->",j,")"), xaxt = "n", xlab = "")
    # points(Gamma_hat2[i,j,], pch = 19)
    axis(1, at = seq(0,12, by = 2), labels = seq(0,24,by = 4))
  }
}
dev.off()

rho = matrix(nrow = n, ncol = N)
for(i in 1:n){
  rho[i,] = stationary(Gamma_cont[,,i])
}


pdf("./figures/stationary_p.pdf", width = 7, height = 4)

par(mfrow = c(1,1))
plot(timeseq, Delta_cont3[,2], type = "l", bty = "n", lwd = 2, ylim = c(0,1), 
     ylab = "Pr(active)", col = color[2], xlab = "time of day", xaxt = "n")
lines(timeseq, rho[,2], lwd = 2, col = "gold2")
lines(timeseq, Delta_cont3[,2], lwd = 2, col = color[2])
points(Delta_hat3[,2], pch = 19, col = color[2])
polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,2,1], rev(DeltaCI[,2,2])), col = alpha(color[2], 0.2), border = FALSE)
axis(1, at = seq(0,12,by = 2), labels = seq(0,24,by=4))

legend("top", bty = "n", lwd = 2, col = c(color[2], "gold2"), legend = c("stationary distribution", "biased approximation"))

dev.off()



dpmf1 = ddwell(1:24, 1, Gamma_hat3)
dpmf2 = ddwell(1:24, 2, Gamma_hat3)

pdf("./figures/dwell_distr.pdf", width = 8, height = 4)
par(mfrow = c(1,2), mar = c(5,4,4,2)+0.1)
plot(1:12, dpmf1[1:12], type = "h", lwd = 2, bty = "n", col = color[1], 
     xaxt = "n", xlab = "dwell time (hours)", ylab = "probabilities", main = "inactive state",
     ylim = c(0, 0.35), xlim = c(0,12))
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))
plot(1:12, dpmf2[1:12], type = "h", lwd = 2, bty = "n", col = color[2], 
     xaxt = "n", xlab = "dwell time (hours)", ylab = "probabilities", main = "active state",
     ylim = c(0, 0.25), xlim = c(0,12))
axis(1, at = seq(0,12, by = 2), labels = seq(0,24, by = 4))

dev.off()







