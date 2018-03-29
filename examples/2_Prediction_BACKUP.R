library("FRK")
library("sp")
library("dplyr")
library("tidyr")
library("ggplot2")
library("Matrix")
library("gridExtra")

rm(list=ls())
source("IDEfunctions.R")

set.seed(1)

spat_points <- SpatialPoints(matrix(runif(100),50,2))
G <- auto_basis(manifold = plane(),
                data = spat_points,
                nres = 1,
                regular = 1)

G_const <- new("Basis", manifold = plane(),
               fn = list(function(s) rep(1, nrow(s))),
               pars = list(),
               df = data.frame(),
               n = 1)

## Spatial decomposition
Y_basis <- auto_basis(manifold = plane(),
                      data = spat_points,
                      regular = 1,
                      nres = 2)
r <- nbasis(Y_basis)

## Simulate Process
T <- 9
beta <- c(0.2,0.2,0.2)
nobs <- 100
sigma2_eta <- 0.01^2
sigma2_eps <- 0.01^2
Sigma_eta <- sigma2_eta * Diagonal(r)
Sigma_eps <- sigma2_eps * Diagonal(nobs * T)
Q_eta <- Sigma_eta %>% solve()
Q_eps <- Sigma_eps %>% solve()

bbox <- matrix(c(0,0,1,1),2,2)
s <- construct_grid(bbox, 41)
Mfun <- construct_M(Y_basis, s)
K_basis <- list(G_const, G_const, G_const, G_const)
M <- Mfun(K_basis, list(150, 0.002, -0.1, 0.1))
alpha <- matrix(0,nbasis(Y_basis),T)
alpha[70,1] <- 1
PHI <- eval_basis(Y_basis, s$s_grid_mat)
s$s_grid_df$Y0 <- (PHI %*% alpha[,1]) %>% as.numeric()
for(i in 1:(T-1)) {
  alpha[,i+1] <- (M %*% alpha[,i]) %>% as.numeric() + sqrt(sigma2_eta)*rnorm(nbasis(Y_basis))
  s$s_grid_df[paste0("Y",i)] <- (PHI %*% alpha[,i+1]) %>% as.numeric()
}
s_long <- gather(s$s_grid_df, time, val, -s1, -s2)
g_truth <- ggplot(s_long) + geom_tile(aes(s1,s2,fill=val)) + facet_wrap(~time) +
  scale_fill_distiller(palette="Spectral", limits = c(-0.1,1.1)) + 
  coord_fixed(xlim=c(0,1), ylim = c(0,1))

## Observe process
zlocs <- data.frame(s1 = runif(nobs),
                    s2 = runif(nobs))
PHI_obs_1 <- eval_basis(Y_basis, zlocs[,1:2] %>% as.matrix())
PHI_obs <- do.call("bdiag", lapply(1:T, function(x) PHI_obs_1))
Xobs <-  cbind(1, do.call("rbind", lapply(1:T, function(x) zlocs))) %>% as.matrix()
Z <- Xobs %*% beta + PHI_obs %*% c(alpha) + sqrt(sigma2_eps) * rnorm(nrow(PHI_obs))
z_df <- data.frame(expand.grid.df(zlocs, data.frame(t = 0:(T-1))))
z_df$z <- Z %>% as.numeric()
g_obs <- ggplot(z_df) + geom_point(aes(s1, s2, colour = z)) + 
  facet_wrap(~t) + 
  scale_colour_distiller(palette = "Spectral") + 
  coord_fixed(xlim=c(0,1), ylim = c(0,1))

## Prediction -- one time point
M <- Mfun(K_basis, list(150, 0.002, -0.1, 0.1))
Q <- construct_Q(Q_eta, M, T)
Qpost <- crossprod(chol(Q_eps) %*% PHI_obs) + Q
Qpostchol <- FRK:::cholPermute(Qpost)
mupost <- FRK:::cholsolve(Q = Qpost, y = t(PHI_obs) %*% Q_eps %*% (Z - Xobs %*% beta), perm = TRUE, 
                          cholQp = Qpostchol$Qpermchol, P = Qpostchol$P)
PHI_pred <- do.call("bdiag", lapply(1:T, function(x) PHI))

ST_grid_df <- s$s_grid_df %>% 
  select(s1,s2) %>%
  expand.grid.df(data.frame(t = 0:(T-1)))
ST_grid_df$Ypred <- (PHI_pred %*% mupost) %>% as.numeric()            

gpred <- ggplot(ST_grid_df) + geom_tile(aes(s1,s2,fill=Ypred)) + facet_wrap(~t) +
  scale_fill_distiller(palette="Spectral", limits = c(-0.1,1.1)) + 
  coord_fixed(xlim=c(0,1), ylim = c(0,1))
arrangeGrob(g_obs, g_truth, gpred, ncol = 3) %>% plot()

calc_betahat <- function(X_obs, PHI_obs, Q_eps, Qpost_cholsolve, Z) {
  tPHIQepsZ <- t(PHI_obs) %*% Q_eps %*% Z
  tXQeps <- (t(Xobs) %*% Q_eps)
  tXQepsPHI <- tXQeps %*% PHI_obs
  Part1 <- solve(tXQeps %*% Xobs - tXQepsPHI %*% Qpost_cholsolve(t(tXQepsPHI)))
  Part2 <- tXQeps %*% Z - tXQepsPHI %*% Qpost_cholsolve(y = tPHIQepsZ)
  betahat <- Part1 %*% Part2
}

## Log likelihood
loglik <- function(theta) {
  ki <- theta[1:4]
  sigma2_eps <- exp(theta[5])
  sigma2_eta <- exp(theta[6])
  
  Q_eps <- 1/sigma2_eps * Diagonal(length(Z))
  Q_eta <- 1/sigma2_eta * Diagonal(r)
  
  ki[1] <- ki[1]*1000
  ki[2] <- exp(ki[2]*10)
  M <- Mfun(K_basis, list(ki[1], ki[2], ki[3], ki[4]))
  
  if(max(abs(eigen(M)$values)) < 1 & abs(ki[3]) < 0.5 & abs(ki[4]) < 0.5) {
    Q <- construct_Q(Q_eta, M, T)
    Qpost <- crossprod(chol(Q_eps) %*% PHI_obs) + Q
    Qpostchol <- FRK:::cholPermute(Qpost)
    
    Qpost_cholsolve <- function(y) {
      FRK:::cholsolve(Q = Qpost, 
                      y = y, 
                      perm = TRUE, 
                      cholQp = Qpostchol$Qpermchol, 
                      P = Qpostchol$P)
    }
    
    betahat <- calc_betahat(X_obs, PHI_obs, Q_eps, Qpost_cholsolve, Z)
    Ztilde <- Z - Xobs %*% betahat
    
    mupost <- Qpost_cholsolve(y = t(PHI_obs) %*% Q_eps %*% Ztilde)
    Qchol <- FRK:::cholPermute(Q)
    -((0.5*logdet(Qchol$Qpermchol) + 
         0.5*logdet(chol(Q_eps)) - 
         0.5*logdet(Qpostchol$Qpermchol) -
         0.5*t(Ztilde) %*% (Q_eps %*% Ztilde) +
         0.5*t(Ztilde) %*% Q_eps %*% PHI_obs %*% mupost) %>% as.numeric())
  } else {
    1e10
  }
}

dfsfs
#O <- optim(c(150, log(0.002), -0.1, 0.1), fn = loglik)
#O <- nlm(f = loglik, p = c(100/1000, log(0.004)/10, 0.1, -0.1))
O <- DEoptim(fn = loglik, 
             lower = c(0,-1,-1,-1,-10,-10), 
             upper = c(1,1,1,1,10,10),
             control = list(parallelType = 1, 
                            packages = c("Matrix","FRK", "sp", "dplyr"), 
                            parVar = c("Z","r","Mfun","K_basis",
                                       "PHI_obs","construct_kernel", 
                                       "rep.col", "construct_Q","Zeromat", 
                                       "T", "logdet", "Xobs", "calc_betahat")))