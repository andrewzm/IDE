library("FRK")
library("sp")
library("DEoptim")
library("dplyr")
library("tidyr")
library("ggplot2")
library("Matrix")
library("gridExtra")
library("spacetime")
devtools::load_all("..")

rm(list=ls())
#source("IDEfunctions.R")

## Simulate data from IDE
SIM <- simIDE(T = 10, nobs = 100)

IDEmodel <- IDE(f = z ~ s1 + s2,
                data = SIM$z_STIDF,
                dt = as.difftime(1, units = "days"),
                process_basis = NULL,
                kernel_basis = NULL,
                grid_size = 41)

#fit_results <- fit.IDE(IDEmodel)
load("./2_Prediction.rda")
show_kernel(fit_results$IDEmodel)

## Prediction grid
s1_pred <- s2_pred <- seq(0,1,length.out = 71)
st_grid <- expand.grid(s1 = s1_pred,
                      s2 = s2_pred,
                      date = unique(time(SIM$z_STIDF))) %>%
           mutate(t = lubridate::day(date))
pred_grid <- STIDF(sp = SpatialPoints(st_grid[,c("s1","s2")]),
                 time = st_grid$date,
                 data = st_grid %>% select(-s1, -s2, -date))

## Predict using prior guesses
ST_grid_df <- predict(IDEmodel,
                      newdata = pred_grid) %>%
              as.data.frame()

gpred <- ggplot(ST_grid_df) + geom_tile(aes(s1,s2,fill=Ypred)) + facet_wrap(~t) +
  scale_fill_distiller(palette="Spectral", limits = c(-0.1,1.1)) +
  coord_fixed(xlim=c(0,1), ylim = c(0,1))
gpredse <- ggplot(ST_grid_df) + geom_tile(aes(s1,s2,fill=Ypredse)) + facet_wrap(~t) +
  scale_fill_distiller(palette="Spectral") +
  coord_fixed(xlim=c(0,1), ylim = c(0,1))
# arrangeGrob(SIM$g_obs, SIM$g_truth, gpredse, gpred, ncol = 2) %>% plot()
#
# ## Optimise log likelihood
# optimfun <- function(theta, IDEmodel) {
#     ki <- theta[1:4]
#     sigma2_eps <- exp(theta[5])
#     sigma2_eta <- exp(theta[6])
#     ki[1] <- ki[1]*1000
#     ki[2] <- exp(ki[2]*10)
#     IDEmodel$set(k =  lapply(1:4,function(i) ki[i]),
#                  sigma2_eps = sigma2_eps,
#                  sigma2_eta = sigma2_eta)
#     IDEmodel$negloglik()
# }
#
# #lower <- c(0,-1,-1,-1,-10,-10)
# #upper <- c(1,1,1,1,10,10)
# P <- IDEmodel$get("plausible_ranges")
# lower = c(P$k1[1]/1000,log(P$k2[1])/10, P$k3[1], P$k4[1],
#           log(P$sigma2_eps[1]), log(P$sigma2_eta[1]))
# upper = c(P$k1[2]/1000,log(P$k2[2])/10, P$k3[2], P$k4[2],
#           log(P$sigma2_eps[2]), log(P$sigma2_eta[2]))
# print("Changed lower and upper")
#
# O <- DEoptim(fn = optimfun,
#              lower = lower,
#              upper = upper,
#              control = list(parallelType = 1,
#                             packages = c("Matrix","FRK", "sp", "dplyr"),
#                             parVar = c("IDEmodel","construct_kernel",
#                                        "rep.col", "construct_Q","Zeromat",
#                                        "logdet")),
#              IDEmodel = IDEmodel)
#
# ## Predict using optimised parameters
# theta <- c(0.071683,   -0.540418,   -0.104860,    0.099691,   -9.026604,   -6.214601)
# ki <- theta[1:4]
# sigma2_eps <- exp(theta[5])
# sigma2_eta <- exp(theta[6])
# ki[1] <- ki[1]*1000
# ki[2] <- exp(ki[2]*10)
# IDEmodel$set(k =  lapply(1:4,function(i) ki[i]),
#              sigma2_eps = sigma2_eps,
#              sigma2_eta = sigma2_eta)
#
# ST_grid_df <- predict(IDEmodel,
#                       newdata = pred_grid) %>%
#   as.data.frame()
# gpred <- ggplot(ST_grid_df) + geom_tile(aes(s1,s2,fill=Ypred)) + facet_wrap(~t) +
#   scale_fill_distiller(palette="Spectral", limits = c(-0.1,1.1)) +
#   coord_fixed(xlim=c(0,1), ylim = c(0,1))
# gpredse <- ggplot(ST_grid_df) + geom_tile(aes(s1,s2,fill=Ypredse)) + facet_wrap(~t) +
#   scale_fill_distiller(palette="Spectral") +
#   coord_fixed(xlim=c(0,1), ylim = c(0,1))
# arrangeGrob(SIM$g_obs, SIM$g_truth, gpredse, gpred, ncol = 2) %>% plot()
#
#
