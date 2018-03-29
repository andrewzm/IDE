library("STRbook")
library("FRK")
library("sp")
library("DEoptim")
library("dplyr")
library("tidyr")
library("ggplot2")
library("Matrix")
library("gridExtra")
library("spacetime")

rm(list=ls())
source("IDEfunctions.R")


data(STObj3)
STObj4 <- STObj3[, "1993-07-01::1993-07-31"]

IDEmodel <- IDE(f = z ~ lon + lat, 
                data = STObj4, 
                dt = as.difftime(1, units = "days"),
                process_basis = NULL,
                kernel_basis = NULL,
                grid_size = 41)

STgrid <- predict(IDEmodel)

## Optimise log likelihood
optimfun <- function(theta, IDEmodel) {
  ki <- theta[1:4]
  sigma2_eps <- exp(theta[5])
  sigma2_eta <- exp(theta[6])
  ki[1] <- ki[1]*1000
  ki[2] <- exp(ki[2]*10)
  IDEmodel$set(k =  lapply(1:4,function(i) ki[i]), 
               sigma2_eps = sigma2_eps,
               sigma2_eta = sigma2_eta)
  IDEmodel$negloglik()
}

P <- IDEmodel$get("plausible_ranges")

lower = c(0,-1,-10,-10,-5,-10) 
upper = c(1,1,10,10,1,8)

O <- DEoptim(fn = optimfun, 
             lower = c(P$k1[1]/1000,log(P$k2[1])/10, P$k3[1], P$k4[1], 
                       log(P$sigma2_eps[1]), log(P$sigma2_eta[1])), 
             upper = c(P$k1[2]/1000,log(P$k2[2])/10, P$k3[2], P$k4[2],
                       log(P$sigma2_eps[2]), log(P$sigma2_eta[2])), 
             control = list(parallelType = 1, 
                            packages = c("Matrix","FRK", "sp", "dplyr"), 
                            parVar = c("IDEmodel","construct_kernel", 
                                       "rep.col", "construct_Q","Zeromat", 
                                       "logdet")),
             IDEmodel = IDEmodel)

theta <- c(0.710391,   -0.458445,   -0.292635,    0.088533,    0.995981,    1.542510)
ki <- theta[1:4]
sigma2_eps <- exp(theta[5])
sigma2_eta <- exp(theta[6])
ki[1] <- ki[1]*1000
ki[2] <- exp(ki[2]*10)
IDEmodel$set(k =  lapply(1:4,function(i) ki[i]), 
             sigma2_eps = sigma2_eps,
             sigma2_eta = sigma2_eta)

ST_grid_df <- predict(IDEmodel) %>%
  as.data.frame()
gpred <- ggplot(ST_grid_df %>% filter(t %in% c(3,8,13,18,23,28))) + geom_tile(aes(lon, lat,fill=Ypred)) + facet_wrap(~t) +
  scale_fill_distiller(palette="Spectral") 
gpredse <- ggplot(ST_grid_df) + geom_tile(aes(s1,s2,fill=Ypredse)) + facet_wrap(~t) +
  scale_fill_distiller(palette="Spectral") + 
  coord_fixed(xlim=c(0,1), ylim = c(0,1))

# data_df <- as.data.frame(STObj4)
# for(i in 1:31) {
#   g <- ggplot(filter(data_df, day == i & !is.na(z))) + geom_point(aes(lon, lat,colour=z))+
#     scale_colour_distiller(palette="Spectral") 
#   ggsave(g, file=paste0("Test",i,".png"), width = 5, height = 5)
# }

