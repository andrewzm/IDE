library("FRK")
library("sp")
library("dplyr")
library("tidyr")
library("ggplot2")
library("Matrix")

rm(list=ls())
source("IDEfunctions.R")

set.seed(1)

spat_points <- SpatialPoints(matrix(runif(100),50,2))
Y_basis <- auto_basis(manifold = plane(),
                data = spat_points,
                nres = 1,
                regular = 1)
G_const <- new("Basis", manifold = plane(),
               fn = list(function(s) rep(1, nrow(s))),
               pars = list(),
               df = data.frame(),
               n = 1)

K_basis <- list(G_const, G_const, G_const, G_const)
K <- construct_kernel(K_basis, list(0.1, 0.01, -0.1, 0.1))

bbox <- matrix(c(-1,1,-1,1),2,2)
s <- construct_grid(bbox, 41)
s$s_grid_df <- s$s_grid_df %>% 
               mutate(Y0 = exp(-sqrt(s1^2 + s2^2)/0.2))

Kmat <- K(s$s_grid_mat, s$s_grid_mat)
s$s_grid_df <- s$s_grid_df %>% 
               mutate(Y1 = (Kmat %*% Y0) %>% as.numeric(),
                      Y2 = Kmat %*% Y1 %>% as.numeric())

ggplot(s$s_grid_df) + geom_tile(aes(s1,s2,fill=Y0))
ggplot(s$s_grid_df) + geom_tile(aes(s1,s2,fill=Y1))