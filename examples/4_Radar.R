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
devtools::load_all("..")

rm(list=ls())
#source("IDEfunctions.R")

s_grid <- expand.grid(s2 = seq(98.75,length.out = 40, by = -2.5),
                      s1 = seq(1.25,length.out = 28, by = 2.5))
s_df <- NULL

time <- seq(as.POSIXct("2000-11-03 08:25:00", tz = "UTC"),
             length.out = 12,
             by = as.difftime(10, units = "mins"))
for(i in 1:12) {
    X <- read.table(paste0("../../data/radar/Z",i,".dat"))
    new_df <- cbind(s_grid, t = time[i], z = as.numeric(as.matrix(X)))
    s_df <- rbind(s_df, new_df)
}

data_STIDF <- STIDF(sp = s_df[,1:2] %>% SpatialPoints(),
                    time = s_df$t,
                    data = s_df %>% select(z))

IDEmodel <- IDE(f = z ~ 1,
                data = data_STIDF,
                dt = as.difftime(10, units = "mins"),
                process_basis = NULL,
                kernel_basis = NULL,
                grid_size = 41,
                forecast = 2)

#fit_results <- fit.IDE(IDEmodel)
load("4_Radar.rda")

ST_grid_df <- predict(fit_results$IDEmodel) %>%
  as.data.frame()

s_df$time <- format(s_df$t, "%H:%M")
ggplot(s_df) + geom_tile(aes(s1,s2,fill=pmin(pmax(z,-20),60))) + facet_wrap(~time) +
    scale_fill_distiller(palette = "Spectral", limits = c(-20,60)) + coord_fixed()
gpred <- ggplot(ST_grid_df) + geom_tile(aes(s1, s2,fill=Ypred)) + facet_wrap(~t) +
  scale_fill_distiller(palette="Spectral", limits = c(-20,60))
gpredse <- ggplot(ST_grid_df) + geom_tile(aes(s1,s2,fill=Ypredse)) + facet_wrap(~t) +
  scale_fill_distiller(palette="Spectral")

# data_df <- as.data.frame(STObj4)
# for(i in 1:31) {
#   g <- ggplot(filter(data_df, day == i & !is.na(z))) + geom_point(aes(lon, lat,colour=z))+
#     scale_colour_distiller(palette="Spectral")
#   ggsave(g, file=paste0("Test",i,".png"), width = 5, height = 5)
# }

