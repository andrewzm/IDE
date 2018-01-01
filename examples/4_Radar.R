library("STRbook")
library("dplyr")
library("ggplot2")
library("sp")
library("spacetime")
devtools::load_all("..")

rm(list=ls())
#source("IDEfunctions.R")

data("radar_STIDF")
IDEmodel <- IDE(f = z ~ 1,
                data = radar_STIDF,
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

