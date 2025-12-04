library(readr)
library(shiny)
library(MASS)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(patchwork)
library(viridis)
library(readxl)
library(gridExtra)
library(tidyr)
library(RColorBrewer)
library(stats)
library(WebPower)
library(matrixStats)
library(glmm)
library(gam)
library(MASS)
library(betareg)
library(stats)
library(gamlss)
library(gamlss.dist)
library(rjags)
library(zoib)
library(lme4)
library(glmmTMB)
library(INLA)
library(inlabru)
library(ggfortify)
library(mgcv)
library(binr)
library(astsa)
library(viridis)
library(ggspatial) # trying to get scale
library(osmdata) # one way to get pretty street map
library(corrplot)
library(readr)
library(shiny)
library(MASS)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(patchwork)
library(viridis)
library(readxl)
library(geobr)
library(sf)
#library(rgdal)
#library(maptools)
library(gridExtra)
library(tidyr)
library(RColorBrewer)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, sf, terra, exactextractr, raster, gstat)
library(simr)
library(lme4)
library(writexl)
library(pwr)
library(CRTSize)




## Get Average Population in Clusters

## population data 
wp_dens_rast2020 <- rast("bra_ppp_2020_UNadj_constrained.tif")

## BH Clusters
regions = st_read("Regional.shp")
aea.proj <- "+proj=longlat +datum=WGS84 +no_defs"
regions = st_transform(regions,aea.proj)

finalschoolsp <- st_read("Clusters_20210107_Final.kml",layer='Clusters_20210107_Final')

new <- st_transform(finalschoolsp, aea.proj)
finalschoolsp_sf <- st_as_sf(finalschoolsp)
finalschoolsp_sf <- st_transform(finalschoolsp_sf, aea.proj)

# ggplot(finalschoolsp_sf)+
#   geom_sf()

finalschoolsp_sf <- finalschoolsp_sf %>%
  dplyr::select(!(Description))

## combine
population <- exact_extract(wp_dens_rast2020, finalschoolsp_sf, fun = "sum", append_cols = "Name") |>
  rename(pop_wp = sum) 

avg_pop <- population %>%
  summarize(median = median(pop_wp), # 10624.33	
            mean = mean(pop_wp)) # 11528.18


write_xlsx(population, "cluster_populations_wp.xlsx")
