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

## Use census tract data to estimate census variables in study clusters

# Population density and income from census tracts
census_tracts = st_read('.','SETOR_CENSITARIO_2010')


options(scipen=999)
pop_by_ct = read_xlsx('./MG_CensusData.xlsx')

# Restrict to those in Belo Horizonte
census_tracts$COD_GEOCOD=as.numeric(as.character(census_tracts$COD_GEOCOD))
pop_by_ct = pop_by_ct[pop_by_ct$Cod_setor %in% census_tracts$COD_GEOCOD,]

# Now get population and income information into shapefile
census_tracts = census_tracts[census_tracts$COD_GEOCOD %in% pop_by_ct$Cod_setor,]
census_tracts$Pop_IBGE = sapply(census_tracts$COD_GEOCOD,function(x) as.numeric(pop_by_ct$Pop_Total[pop_by_ct$Cod_setor==x]))
census_tracts$MeanHHIncome_10plus_WithIncome = sapply(census_tracts$COD_GEOCOD,function(x) as.numeric(pop_by_ct$MeanHHIncome_10plus_WithIncome[pop_by_ct$Cod_setor==x]))


census_tracts$CTArea = st_area(census_tracts)/1000000
census_tracts$PopDens = census_tracts$Pop_IBGE/census_tracts$CTArea
finalschoolsp_sf = st_transform(finalschoolsp_sf,st_crs(census_tracts))

finalschoolsp_sf$WeightedIncome = sapply(1:nrow(finalschoolsp_sf),function(x) {
  print(x)
  ct_int = st_intersection(census_tracts,finalschoolsp_sf[x,])
  ct_int$newarea = st_area(ct_int)/1000000
  ct_int$newinc = ct_int$MeanHHIncome_10plus_WithIncome*ct_int$newarea/ct_int$Ar_sqkm
  return(sum(ct_int$newinc))
})





