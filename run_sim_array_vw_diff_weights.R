# libraries #################################################################################################
library(here)
library(rio)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(survival)

# command line arguments ####################################################################################
args <- commandArgs(TRUE)
if (length(args) < 3) stop("Usage: Rscript run_sim_array.R <comm_weight_kids> <comm_weight_adults> <run_chunk>")

comm_weight_kids <- as.numeric(args[1])
comm_weight_adults <- as.numeric(args[2])
chunk_id <- as.integer(args[3])           # 1 … 10

# limit the number of sims per chunk
sims_per_chunk <- 5

# load functions ############################################################################################
source("code/src/sim_functions_vw_diff_weights.R")
source("code/src/sim_constants_vary_weight.R")

# calculate parameters ######################################################################################
alpha <- -log((block_R0 + 1 - n_per_block) / (1 - n_per_block))
beta  <- -log((cluster_R0 - n_per_cluster + n_per_block) /
                (n_per_block - n_per_cluster))

# import data ###############################################################################################
## SINAVE case data, only dengue
sinave_dengue_cases <- here::here("data", "SINAVE_Data", 
                                  "dengue_abierto.csv") %>%
  rio::import() %>% 
  dplyr::as_tibble() %>% 
  # only want data for Mérida + Kanasín (050, 041 MUNICIPO_RES)
  filter(MUNICIPIO_RES == 50 | MUNICIPIO_RES == 41) %>% 
  # change symptom onset date to match cases dataset (i.e., weeks since 
  #  start of year)
  mutate(Week_onset = ymd(FECHA_SIGN_SINTOMAS),
         Week_onset = epiweek(Week_onset),
         Year = year(FECHA_SIGN_SINTOMAS)) %>%
  relocate(c(FECHA_SIGN_SINTOMAS, Week_onset), .before = ID_REGISTRO) %>% 
  filter(ESTATUS_CASO == 2) %>%  # 1 likely; 2 confirmed; 3 discarded
  mutate(age_grp = as.factor(ifelse(EDAD_ANOS < 18, "Child", "Adult"))) #%>% 

# incidence in merida
city_cases_by_week <- sinave_dengue_cases %>% 
  count(Week_onset, name = "cases_per_week")

## divided by population...
merida_pop <- read.csv("data/yuc_pop_all.csv") %>% # all of yucatan
  filter(MUN == 50 | MUN == 41 | MUN == 101) %>% # only metro area of merida
  filter(NOM_LOC == "Total del municipio") # only total pops

I_m <- city_cases_by_week$cases_per_week #/ sum(merida_pop$POBTOT)

# make population and generation time ######################################################################
set.seed(12345)
gen_prob <- gen_time(Gmu, Gs2, max_gen_time)
trial_pop <- sim_pop(n_clusters, n_per_cluster, n_per_block)

## run simulation ###########################################################################################
offset_per_chunk <- (chunk_id - 1) * sims_per_chunk # chunk1 is runs 1-10, chunk2 is 11-20

sim_out <- run_sims(
  num_sims = sims_per_chunk, 
  trial_pop = trial_pop, 
  gen_prob = gen_prob, 
  initially_infected = 25, 
  alpha = alpha, 
  beta = beta, 
  gamma = gamma, 
  r = 0.039,
  comm_weight_kids = comm_weight_kids,
  comm_weight_adults = comm_weight_adults,
  seed_offset = offset_per_chunk,
  chunk_id = chunk_id,
  sims_per_chunk = sims_per_chunk
)