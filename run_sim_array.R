# libraries #################################################################################################
library(here)
library(rio)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(survival)

# command line arguments ####################################################################################
args=(commandArgs(TRUE))
if (length(args) < 3) { #check to see if arguments are passed. 
  stop("Usage: Rscript run_sim.R <block_R0> <cluster_R0> <gamma>")
}

block_R0   <- as.numeric(args[1])
cluster_R0 <- as.numeric(args[2])
gamma      <- as.numeric(args[3])

# load functions ############################################################################################
source("code/src/sim_functions.R")
source("code/src/sim_constants.R")

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
# # only want cases that are during the ABD season -- week 26 (July 1 2023)
# filter(Week_onset >= 26)

# incidence in merida??
city_cases_by_week <- sinave_dengue_cases %>% 
  count(Week_onset, name = "cases_per_week")

## divided by population...
yuc_pop <- read.csv(here("data", "SINAVE_Data", "yuc_data.csv"))
I_m <- city_cases_by_week$cases_per_week #/ yuc_pop$POBTOT[1]

# make population and generation time ######################################################################
gen_prob <- gen_time(Gmu, Gs2, max_gen_time)
trial_pop <- sim_pop(n_clusters, n_per_cluster, n_per_block)

## run simulation ###########################################################################################
sim_out <- run_sims(
  num_sims = 5, 
  trial_pop = trial_pop, 
  gen_prob = gen_prob, 
  initially_infected = 5, 
  alpha = alpha, 
  beta = beta, 
  gamma = gamma, 
  r = 0.5
)

## print sim results #######################################################################################
snum <- function(x) gsub("\\.", "p", sprintf("%.2f", x)) # p for decimal point
file_stem <- paste0(
  "block", snum(block_R0),
  "_cluster", snum(cluster_R0),
  "_gamma", snum(gamma)
)
outfile <- file.path("output/sim_results", paste0(file_stem, ".rds"))

saveRDS(sim_out, file = outfile)
cat("Finished -- results written to ", outfile, "\n")