library(here)
library(rio)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(survival)

# set seed
set.seed(12345)

# import incidence data ############################################################################################################
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

# functions ###############################################################
## trial population structure ################################################
sim_pop <- function(n_clusters = 50, n_per_cluster, n_per_block, 
                    min_age = 1, max_age = 82){
  
  n <- n_clusters * n_per_cluster 
  
  trial_pop_df <- data.frame(
    "ID" = 1:n,
    "Cluster" = rep(1:n_clusters, each = n_per_cluster),
    "Block" = rep_len(1:round(n_per_cluster/n_per_block, 0), n_clusters),
    "TIRS" = rep(c(0, 1), each = n_per_cluster, length.out = n),
    "Age" = sample(min_age:max_age, n, replace = TRUE)
  )
  
  trial_pop_df$cluster_block <- paste(trial_pop_df$Cluster, 
                                      trial_pop_df$Block, sep = "")
  return(trial_pop_df)
}

## infect population ##########################################################################################################
### generation time ###########################################################################################################
# vector of generation time pmf given the parameters; 
# evaluated at each time point
gen_time <- function(Gmu, Gs2, max_gen_time) {
  shape <- Gmu^2 / Gs2
  scale <- Gs2 / Gmu
  
  i <- 0:max_gen_time
  
  genProb <- pgamma(i + 1, shape = shape, scale = scale) -
    pgamma(i,     shape = shape, scale = scale)
  
  return(genProb)
}

### infections ####################################################################################################################
infect_pop <- function(trial_pop, gen_prob, initially_infected = 1, 
                       alpha, beta, gamma, r) {
  
  N <- nrow(trial_pop)
  trial_pop$infected <- 0
  trial_pop$tau <- NA
  trial_pop$infector <- NA
  trial_pop$source <- NA
  
  # susceptibility
  trial_pop$susceptible <- 1
  adult_ids <- which(trial_pop$Age >= 18)
  immune_adults <- sample(adult_ids, size = round(0.05 * length(adult_ids))) # 5% adults immune??
  trial_pop$susceptible[immune_adults] <- 0
  
  # tracking infectors
  susceptible_ids <- which(trial_pop$susceptible == 1)
  initial_ids <- sample(susceptible_ids, initially_infected)
  trial_pop$infected[initial_ids] <- 1
  trial_pop$tau[initial_ids] <- trial_start
  trial_pop$index_case <- 0
  trial_pop$index_case[initial_ids] <- 1
  
  for (t  in (trial_start + 1):trial_end){
    
    new_infections <- c()
    
    for (j in 1:N){
      if (trial_pop$infected[j] == 1) { next }
      if (trial_pop$susceptible[j] == 0) { next }
      
      lambda_jb <- 0
      lambda_jc <- 0
      lambda_jm <- 0
      
      # block FOI
      same_block <- which(trial_pop$cluster_block ==
                            trial_pop$cluster_block[j] & (1:N) != j)
      for (i in same_block){
        if (trial_pop$infected[i] == 1) {
          dt <- t - trial_pop$tau[i]
          if (dt <= max_gen_time){ 
            if (trial_pop$TIRS[j] == 1) {
              lambda_jb <- lambda_jb + r * alpha * gen_prob[dt]
              # treatment reduction * hazard from i infecting j * prob that gen time is dt
            } else {
              lambda_jb <- lambda_jb + alpha * gen_prob[dt]
              # hazard from i infecting j * prob that gen time is dt
            }
          }
        }
      }
      
      # cluster FOI
      same_cluster <- which(trial_pop$Cluster == trial_pop$Cluster[j] 
                            & (1:N) != j)
      diff_block <- setdiff(same_cluster, same_block)
      for (i in diff_block) {
        if (trial_pop$infected[i] == 1) {
          dt <- t - trial_pop$tau[i]
          if (dt <= max_gen_time){ 
            if (trial_pop$TIRS[j] == 1) {
              lambda_jc <- lambda_jc + r * beta * gen_prob[dt]
            } else {
              lambda_jc <- lambda_jc + beta * gen_prob[dt]
              # cum. hazard from i infecting j * prob that gen time is dt
            }
          }
        }
      }
      
      # merida FOI
      for (z in trial_start:(t-1)) {
        dt <- t - z
        if (dt <= max_gen_time) {
          lambda_jm <- lambda_jm + gamma * (I_m[z] / yuc_pop$POBTOT[1]) * gen_prob[dt]
          # hazard from a case in merid * cases rate in merida * prob gen time is dt
        }
      }
      
      # total FOI on j
      lambda_j <- lambda_jb + lambda_jc + lambda_jm
      
      # probability of infection given hazard
      p_j <- 1 - exp(-lambda_j)
      
      # j becomes a new infection
      if (runif(1) < p_j) {
        # possible infectors of j
        sources <- c(same_block, diff_block, 0)
        weights <- c(sapply(same_block, function(i) {
          if (trial_pop$infected[i] == 1) {
            dt <- t - trial_pop$tau[i]
            if (dt <= max_gen_time){ 
              base <- alpha * gen_prob[dt]
              # hazard from i infecting j * prob that gen time is dt
              if (trial_pop$TIRS[i] == 1) { base <- r * base }
              return(base)
            }
          }
          return(0)
        }), sapply(diff_block, function(i){
          if (trial_pop$infected[i] == 1) {
            dt <- t - trial_pop$tau[i]
            if (dt <= max_gen_time){ 
              base <- beta * gen_prob[dt]
              if (trial_pop$TIRS[i] == 1) { base <- base * r }
              return(base)
            }
          }
          return(0)
        }), lambda_jm)
        
        # pick infector
        if (sum(weights) > 0) {
          infector <- sample(sources, 1, prob = weights)
        } else {
          infector <- 0 # infection is imported
        }
        
        new_infections <- c(new_infections, j)
        if (infector == 0) {
          trial_pop$infector[j] <- 0 # imported infection
          trial_pop$source[j] <- 0
        } else {
          trial_pop$infector[j] <- trial_pop$ID[infector]
          if (infector %in% same_block){
            trial_pop$source[j] <- 1
          } else {
            trial_pop$source[j] <- 2
          }
        }
      }
    }
    
    if (length(new_infections) > 0) {
      trial_pop$infected[new_infections] <- 1
      trial_pop$tau[new_infections] <- t
    }
  }
  
  return(trial_pop = trial_pop)
}

## run simulations ##########################################################################################################
run_sims <- function(num_sims, trial_pop, gen_prob, initially_infected, alpha, 
                     beta, gamma, r) {
  
  sim_list <- list()
  offspring_index <- c()
  rr_peak <- c()
  time_diff_at_peaks <- c()
  control_list <- list()
  tirs_list <- list()
  prop_infected <- c()
  control_infected <- c()
  tirs_infected <- c()
  pois_efficacy <- c()
  pois_popoff_eff <- c()
  pois_ptoff_eff <- c()
  cox_efficacy <- c()
  source_community <- c()
  source_block <- c()
  source_cluster <- c()
  
  real_df <- data.frame(
    time = 1:trial_end,
    case_rate = city_cases_by_week$cases_per_week / yuc_pop$POBTOT[1])
  
  for (s in 1:num_sims) {
    infected_pop <- infect_pop(trial_pop, gen_prob, 1, alpha = alpha,
                               beta = beta, gamma = gamma, r = r)
    
    # epidemic curve
    ## average incidence across simulations
    incidence <- infected_pop %>% 
      filter(is.finite(tau)) %>% 
      count(tau, name = "cases") %>% 
      complete(tau = trial_start:trial_end, fill = list(cases = 0)) %>%
      rename(time = tau) %>% 
      mutate(case_rate = cases / (n_clusters*n_per_cluster),
             sim = s)
    
    sim_list[[s]] <- incidence
    
    ## by trial arm
    control <- infected_pop %>% 
      filter(is.finite(tau)) %>% 
      filter(TIRS == 0) %>% 
      count(tau, name = "cases") %>% 
      complete(tau = trial_start:trial_end, fill = list(cases = 0)) %>%
      rename(time = tau) %>% 
      mutate(sim = s)
    control_list[[s]] <- control
    
    tirs <- infected_pop %>% 
      filter(is.finite(tau)) %>% 
      filter(TIRS == 1) %>% 
      count(tau, name = "cases") %>% 
      complete(tau = trial_start:trial_end, fill = list(cases = 0)) %>%
      rename(time = tau) %>% 
      mutate(sim = s)
    tirs_list[[s]] <- tirs
    
    # count offspring
    index_ids <- infected_pop$ID[infected_pop$index_case == 1]
    count_offspring <- sapply(index_ids, function(id) {
      sum(infected_pop$infector == id, na.rm = TRUE)
    })
    offspring_index[s] <- mean(count_offspring)
    
    # find peak incidence in simulation
    study_peak <- max(incidence$case_rate)
    real_peak <- max(real_df$case_rate)
    
    rr_peak[s] <- study_peak / real_peak # difference in magnitude
    time_diff_at_peaks[s] <- incidence$time[which.max(incidence$case_rate)] - 
      real_df$time[which.max(real_df$case_rate)]
    
    # proportion of population infected at end of season
    prop_infected[s] <- sum(incidence$cases) / (n_clusters*n_per_cluster)
    
    # num infected by trial arm
    control_infected[s] <- sum(control$cases)
    tirs_infected[s] <- sum(tirs$cases)
    
    # Poisson regression (no offset)
    cases_by_cluster <- infected_pop %>% 
      group_by(Cluster, TIRS) %>% 
      summarize(num_cases = sum(infected), .groups = "drop")
    
    pois_fit <- glm(num_cases ~ TIRS, family = poisson(link = "log"), 
                    data = cases_by_cluster)
    pois_efficacy[s] <- 1 - exp(coef(pois_fit)["TIRS"])
    
    # Poisson regression (offset by population)
    cases_pop_by_cluster <- cbind(cases_by_cluster, 
                                  n_per_cluster = rep(n_per_cluster, 
                                                      nrow(cases_by_cluster)))
    
    pois_ofit_pop <- glm(num_cases ~ TIRS + offset(log(n_per_cluster)), 
                         family = poisson(link = "log"), 
                         data = cases_pop_by_cluster)
    pois_popoff_eff[s] <- 1 - exp(coef(pois_ofit_pop)["TIRS"])
    
    # Poisson regression (offset by person-time?)
    infected_pop$person_time <- ifelse(infected_pop$infected == 1, 
                                       infected_pop$tau, trial_end)
    cases_per_cluster_pt <- infected_pop %>% 
      group_by(Cluster, TIRS) %>% 
      summarize(num_cases = sum(infected),
                person_time = sum(person_time),
                .groups = "drop")
    
    pois_ofit_pt <- glm(num_cases ~ TIRS + offset(log(person_time)), 
                        family = poisson(link = "log"), 
                        data = cases_per_cluster_pt)
    pois_ptoff_eff[s] <- 1 - exp(coef(pois_ofit_pt)["TIRS"])
    
    # Cox model (as done in paper)
    # first remove index case
    cox_dat <- infected_pop[infected_pop$index_case != 1,]
    # itt.cox <- coxph(Surv(entry, exit, Pos_All) ~ arm + cluster(cluster),
    #                   data = data)
    cox_fit <- coxph(Surv(rep(26, nrow(cox_dat)), 
                          tau, infected) ~ TIRS + cluster(Cluster),
                     data = cox_dat)
    cox_efficacy[s] <- 1 - exp(coef(cox_fit))
    
    # infector source: num infected per source / total pop infected
    source_community[s] <- sum(infected_pop$source == 0, na.rm = TRUE) / 
      sum(infected_pop$infected)
    source_block[s] <- sum(infected_pop$source == 1, na.rm = TRUE) / 
      sum(infected_pop$infected)
    source_cluster[s] <- sum(infected_pop$source == 2, na.rm = TRUE) / 
      sum(infected_pop$infected)
  }
  
  # R_0
  avg_offspring <- mean(offspring_index)
  var_offspring <- var(offspring_index)
  
  # rate ratio between study population and merida
  avg_rr_peak <- mean(rr_peak, na.rm = TRUE)
  var_rr_peak <- var(rr_peak, na.rm = TRUE)
  avg_time_diff <- mean(time_diff_at_peaks, na.rm = TRUE)
  
  # proportion of population infected
  avg_prop_infected <- mean(prop_infected, na.rm = TRUE)
  
  # avg infected by trial arm
  avg_control_infected <- mean(control_infected, na.rm = TRUE)
  avg_tirs_infected <- mean(tirs_infected, na.rm = TRUE)
  
  # avg efficacy estimates
  avg_pois_estimate <- mean(pois_efficacy, na.rm = TRUE)
  avg_pois_pop_estimate <- mean(pois_popoff_eff, na.rm = TRUE)
  avg_pois_pt_estimate <- mean(pois_ptoff_eff, na.rm = TRUE)
  avg_cox_estimate <- mean(cox_efficacy, na.rm = TRUE)
  
  # avg proportion of infector sources
  avg_infected_within_block <- mean(source_block, na.rm = TRUE)
  avg_infected_within_cluster <- mean(source_cluster, na.rm = TRUE)
  avg_infected_from_community <- mean(source_community, na.rm = TRUE)
  
  return(list(avg_offspring = avg_offspring, var_offspring = var_offspring,
              avg_rr_peak = avg_rr_peak, var_rr_peak = var_rr_peak,
              avg_time_diff = avg_time_diff, 
              avg_prop_infected = avg_prop_infected, 
              avg_control_infected = avg_control_infected,
              avg_tirs_infected = avg_tirs_infected,
              avg_pois_estimate = avg_pois_estimate,
              avg_pois_pop_estimate = avg_pois_pop_estimate,
              avg_pois_pt_estimate = avg_pois_pt_estimate,
              avg_cox_estimate = avg_cox_estimate,
              avg_infected_within_block = avg_infected_within_block,
              avg_infected_within_cluster = avg_infected_within_cluster,
              avg_infected_from_community = avg_infected_from_community))
}

# run simulations ##########################################################################################################
## trial population #############################################################
min_age <- 1 # youngest age
max_age <- 82 # oldest age

trial_start <- 26
trial_end <- 50 # 52 is the last week of the year but we don't have city incidence data for that??? so 50

# parameters
Gmu <- 2; Gs2 <- 1; max_gen_time <- 90 # 17 days?; 2 days?
# probability that the generation time = [index]
gen_prob <- gen_time(Gmu, Gs2, max_gen_time)

# Population size
n_clusters <- 50; n_per_cluster <- 2500; n_per_block <- 250
trial_pop <- sim_pop(n_clusters, n_per_cluster, n_per_block)

## parameters to test #############################################################
param_grid <- expand.grid(block_R0   = c(1.6, 1.9, 2.2, 2.5), 
                          cluster_R0 = c(1.6, 1.9, 2.2, 2.5), 
                          gamma = c(0.5, 1, 2)) %>%
  filter(block_R0 >= cluster_R0) %>%
  mutate(alpha = -log((block_R0 + 1 - n_per_block) / (1 - n_per_block)), 
         beta  = -log((cluster_R0 - n_per_cluster + n_per_block) /
                        (n_per_block - n_per_cluster)))

grid_results <- vector("list", nrow(param_grid))

## run grid search ###########################################################
for (p in 1:nrow(param_grid)) {
  a <- param_grid$alpha[p]
  b <- param_grid$beta[p]
  g <- param_grid$gamma[p]
  
  sim_out <- run_sims(25, trial_pop, gen_prob, 1, a, b, g, r = 0.5)
  
  grid_results[[p]] <- data.frame(alpha = a,
                                  beta = b,
                                  gamma = g,
                                  avg_R0 = round(sim_out$avg_offspring, 2),
                                  var_R0 = round(sim_out$var_offspring, 2),
                                  avg_RR_peak = round(sim_out$avg_rr_peak, 2),
                                  avg_time_diff = round(sim_out$avg_time_diff, 2),
                                  avg_prop_cases = round(sim_out$avg_prop_infected, 2),
                                  avg_cases_control = round(sim_out$avg_control_infected, 2),
                                  avg_cases_tirs = round(sim_out$avg_tirs_infected, 2),
                                  avg_eff_estimate = round(sim_out$avg_pois_estimate, 2)
  )
}

big_pop_25ea <- bind_rows(grid_results)
saveRDS(big_pop_25ea, file = "output/sim_results/param_grid_search/big_pop_25ea.rds")