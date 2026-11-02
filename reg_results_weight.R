##############################################################################
# 1. libraries
##############################################################################
library(tidyverse)
library(survival)

##############################################################################
# 2. directory location, i.e., where all the populations are saved
##############################################################################
infected_populations <- list.files("output/populations", 
                                   pattern = "\\.csv$", full.names = TRUE)

# infected_populations <- list.files("sim_results/populations/41767817_base_no_weight_0p3_0p3_1p0", 
#                                    pattern = "\\.csv$", full.names = TRUE)

##############################################################################
# 3. function; get parameter set from the file name
##############################################################################
param_set_names <- function(name) {
  # e.g. blockRe0.3000_clusterRe0.3000_gamma0.5000_run001.csv
  m <- str_match(
    basename(name),
    "blockRe([0-9]*\\.?[0-9]+)_clusterRe([0-9]*\\.?[0-9]+)_gamma([0-9]*\\.?[0-9]+)_comWeightKids([0-9]*\\.?[0-9]+)_comWeightAdults([0-9]*\\.?[0-9]+)_run[0-9]+\\.csv"
  )
  tibble(
    block = as.numeric(m[2]),
    cluster = as.numeric(m[3]),
    gamma = as.numeric(m[4]),
    commWeightKids = as.numeric(m[5]),
    commWeightAdults = as.numeric(m[6])
  )
}

param_set <- infected_populations %>% 
  map_dfr(param_set_names)

##############################################################################
# 4. functions; some summary statistics
##############################################################################
cases_tables <- function(infected_pop) {
  cases_arm <- infected_pop %>% 
    group_by(TIRS) %>% 
    summarize(true_cases_arm = sum(infected, na.rm = TRUE),
              reported_cases_arm = sum(reported, na.rm = TRUE),
              .groups = "drop")
  
  cases_cluster <- infected_pop %>% 
    group_by(Cluster, TIRS) %>% 
    summarize(true_cases_cluster = sum(infected, na.rm = TRUE),
              reported_cases_cluster = sum(reported, na.rm = TRUE),
              .groups = "drop")
  
  list(cases_arm = cases_arm,
       cases_cluster = cases_cluster)
}

secondary_cases <- function(infected_pop) {
  
  # all secondary cases from each index case
  secondary <- infected_pop %>% 
    filter(index_case == 1) %>% 
    group_by(ID, TIRS) %>% 
    summarize(secondary_cases = sum(infected_pop$infector == ID, na.rm = TRUE),
              .groups = "drop")
  
  # average by arm
  avg_secondary_by_arm <- secondary %>% 
    group_by(TIRS) %>% 
    summarize(avg_secondary = mean(secondary_cases, na.rm = TRUE),
              .groups = "drop")
  
  list(secondary = secondary,
       avg_secondary_by_arm = avg_secondary_by_arm)
}

infector_sources <- function(infected_pop) {
  source_community <- sum(infected_pop$source == 0, na.rm = TRUE) /
    sum(infected_pop$infected)
  source_block <- sum(infected_pop$source == 1, na.rm = TRUE) /
    sum(infected_pop$infected)
  source_cluster <- sum(infected_pop$source == 2, na.rm = TRUE) /
    sum(infected_pop$infected)
  
  list(source_community = source_community,
       source_block = source_block,
       source_cluster = source_cluster)
}

prop_infected <- function(infected_pop) {
  ttl_infec <- sum(infected_pop$infected == 1) / (50 * 2500)
  rprt_infec <- sum(infected_pop$reported == 1) / (50 * 2500)
  
  list(ttl_infec = ttl_infec,
       rprt_infec = rprt_infec)
}

##############################################################################
# 5. functions; cox and poisson regressions - get point estimate and CI
##############################################################################
# function for cox model
cox_model <- function(infected_pop) {
  
  # filter to only TIRS population (and no index cases)
  tirs_dat_r <- infected_pop[infected_pop$index_case != 1,]
  tirs_dat_r <- tirs_dat_r[tirs_dat_r$TIRS_cohort == 1,]
  
  # exit time definition
  # if they reported then exit time is the time they reported
  tirs_dat_r$exit_time <- ifelse(tirs_dat_r$reported == 1, tirs_dat_r$tau, 52)
  
  # cox model
  cox_tirsfit_r <- coxph(
    Surv(rep(26, nrow(tirs_dat_r)), exit_time, reported) ~ 
      TIRS + cluster(Cluster), 
    data = tirs_dat_r)
  
  # return efficacy and confidence interval
  hazard_ratio <- exp(coef(cox_tirsfit_r)["TIRS"])
  conf_int <- exp(confint(cox_tirsfit_r)["TIRS",])
  return(list(eff = 1 - hazard_ratio,
              lower = 1 - conf_int[2],
              upper = 1 - conf_int[1]))
}

# function for poisson model
pois_model <- function(infected_pop) {
  
  # filter to only community, no TIRS kids
  community_group <- infected_pop[infected_pop$TIRS_cohort == 0,]
  
  # number of cases per cluster and trial arm
  comm_reported_by_cluster <- community_group %>% 
    group_by(Cluster, TIRS) %>% 
    summarize(reported_cases = sum(reported), .groups = "drop")
  
  # Does any ARM have 0 zero counts?
  if (any(
    comm_reported_by_cluster %>%
    group_by(TIRS) %>%
    summarise(n_nonzero = sum(reported_cases > 0),
              .groups = "drop") %>%
    pull(n_nonzero) == 0
  )) {
    warning("An entire arm has no cases – IRR cannot be estimated.")
    return(list(eff = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  
  # model
  pois_reported_comm <- glm(reported_cases ~ TIRS, family = poisson(link = "log"),
                            data = comm_reported_by_cluster)
  
  # return IRR and confint
  irr <- exp(coef(pois_reported_comm)["TIRS"])
  conf_int <- exp(confint(pois_reported_comm)["TIRS", ])
  return(list(eff = 1 - irr,
              lower = 1 - conf_int[2],
              upper = 1 - conf_int[1]))
}

##############################################################################
# 6. go through all files and store the results
##############################################################################
results_list <- map_dfr(infected_populations, function(name) {
  # read CSV
  pop <- read_csv(name, show_col_types = FALSE)
  
  # summary stats
  case_tallies <- cases_tables(pop)
  secondary <- secondary_cases(pop)
  source <- infector_sources(pop)
  prop_infec <- prop_infected(pop)
  
  # regression models
  cox_reg <- cox_model(pop)
  pois_reg <- pois_model(pop)
  
  tibble(cox_eff = cox_reg$eff,
         cox_lower = cox_reg$lower,
         cox_upper = cox_reg$upper,
         pois_eff = pois_reg$eff,
         pois_lower = pois_reg$lower,
         pois_upper = pois_reg$upper,
         arm_counts = list(case_tallies$cases_arm), # one row, 2 cols 
         cluster_counts = list(case_tallies$cases_cluster), # 1 row per cluster
         secondary_by_index = list(secondary$secondary),
         secondary_by_arm = list(secondary$avg_secondary_by_arm),
         source_comm = list(source$source_community),
         source_block = list(source$source_block),
         source_cluster = list(source$source_cluster),
         total_infected = list(prop_infec$ttl_infec),
         reported_infected = list(prop_infec$rprt_infec))
}) 

##############################################################################
# 7. combine results with names
##############################################################################
full_df <- bind_cols(param_set, results_list)

# only keep reasonable estimates
is_reasonable <- function(x) !is.na(x) & abs(x) <= 1e6

# true treatment effect
true_effect <- 0.25

# summaries
summary_df <- full_df %>% 
  nnest_longer(arm_counts) %>% 
  unnest_wider(arm_counts) %>% 
  pivot_wider(names_from = TIRS,
              values_from = c(true_cases_arm,
                              reported_cases_arm),
              names_glue = "{.value}_{TIRS}") %>% 
  unnest_longer(secondary_by_arm) %>% 
  unnest_wider(secondary_by_arm) %>% 
  pivot_wider(names_from = TIRS,
              values_from = avg_secondary,
              names_glue = "{.value}_{TIRS}") %>% 
  group_by(block, cluster, gamma, commWeightKids, commWeightAdults) %>% 
  summarize(
    # Cox stats
    cox_mean = mean(cox_eff[is_reasonable(cox_eff)]),
    cox_lwr_mean = mean(cox_lower[is_reasonable(cox_eff)]),
    cox_upr_mean = mean(cox_upper[is_reasonable(cox_eff)]),
    cox_median = median(cox_eff, na.rm = TRUE),
    cox_bias = cox_mean - true_effect, 
    cox_count = sum(is_reasonable(cox_eff), na.rm = TRUE),
    cox_ci_cntns_tru = sum(true_effect >= cox_lower & true_effect <= cox_upper, na.rm = TRUE),
    prop_cox_yes = cox_ci_cntns_tru / cox_count,
    
    # Poisson stats
    pois_mean = mean(pois_eff[is_reasonable(pois_eff)]),
    pois_lwr_mean = mean(pois_lower[is_reasonable(pois_eff)]),
    pois_upr_mean = mean(pois_upper[is_reasonable(pois_eff)]),
    pois_median = median(pois_eff, na.rm = TRUE),
    pois_bias = pois_mean - true_effect, 
    pois_count = sum(is_reasonable(pois_eff), na.rm = TRUE),
    pois_ci_cntns_tru = sum(true_effect >= pois_lower & true_effect <= pois_upper, na.rm = TRUE),
    prop_pois_yes = pois_ci_cntns_tru / pois_count,
    
    # Source Proportions
    source_comm = mean(as.numeric(source_comm), na.rm = TRUE),
    source_block = mean(as.numeric(source_block), na.rm = TRUE),
    source_cluster = mean(as.numeric(source_cluster), na.rm = TRUE),
    
    # Reported Infected
    total_infected = mean(as.numeric(total_infected), na.rm = TRUE),
    reported_infected = mean(as.numeric(reported_infected), na.rm = TRUE),
    
    # Case Counts
    true_cases_ctrl = mean(true_cases_arm_0, na.rm = TRUE),
    true_cases_tirs = mean(true_cases_arm_1, na.rm = TRUE),
    rprtd_cases_ctrl = mean(reported_cases_arm_0, na.rm = TRUE),
    rprtd_cases_tirs = mean(reported_cases_arm_1, na.rm = TRUE),
    
    # Secondary Cases
    secondary_ctrl = mean(avg_secondary_0, na.rm = TRUE),
    secondary_tirs = mean(avg_secondary_1, na.rm = TRUE)
  ) %>% ungroup()


##############################################################################
# 8. write to csv
##############################################################################
saveRDS(full_df, "full_df.rds")
write_csv(summary_df, "summary_results.csv")