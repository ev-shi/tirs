##############################################################################
# 1. libraries
##############################################################################
library(tidyverse)
library(survival)

##############################################################################
# 2. directory location, i.e., where all the populations are saved
##############################################################################
# infected_populations <- list.files("../../../../../../../../orange/sberube1/ev.shi/populations/11alphabeta/no_weights/base",
#                                    pattern = "\\.csv$", full.names = TRUE)

infected_populations <- list.files("../../../../../../../../orange/sberube1/ev.shi/populations/vary_weight/4_diff_weights_diff_sus/11_0039",
                                   pattern = "\\.csv$", full.names = TRUE)

# infected_populations <- list.files("output/populations",
#                                    pattern = "\\.csv$", full.names = TRUE)
# 
# infected_populations <- list.files("sim_results/populations/41767817_base_no_weight_0p3_0p3_1p0",
#                                    pattern = "\\.csv$", full.names = TRUE)

##############################################################################
# 3. function; get parameter set from the file name
##############################################################################
param_set_names <- function(name) {
  # e.g. blockRe0.3000_clusterRe0.3000_gamma0.5000_run001.csv
  m <- str_match(
    basename(name),
    "blockRe([0-9]*\\.?[0-9]+)_clusterRe([0-9]*\\.?[0-9]+)_gamma([0-9]*\\.?[0-9]+)_run[0-9]+\\.csv"
  )
  tibble(
    block = as.numeric(m[2]),
    cluster = as.numeric(m[3]),
    gamma = as.numeric(m[4])
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
  
  # TIRS ONLY
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
  
  # ENTIRE POP
  ## filter to no index cases
  infected_pop <- infected_pop[infected_pop$index_case != 1,]
  ## exit time definition
  infected_pop$exit_time <- ifelse(infected_pop$reported == 1, infected_pop$tau, 52)
  ## cox model
  cox_r <- coxph(Surv(rep(26, nrow(infected_pop)), exit_time, reported) ~ TIRS + cluster(Cluster), data = infected_pop)
  ## return efficacy and confidence interval
  hr <- exp(coef(cox_r)["TIRS"])
  ci <- exp(confint(cox_r)["TIRS",])
  
  # ADULTS
  ## already filtered to no index cases, adults only
  adults <- infected_pop[infected_pop$Age >= 18,]
  ## cox model
  cox_adults <- coxph(Surv(rep(26, nrow(adults)), exit_time, reported) ~ TIRS + cluster(Cluster), data = adults)
  ## return efficacy and confidence interval
  hr_adults <- exp(coef(cox_adults)["TIRS"])
  ci_adults <- exp(confint(cox_adults)["TIRS",])
  
  # KIDS
  ## already filtered to no index cases, kids only
  kids <- infected_pop[infected_pop$Age < 18,]
  ## cox model
  cox_kids <- coxph(Surv(rep(26, nrow(kids)), exit_time, reported) ~ TIRS + cluster(Cluster), data = kids)
  ## return efficacy and confidence interval
  hr_kids <- exp(coef(cox_kids)["TIRS"])
  ci_kids <- exp(confint(cox_kids)["TIRS",])
  
  return(list(eff_tirs = 1 - hazard_ratio,
              lower_tirs = 1 - conf_int[2],
              upper_tirs = 1 - conf_int[1],
              eff = 1 - hr,
              lower = 1 - ci[2],
              upper = 1 - ci[1],
              eff_adults = 1 - hr_adults,
              lower_adults = 1 - ci_adults[2],
              upper_adults = 1 - ci_adults[1],
              eff_kids = 1 - hr_kids,
              lower_kids = 1 - ci_kids[2],
              upper_kids = 1 - ci_kids[1]))
}

# function for poisson model
pois_model <- function(infected_pop) {
  
  # COMMUNITY ONLY
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
  irr_comm <- exp(coef(pois_reported_comm)["TIRS"])
  conf_int_comm <- exp(confint(pois_reported_comm)["TIRS", ])
  
  # ENTIRE POP
  ## number of cases per cluster and trial arm
  reported_by_cluster <- infected_pop %>% 
    group_by(Cluster, TIRS) %>% 
    summarize(reported_cases = sum(reported), .groups = "drop")
  # model
  pois_r <- glm(reported_cases ~ TIRS, family = poisson(link = "log"), data = reported_by_cluster)
  # return IRR and confint
  irr <- exp(coef(pois_r)["TIRS"])
  ci <- exp(confint(pois_r)["TIRS", ])
  
  # ADULTS
  ## number of cases per cluster and trial arm
  adults <- infected_pop %>% 
    filter(Age >= 18) %>% 
    group_by(Cluster, TIRS) %>% 
    summarize(reported_cases = sum(reported), .groups = "drop")
  # Does any ARM have 0 zero counts?
  if (any(
    adults %>%
    group_by(TIRS) %>%
    summarise(n_nonzero = sum(reported_cases > 0),
              .groups = "drop") %>%
    pull(n_nonzero) == 0
  )) {
    warning("An entire arm has no cases – IRR cannot be estimated.")
    return(list(eff = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  # model
  pois_adults <- glm(reported_cases ~ TIRS, family = poisson(link = "log"), data = adults)
  # return IRR and confint
  irr_adults <- exp(coef(pois_adults)["TIRS"])
  ci_adults <- exp(confint(pois_adults)["TIRS", ])
  
  # KIDS
  ## number of cases per cluster and trial arm
  kids <- infected_pop %>% 
    filter(Age < 18) %>% 
    group_by(Cluster, TIRS) %>% 
    summarize(reported_cases = sum(reported), .groups = "drop")
  # Does any ARM have 0 zero counts?
  if (any(
    kids %>%
    group_by(TIRS) %>%
    summarise(n_nonzero = sum(reported_cases > 0),
              .groups = "drop") %>%
    pull(n_nonzero) == 0
  )) {
    warning("An entire arm has no cases – IRR cannot be estimated.")
    return(list(eff = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  # model
  pois_kids <- glm(reported_cases ~ TIRS, family = poisson(link = "log"), data = kids)
  # return IRR and confint
  irr_kids <- exp(coef(pois_kids)["TIRS"])
  ci_kids <- exp(confint(pois_kids)["TIRS", ])
  
  return(list(eff_comm = 1 - irr_comm,
              lower_comm = 1 - conf_int_comm[2],
              upper_comm = 1 - conf_int_comm[1],
              eff = 1 - irr,
              lower = 1 - ci[2],
              upper = 1 - ci[1],
              eff_adults = 1 - irr_adults,
              lower_adults = 1 - ci_adults[2],
              upper_adults = 1 - ci_adults[1],
              eff_kids = 1 - irr_kids,
              lower_kids = 1 - ci_kids[2],
              upper_kids = 1 - ci_kids[1]))
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
  
  tibble(cox_eff_tirs = cox_reg$eff_tirs,
         cox_lower_tirs = cox_reg$lower_tirs,
         cox_upper_tirs = cox_reg$upper_tirs,
         cox_eff_full = cox_reg$eff,
         cox_lower_full = cox_reg$lower,
         cox_upper_full = cox_reg$upper,
         cox_eff_adults = cox_reg$eff_adults,
         cox_lower_adults = cox_reg$lower_adults,
         cox_upper_adults = cox_reg$upper_adults,
         cox_eff_kids = cox_reg$eff_kids,
         cox_lower_kids = cox_reg$lower_kids,
         cox_upper_kids = cox_reg$upper_kids,
         # pois
         pois_eff_comm = pois_reg$eff_comm,
         pois_lower_comm = pois_reg$lower_comm,
         pois_upper_comm = pois_reg$upper_comm,
         pois_eff_full = pois_reg$eff,
         pois_lower_full = pois_reg$lower,
         pois_upper_full = pois_reg$upper,
         pois_eff_adults = pois_reg$eff_adults,
         pois_lower_adults = pois_reg$lower_adults,
         pois_upper_adults = pois_reg$upper_adults,
         pois_eff_kids = pois_reg$eff_kids,
         pois_lower_kids = pois_reg$lower_kids,
         pois_upper_kids = pois_reg$upper_kids,
         # other
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
  unnest_longer(arm_counts) %>% 
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
  group_by(block, cluster, gamma) %>% 
  summarize(
    # Cox stats
    ## tirs
    cox_mean_tirs = mean(cox_eff_tirs[is_reasonable(cox_eff_tirs)]),
    cox_lwr_mean_tirs = mean(cox_lower_tirs[is_reasonable(cox_eff_tirs)]),
    cox_upr_mean_tirs = mean(cox_upper_tirs[is_reasonable(cox_eff_tirs)]),
    cox_median_tirs = median(cox_eff_tirs, na.rm = TRUE),
    cox_bias_tirs = cox_mean_tirs - true_effect, 
    cox_count_tirs = sum(is_reasonable(cox_eff_tirs), na.rm = TRUE),
    cox_ci_cntns_tru_tirs = sum(true_effect >= cox_lower_tirs & true_effect <= cox_upper_tirs, na.rm = TRUE),
    prop_cox_yes_tirs = cox_ci_cntns_tru_tirs / cox_count_tirs,
    ## entire population
    cox_mean_full = mean(cox_eff_full[is_reasonable(cox_eff_full)]),
    cox_lwr_mean_full = mean(cox_lower_full[is_reasonable(cox_eff_full)]),
    cox_upr_mean_full = mean(cox_upper_full[is_reasonable(cox_eff_full)]),
    cox_median_full = median(cox_eff_full, na.rm = TRUE),
    cox_bias_full = cox_mean_full - true_effect, 
    cox_count_full = sum(is_reasonable(cox_eff_full), na.rm = TRUE),
    cox_ci_cntns_tru_full = sum(true_effect >= cox_lower_full & true_effect <= cox_upper_full, na.rm = TRUE),
    prop_cox_yes_full = cox_ci_cntns_tru_full / cox_count_full,
    ## adults
    cox_mean_adults = mean(cox_eff_adults[is_reasonable(cox_eff_adults)]),
    cox_lwr_mean_adults = mean(cox_lower_adults[is_reasonable(cox_eff_adults)]),
    cox_upr_mean_adults = mean(cox_upper_adults[is_reasonable(cox_eff_adults)]),
    cox_median_adults = median(cox_eff_adults, na.rm = TRUE),
    cox_bias_adults = cox_mean_adults - true_effect, 
    cox_count_adults = sum(is_reasonable(cox_eff_adults), na.rm = TRUE),
    cox_ci_cntns_tru_adults = sum(true_effect >= cox_lower_adults & true_effect <= cox_upper_adults, na.rm = TRUE),
    prop_cox_yes_adults = cox_ci_cntns_tru_adults / cox_count_adults,
    ## kids
    cox_mean_kids = mean(cox_eff_kids[is_reasonable(cox_eff_kids)]),
    cox_lwr_mean_kids = mean(cox_lower_kids[is_reasonable(cox_eff_kids)]),
    cox_upr_mean_kids = mean(cox_upper_kids[is_reasonable(cox_eff_kids)]),
    cox_median_kids = median(cox_eff_kids, na.rm = TRUE),
    cox_bias_kids = cox_mean_kids - true_effect, 
    cox_count_kids = sum(is_reasonable(cox_eff_kids), na.rm = TRUE),
    cox_ci_cntns_tru_kids = sum(true_effect >= cox_lower_kids & true_effect <= cox_upper_kids, na.rm = TRUE),
    prop_cox_yes_kids = cox_ci_cntns_tru_kids / cox_count_kids,
    
    # Poisson stats
    ## community
    pois_mean_comm = mean(pois_eff_comm[is_reasonable(pois_eff_comm)]),
    pois_lwr_mean_comm = mean(pois_lower_comm[is_reasonable(pois_eff_comm)]),
    pois_upr_mean_comm = mean(pois_upper_comm[is_reasonable(pois_eff_comm)]),
    pois_median_comm = median(pois_eff_comm, na.rm = TRUE),
    pois_bias_comm = pois_mean_comm - true_effect, 
    pois_count_comm = sum(is_reasonable(pois_eff_comm), na.rm = TRUE),
    pois_ci_cntns_tru_comm = sum(true_effect >= pois_lower_comm & true_effect <= pois_upper_comm, na.rm = TRUE),
    prop_pois_yes_comm = pois_ci_cntns_tru_comm / pois_count_comm,
    ## entire population
    pois_mean_full = mean(pois_eff_full[is_reasonable(pois_eff_full)]),
    pois_lwr_mean_full = mean(pois_lower_full[is_reasonable(pois_eff_full)]),
    pois_upr_mean_full = mean(pois_upper_full[is_reasonable(pois_eff_full)]),
    pois_median_full = median(pois_eff_full, na.rm = TRUE),
    pois_bias_full = pois_mean_full - true_effect, 
    pois_count_full = sum(is_reasonable(pois_eff_full), na.rm = TRUE),
    pois_ci_cntns_tru_full = sum(true_effect >= pois_lower_full & true_effect <= pois_upper_full, na.rm = TRUE),
    prop_pois_yes_full = pois_ci_cntns_tru_full / pois_count_full,
    ## adults
    pois_mean_adults = mean(pois_eff_adults[is_reasonable(pois_eff_adults)]),
    pois_lwr_mean_adults = mean(pois_lower_adults[is_reasonable(pois_eff_adults)]),
    pois_upr_mean_adults = mean(pois_upper_adults[is_reasonable(pois_eff_adults)]),
    pois_median_adults = median(pois_eff_adults, na.rm = TRUE),
    pois_bias_adults = pois_mean_adults - true_effect, 
    pois_count_adults = sum(is_reasonable(pois_eff_adults), na.rm = TRUE),
    pois_ci_cntns_tru_adults = sum(true_effect >= pois_lower_adults & true_effect <= pois_upper_adults, na.rm = TRUE),
    prop_pois_yes_adults = pois_ci_cntns_tru_adults / pois_count_adults,
    ## kids
    pois_mean_kids = mean(pois_eff_kids[is_reasonable(pois_eff_kids)]),
    pois_lwr_mean_kids = mean(pois_lower_kids[is_reasonable(pois_eff_kids)]),
    pois_upr_mean_kids = mean(pois_upper_kids[is_reasonable(pois_eff_kids)]),
    pois_median_kids = median(pois_eff_kids, na.rm = TRUE),
    pois_bias_kids = pois_mean_kids - true_effect, 
    pois_count_kids = sum(is_reasonable(pois_eff_kids), na.rm = TRUE),
    pois_ci_cntns_tru_kids = sum(true_effect >= pois_lower_kids & true_effect <= pois_upper_kids, na.rm = TRUE),
    prop_pois_yes_kids = pois_ci_cntns_tru_kids / pois_count_kids,
    
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
    secondary_tirs = mean(avg_secondary_1, na.rm = TRUE),
    
    # IRR from Secondary Cases
    irr_secondary = 1 - (secondary_tirs / secondary_ctrl)
  ) %>% ungroup()


##############################################################################
# 8. write to csv
##############################################################################
saveRDS(full_df, "full_df.rds")
write_csv(summary_df, "summary_results.csv")