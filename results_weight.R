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
# 4. functions; cox and poisson regressions - get point estimate and CI
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
# 5. go through all files and store the results
##############################################################################
results_list <- map_dfr(infected_populations, function(name) {
  # read CSV
  pop <- read_csv(name, show_col_types = FALSE)
  
  # regression models
  cox_reg <- cox_model(pop)
  pois_reg <- pois_model(pop)
  
  tibble(cox_eff = cox_reg$eff,
         cox_lower = cox_reg$lower,
         cox_upper = cox_reg$upper,
         pois_eff = pois_reg$eff,
         pois_lower = pois_reg$lower,
         pois_upper = pois_reg$upper)
}) 

##############################################################################
# 6. combine results with names
##############################################################################
full_df <- bind_cols(param_set, results_list)

# only keep reasonable estimates
is_reasonable <- function(x) !is.na(x) & abs(x) <= 1e6

# true treatment effect
true_effect <- 0.039

# summaries
summary_df <- full_df %>% 
  group_by(block, cluster, gamma) %>% 
  summarize(
    # Cox stats
    cox_mean = mean(cox_eff[is_reasonable(cox_eff)]),
    cox_lwr_mean = mean(cox_lower[is_reasonable(cox_eff)]),
    cox_upr_mean = mean(cox_upper[is_reasonable(cox_eff)]),
    cox_median = median(cox_eff, na.rm = TRUE),
    cox_bias = true_effect - cox_mean, 
    cox_count = sum(is_reasonable(cox_eff), na.rm = TRUE),
    cox_ci_cntns_tru = sum(true_effect >= cox_lower & true_effect <= cox_upper, na.rm = TRUE),
    prop_cox_yes = cox_ci_cntns_tru / cox_count,
    
    # Poisson stats
    pois_mean = mean(pois_eff[is_reasonable(pois_eff)]),
    pois_lwr_mean = mean(pois_lower[is_reasonable(pois_eff)]),
    pois_upr_mean = mean(pois_upper[is_reasonable(pois_eff)]),
    pois_median = median(pois_eff, na.rm = TRUE),
    pois_bias = true_effect - pois_mean, 
    pois_count = sum(is_reasonable(pois_eff), na.rm = TRUE),
    pois_ci_cntns_tru = sum(0.039 >= pois_lower & 0.039 <= pois_upper, na.rm = TRUE),
    prop_pois_yes = pois_ci_cntns_tru / pois_count,
    .groups = "drop"
  ) 


##############################################################################
# 7. write to csv
##############################################################################
write_csv(summary_df, "reg_estimations.csv")

# progress
cat(sprintf("Processed %d runs; produced %d parameter‑set rows.\n",
            nrow(full_df), nrow(summary_df)))