## trial population structure ###############################################################################
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

## generation time ##########################################################################################
gen_time <- function(Gmu, Gs2, max_gen_time) {
  shape <- Gmu^2 / Gs2
  scale <- Gs2 / Gmu
  
  i <- 0:max_gen_time
  
  genProb <- pgamma(i + 1, shape = shape, scale = scale) -
    pgamma(i,     shape = shape, scale = scale)
  
  return(genProb)
}

## infections ##############################################################################################
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
  immune_adults <- sample(adult_ids, size = round(0.3 * length(adult_ids)))
  trial_pop$susceptible[immune_adults] <- 0
  
  kid_ids <- which(trial_pop$Age < 18) 
  immune_kids <- sample(kid_ids, size = round(0.3 * length(kid_ids)))
  trial_pop$susceptible[immune_kids] <- 0
  
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
              # hazard from i infecting j * prob that gen time is dt
            }
          }
        }
      }
      
      # merida FOI
      for (z in trial_start:(t-1)) {
        dt <- t - z
        if (dt <= max_gen_time) {
          lambda_jm <- lambda_jm + gamma * (I_m[z] / sum(merida_pop$POBTOT)) * gen_prob[dt]
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
  
  ## TIRS kids
  trial_pop$TIRS_cohort <- 0
  trial_pop$reported <- 0
  
  # TIRS cohort 15% of total kids
  TIRS_ids <- sample(kid_ids, size = round(0.15 * length(kid_ids))) 
  trial_pop$TIRS_cohort[TIRS_ids] <- 1
  
  # reporting
  infected_TIRS_ids <- which(trial_pop$infected == 1 & 
                               trial_pop$TIRS_cohort == 1)
  reported_tirs_ids <- sample(infected_TIRS_ids, 
                              size = round(0.8 * length(infected_TIRS_ids)))
  infected_sinave_ids <- which(trial_pop$infected == 1 & 
                                 trial_pop$TIRS_cohort == 0)
  reported_sinave_ids <- sample(infected_sinave_ids, 
                                size = round(0.1 * length(infected_sinave_ids)))
  trial_pop$reported[reported_tirs_ids] <- 1
  trial_pop$reported[reported_sinave_ids] <- 1
  
  return(trial_pop = trial_pop)
}

## extra ##################################################################################################
format <- function(x) gsub("\\.", "p", sprintf("%.2f", x)) 

## run simulation ##########################################################################################
run_sims <- function(num_sims, trial_pop, gen_prob, initially_infected, alpha, 
                     beta, gamma, r) {
  
  sim_list <- list()
  offspring_index_control <- c()
  offspring_index_tirs <- c()
  rr_peak <- c()
  time_diff_at_peaks <- c()
  control_list <- list()
  tirs_list <- list()
  prop_infected <- c()
  prop_infected_sus <- c()
  control_infected <- c()
  tirs_infected <- c()
  pois_efficacy <- c()
  reported_pois <- c()
  pois_popoff_eff <- c()
  pois_ptoff_eff <- c()
  cox_efficacy <- c()
  reported_cox <- c()
  source_community <- c()
  source_block <- c()
  source_cluster <- c()
  cox_tirs <- c()
  cox_tirs_r <- c()
  pois_community <- c()
  pois_community_r <- c()
  
  real_df <- data.frame(
    time = 1:trial_end,
    case_rate = city_cases_by_week$cases_per_week / sum(merida_pop$POBTOT))
  
  for (s in 1:num_sims) {
    infected_pop <- infect_pop(trial_pop, gen_prob, initially_infected, 
                               alpha = alpha, beta = beta, gamma = gamma, r = r)
    
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
    index_cases <- infected_pop %>%
      filter(index_case == 1) %>%
      mutate(
        offspring = sapply(ID, function(id) {
          sum(infected_pop$infector == id, na.rm = TRUE)
        })
      ) %>%
      group_by(TIRS) %>%
      summarize(mean_offspring = mean(offspring), .groups = "drop")
    
    offspring_index_control[s] <- ifelse(
      any(index_cases$TIRS == 0),
      index_cases$mean_offspring[index_cases$TIRS == 0],
      NA
    )
    
    offspring_index_tirs[s] <- ifelse(
      any(index_cases$TIRS == 1),
      index_cases$mean_offspring[index_cases$TIRS == 1],
      NA
    )
    
    # find peak incidence in simulation
    study_peak <- max(incidence$case_rate)
    real_peak <- max(real_df$case_rate)
    
    rr_peak[s] <- study_peak / real_peak # difference in magnitude
    time_diff_at_peaks[s] <- incidence$time[which.max(incidence$case_rate)] - 
      real_df$time[which.max(real_df$case_rate)]
    
    # proportion of population infected at end of season
    prop_infected[s] <- sum(incidence$cases) / (n_clusters*n_per_cluster)
    
    # proportion of population infected at end of season out of susceptible pop
    prop_infected_sus[s] <- sum(incidence$cases) / sum(infected_pop$susceptible)
    
    # num infected by trial arm
    control_infected[s] <- sum(control$cases)
    tirs_infected[s] <- sum(tirs$cases)
    
    # Poisson regression (no offset) - true efficacy
    cases_by_cluster <- infected_pop %>% 
      group_by(Cluster, TIRS) %>% 
      summarize(num_cases = sum(infected), .groups = "drop")
    
    pois_fit <- glm(num_cases ~ TIRS, family = poisson(link = "log"), 
                    data = cases_by_cluster)
    pois_efficacy[s] <- 1 - exp(coef(pois_fit)["TIRS"])
    
    # Poisson regression (no offset) - REPORTED efficacy
    reported_by_cluster <- infected_pop %>% 
      group_by(Cluster, TIRS) %>% 
      summarize(reported_cases = sum(reported), .groups = "drop")
    
    pois_reported <- glm(reported_cases ~ TIRS, family = poisson(link = "log"),
                         data = reported_by_cluster)
    reported_pois[s] <- 1 - exp(coef(pois_reported)["TIRS"])
    
    # Poisson regression (no offset, community effect) - all infections
    community_group <- infected_pop[infected_pop$TIRS_cohort == 0,]
    comm_cases_by_cluster <- community_group %>% 
      group_by(Cluster, TIRS) %>% 
      summarize(num_cases = sum(infected), .groups = "drop")
    
    pois_fit_comm <- glm(num_cases ~ TIRS, family = poisson(link = "log"), 
                         data = comm_cases_by_cluster)
    pois_community[s] <- 1 - exp(coef(pois_fit_comm)["TIRS"])
    
    # Poisson regression (no offset, community effect) - REPORTED efficacy
    comm_reported_by_cluster <- community_group %>% 
      group_by(Cluster, TIRS) %>% 
      summarize(reported_cases = sum(reported), .groups = "drop")
    
    pois_reported_comm <- glm(reported_cases ~ TIRS, family = poisson(link = "log"),
                              data = comm_reported_by_cluster)
    pois_community_r[s] <- 1 - exp(coef(pois_reported_comm)["TIRS"])
    
    # Poisson regression (offset by population) - true efficacy
    cases_pop_by_cluster <- cbind(cases_by_cluster, 
                                  n_per_cluster = rep(n_per_cluster, 
                                                      nrow(cases_by_cluster)))
    
    pois_ofit_pop <- glm(num_cases ~ TIRS + offset(log(n_per_cluster)), 
                         family = poisson(link = "log"), 
                         data = cases_pop_by_cluster)
    pois_popoff_eff[s] <- 1 - exp(coef(pois_ofit_pop)["TIRS"])
    
    # Poisson regression (offset by person-time?) - true efficacy
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
    
    # Cox model (as done in paper) - true efficacy
    # first remove index case
    cox_dat <- infected_pop[infected_pop$index_case != 1,]
    # itt.cox <- coxph(Surv(entry, exit, Pos_All) ~ arm + cluster(cluster),
    #                   data = data)
    cox_dat$exit_time <- ifelse(is.na(cox_dat$tau), trial_end, cox_dat$tau)
    cox_fit <- coxph(Surv(rep(26, nrow(cox_dat)), 
                          exit_time, infected) ~ TIRS + cluster(Cluster),
                     data = cox_dat)
    cox_efficacy[s] <- 1 - exp(coef(cox_fit))
    
    # Cox model, using only reported cases
    # first remove index cases
    reported_dat <- infected_pop[infected_pop$index_case != 1,]
    reported_dat$exit_time <- ifelse(reported_dat$reported == 1, 
                                     reported_dat$tau, trial_end)
    cox_reported <- coxph(Surv(rep(26, nrow(reported_dat)),
                               exit_time, reported) ~ TIRS + cluster(Cluster),
                          data = reported_dat)
    reported_cox[s] <- 1 - exp(coef(cox_reported))
    
    # Cox model for the all TIRS trial kids
    # first remove index cases
    tirs_dat <- infected_pop[infected_pop$index_case != 1,]
    # now only TIRS population
    tirs_dat <- tirs_dat[tirs_dat$TIRS_cohort == 1,]
    tirs_dat$exit_time <- ifelse(is.na(tirs_dat$tau), trial_end, tirs_dat$tau)
    cox_tirsfit <- coxph(Surv(rep(26, nrow(tirs_dat)),
                              exit_time, infected) ~ TIRS + cluster(Cluster),
                         data = tirs_dat)
    cox_tirs[s] <- 1 - exp(coef(cox_tirsfit))    
    
    # Cox model for the reported TIRS trial kids, just like in paper
    # first remove index cases
    tirs_dat_r <- infected_pop[infected_pop$index_case != 1,]
    # now only TIRS population
    tirs_dat_r <- tirs_dat_r[tirs_dat_r$TIRS_cohort == 1,]
    tirs_dat_r$exit_time <- ifelse(tirs_dat_r$reported == 1, 
                                   tirs_dat_r$tau, trial_end)
    cox_tirsfit_r <- coxph(Surv(rep(26, nrow(tirs_dat_r)),
                                exit_time, reported) ~ TIRS + cluster(Cluster),
                           data = tirs_dat_r)
    cox_tirs_r[s] <- 1 - exp(coef(cox_tirsfit_r))    
    
    # infector source: num infected per source / total pop infected
    source_community[s] <- sum(infected_pop$source == 0, na.rm = TRUE) / 
      sum(infected_pop$infected)
    source_block[s] <- sum(infected_pop$source == 1, na.rm = TRUE) / 
      sum(infected_pop$infected)
    source_cluster[s] <- sum(infected_pop$source == 2, na.rm = TRUE) / 
      sum(infected_pop$infected)
  }
  
  # R_0
  avg_offspring_c <- mean(offspring_index_control, na.rm = TRUE)
  var_offspring_c <- var(offspring_index_control, na.rm = TRUE)
  
  avg_offspring_t <- mean(offspring_index_tirs, na.rm = TRUE)
  var_offspring_t <- var(offspring_index_tirs, na.rm = TRUE)
  
  # rate ratio between study population and merida
  avg_rr_peak <- mean(rr_peak, na.rm = TRUE)
  var_rr_peak <- var(rr_peak, na.rm = TRUE)
  avg_time_diff <- mean(time_diff_at_peaks, na.rm = TRUE)
  
  # proportion of population infected
  avg_prop_infected <- mean(prop_infected, na.rm = TRUE)
  # proportion of susceptibles infected
  avg_prop_inf_sus <- mean(prop_infected_sus, na.rm = TRUE)
  
  # avg infected by trial arm
  avg_control_infected <- mean(control_infected, na.rm = TRUE)
  avg_tirs_infected <- mean(tirs_infected, na.rm = TRUE)
  
  # avg efficacy estimates
  avg_pois_estimate <- mean(pois_efficacy, na.rm = TRUE)
  avg_pois_pop_estimate <- mean(pois_popoff_eff, na.rm = TRUE)
  avg_pois_pt_estimate <- mean(pois_ptoff_eff, na.rm = TRUE)
  avg_cox_estimate <- mean(cox_efficacy, na.rm = TRUE)
  
  avg_pois_reported_est <- mean(reported_pois, na.rm = TRUE)
  avg_cox_reported_est <- mean(reported_cox, na.rm = TRUE)
  
  avg_pois_comm <- mean(pois_community, na.rm = TRUE)
  avg_pois_comm_r <- mean(pois_community_r, na.rm = TRUE)
  avg_cox_tirs <- mean(cox_tirs, na.rm = TRUE)
  avg_cox_tirs_r <- mean(cox_tirs_r, na.rm = TRUE)
  
  # avg proportion of infector sources
  avg_infected_within_block <- mean(source_block, na.rm = TRUE)
  avg_infected_within_cluster <- mean(source_cluster, na.rm = TRUE)
  avg_infected_from_community <- mean(source_community, na.rm = TRUE)
  
  return(list(avg_offspring_c = avg_offspring_c, 
              var_offspring_c = var_offspring_c,
              avg_offspring_t = avg_offspring_t, 
              var_offspring_t = var_offspring_t,
              avg_rr_peak = avg_rr_peak, var_rr_peak = var_rr_peak,
              avg_time_diff = avg_time_diff, 
              avg_prop_infected = avg_prop_infected,
              avg_prop_inf_sus = avg_prop_inf_sus,
              avg_control_infected = avg_control_infected,
              avg_tirs_infected = avg_tirs_infected,
              avg_pois_estimate = avg_pois_estimate,
              avg_pois_reported_est = avg_pois_reported_est,
              avg_pois_pop_estimate = avg_pois_pop_estimate,
              avg_pois_pt_estimate = avg_pois_pt_estimate,
              avg_cox_estimate = avg_cox_estimate,
              avg_cox_reported_est = avg_cox_reported_est,
              avg_pois_comm = avg_pois_comm,
              avg_pois_comm_r = avg_pois_comm_r,
              avg_cox_tirs = avg_cox_tirs,
              avg_cox_tirs_r = avg_cox_tirs_r,
              avg_infected_within_block = avg_infected_within_block,
              avg_infected_within_cluster = avg_infected_within_cluster,
              avg_infected_from_community = avg_infected_from_community))
}