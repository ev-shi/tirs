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
                       alpha, beta, gamma, r, comm_weight) {
  
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
  # 70% of adults are susceptible
  
  kid_ids <- which(trial_pop$Age < 18) 
  immune_kids <- sample(kid_ids, size = round(0 * length(kid_ids)))
  trial_pop$susceptible[immune_kids] <- 0
  # 100% of kids are susceptible
  
  # # time spent in the area
  # trial_pop$block_weight <- ifelse(trial_pop$Age < 18,
  #                                  0.25, 0.4)
  # trial_pop$cluster_weight <- ifelse(trial_pop$Age < 18,
  #                                    0.15, 0.3)
  # # kids spend 40% of time in block/cluster, adults 70%
  
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
              lambda_jb <- lambda_jb + (1 - r) * alpha * gen_prob[dt]
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
              lambda_jc <- lambda_jc + (1 - r) * beta * gen_prob[dt]
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
          # hazard from a case in merida * cases rate in merida * prob gen time is dt
        }
      }
      
      # how long do they typically spend in the neighborhood?
      lambda_jb <- lambda_jb * ((1 - comm_weight)/2)
      lambda_jc <- lambda_jc * ((1 - comm_weight)/2)
      lambda_jm <- lambda_jm * comm_weight
      
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
              if (trial_pop$TIRS[i] == 1) { base <- (1 - r) * base }
              return(base)
            }
          }
          return(0)
        }), sapply(diff_block, function(i){
          if (trial_pop$infected[i] == 1) {
            dt <- t - trial_pop$tau[i]
            if (dt <= max_gen_time){ 
              base <- beta * gen_prob[dt]
              if (trial_pop$TIRS[i] == 1) { base <- base * (1 - r) }
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
                              size = round(0.6 * length(infected_TIRS_ids)))
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
                     beta, gamma, r, comm_weight, seed_offset = 0, chunk_id, sims_per_chunk = 10) {
  
  # each call should have distinct seed
  set.seed(12345 + seed_offset)
  
  real_df <- data.frame(
    time = 1:trial_end,
    case_rate = city_cases_by_week$cases_per_week / sum(merida_pop$POBTOT))
  
  for (s in 1:num_sims) {
    to_hund <- (chunk_id - 1) * sims_per_chunk + s
    
    infected_pop <- infect_pop(trial_pop, gen_prob, initially_infected, 
                               alpha = alpha, beta = beta, gamma = gamma, r = r,
                               comm_weight = comm_weight)
    
    # save this population
    suffix <- sprintf("_run%03d", to_hund)
    name <- paste0("blockRe", formatC((n_per_block - 1) * (1 - exp(-alpha)), 
                                      format = "f", digits = 2),
                   "_clusterRe", formatC((n_per_cluster - n_per_block) * 
                                           (1 - exp(-beta)), format = "f", 
                                         digits = 2),
                   "_gamma", formatC(gamma, format = "f", digits = 2),
                   "_comWeight", formatC(comm_weight, format = "f", digits = 2),
                   suffix)
    outfile <- file.path(here::here("output", "populations"),
                         paste0(name, ".csv"))
    write.csv(infected_pop, file = outfile)
  }
}