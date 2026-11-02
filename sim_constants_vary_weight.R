## trial population #############################################################
min_age <- 1 # youngest age
max_age <- 82 # oldest age

trial_start <- 26
trial_end <- 50 # 52 is the last week of the year but we don't have city incidence data for that??? so 50

# parameters
Gmu <- 2; Gs2 <- 1; max_gen_time <- 5 # 17 days?; 2 days?
# probability that the generation time = [index]

# Population size
n_clusters <- 50; n_per_cluster <- 2500; n_per_block <- 250

# Transmission parameters
block_R0 <- 0.9; cluster_R0 <- 0.9; gamma <- 2