# Clean environment
rm(list=ls())

# Set up basic path
library(tidyverse)
library(dplyr)
library(here)
#setwd("..")
here()
set.seed(23)
batches <- 100

# Reading mcmc results
mcmc_alpha_s1_q1 <- read.table(
  here("code-Molina2023","output","mcmc.alpha-s1.txt"),
  sep=" ",
  header = TRUE
)

nSim <- as.integer(length(mcmc_alpha_s1_q1[,1])/batches) # 20000 were used to obtained the results presented in the article
if(nSim < 10000){burnin <- round(nSim/4)} else{burnin <- 5000}
if((nSim-burnin) < 500){selsamp <- (nSim-burnin)} else{selsamp <- 500}

nParam <- 4
mcmc_alpha_all_s1_q1 <- mcmc_alpha_s1_q1[mcmc_alpha_s1_q1$nSim %in% ((burnin+1):nSim),] # All values from the posterior
mcmc_alpha_sel_s1_q1 <- mcmc_alpha_all_s1_q1[1,] # Random selection of values from the posterior
for(i in 1:batches){
  index_selection <- sample(1:(nSim-burnin), selsamp)
  param_options <- mcmc_alpha_all_s1_q1[((i-1)*(nSim-burnin)+1):(i*(nSim-burnin)),]
  param_selection <- param_options[index_selection,]
  mcmc_alpha_sel_s1_q1 <- rbind(mcmc_alpha_sel_s1_q1, param_selection)
}
mcmc_alpha_sel_s1_q1 <- mcmc_alpha_sel_s1_q1[-1,]
write.table(mcmc_alpha_sel_s1_q1[,2:(nParam+1)], 
            file = here("code-Molina2023","output","posterior.params.alpha-s1.txt"), 
            row.names = FALSE, col.names = FALSE, sep = " ")


infection_times <- read.table(here("code-Molina2023","output","infection.times.alpha.txt"),
                              sep=" ", header = FALSE)
infection_times_2_good <- infection_times[1,]
infection_times_2_good <- infection_times_2_good[-1,]

for(i in 1:batches){
  x <- infection_times
  x <- x %>% slice(rep(i,selsamp))
  infection_times_2_good <- rbind(infection_times_2_good, x)
}

write.table(infection_times_2_good, file = here("code-Molina2023","output","infection.times.alpha.2.txt"),
            row.names = FALSE, col.names = FALSE, sep = " ")


