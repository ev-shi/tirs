#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "headerTIRS.h"

/* Reads global variables, the number of individuals, and the sigmas for the MCMC parameters */
void ReadFiles1(char *address_global, Individuals *indiv, Global *global, Incidence *incidence, Parameters *param){
  
  FILE *globalFile;
  globalFile = fopen(address_global, "r");
  if(globalFile == NULL){
    fprintf(stderr, "\nERROR opening the file %s for reading\n", address_global);
    exit(2);
  }
  while(!feof(globalFile)){fscanf(globalFile, "%d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf", 
        &indiv->Nindiv,
        &global->Nclusters,
        &global->Nweeks,
        &global->Week0,
        &global->WeekL,
        &global->MAX_GENTIME,
        &incidence->POBTOT,
        &param->sigma_alpha,
        &param->sigma_beta,
        &param->sigma_gamma,
        &param->sigma_r,
        &param->sigma_Gta,
        &param->sigma_Gtb);
  }
  fclose(globalFile);
  
  // initialization of the parameters
  param->alpha = 0.; param->beta = 0.; param->gamma = 0.; param->Gta = 0.; param->Gtb = 0.;
}

/* Reads individuals data and city incidence */
void ReadFiles2(char *address_indiv, char *address_incidence, Individuals *indiv,
                Global global, Incidence *incidence){
  int i;
  
  // Data on individuals
  i = 0;
  FILE *indivFile;
  indivFile = fopen(address_indiv, "r");
  if(indivFile == NULL){
    fprintf(stderr, "\nERROR opening the file %s for reading\n", address_indiv);
    exit(2);
  }
  while(!feof(indivFile)){
    fscanf(indivFile,"%d %d %d %d %d", 
           &indiv->IDperson[i], 
           &indiv->IDblock[i],
           &indiv->IDcluster[i],
           &indiv->TIRS[i],
           &indiv->symptonset[i]); 
           i++;
  }
  fclose(indivFile);
  
  // Incidence in city
  i = 0;
  FILE *commFile;
  commFile = fopen(address_incidence, "r");
  if(commFile == NULL){
    fprintf(stderr, "\nERROR opening the file %s for reading\n", address_incidence);
    exit(2);
  }
  
  while(!feof(commFile)){
    fscanf(commFile,"%d %lf", &incidence->date[i], &incidence->incidence[i]);
    i++;
  }
  fclose(commFile);
  
  // scale by max incidence
  double max_incidence = 0.; // find max of all incidence vals over the weeks
  for(i=0; i<global.Nweeks; i++){
    if(incidence->incidence[i] > max_incidence){
      max_incidence = incidence->incidence[i];
    }
  }
  for(i=0; i<global.Nweeks; i++){ // scale incidence by max
    incidence->incidence[i] = incidence->incidence[i] / max_incidence;
  }
}

/* Initialize complete data for model with observed symptom onset times.*/
void initInfecData(Individuals *indiv){
  
  for(int i=0; i < indiv->Nindiv; i++){
    
    double incubation = 3.0 + gsl_rng_uniform(RNG) * (10.0 - 3.0); // Uniform(3, 10)
    
    // infection time = onset - incubation
    double inf_time = (double)indiv->symptonset[i] - incubation;
    indiv->infecTime[i] = inf_time;
  }
}

/* Calculate instantaneous FOI for given individual j and week t*/
long double foi_at_t(int j, int t, Individuals *indiv, Global global, Incidence incidence, Parameters param){
  
  // for convenience
  const int N = indiv->Nindiv;
  const int max_genTime = global.MAX_GENTIME;
  const int *block = indiv->IDblock;
  const int *cluster = indiv->IDcluster;
  
  double lambda = 0.0;
  
  // contribution from previously infected individuals i in the same block and cluster as j
  for (int i=0; i<N; i++) {
    
    // for indiv i to exert FOI on j, i must have been infected prior to j
    double ti = indiv->infecTime[i];
    if (ti >= (double)t) continue; 
      
    // time between infections i and j
    int dt = t - (int)ti;
    // generation time must be less than max_genTime
    if (dt <= 0 || dt > max_genTime) continue;
      
    // check for same block/cluster
    bool same_block = (block[i] == block[j]);
    bool same_cluster = ((cluster[i] == cluster[j]) && (block[i] != block[j]));
      
    // base hazard before adjusting for treatment
    double base = 0.0;
    if (same_block) base = param.alpha * param.genDensity[dt];
    else if (same_cluster) base = param.beta * param.genDensity[dt];
    else continue;
      
    // if in treatment arm, multiply by reduction to hazard
    if (indiv->TIRS[j] == 1) base *= (1.0 - param.r);
    
    // add to total FOI  
    if (same_block) lambda += base;
    else lambda += base;
      
    } // end block and cluster for loop
    
    // Merida/community FOI
    for (int z = 0; z <= t; z++) {
      
      // check generation time
      int dt = t - z;
      if (dt < 0 || dt > max_genTime) continue;
      
      lambda += param.gamma * (incidence.incidence[z] / incidence.POBTOT) * param.genDensity[dt];
      
    } // end community FOI
    
  return lambda;
}

/* Calculate Likelihood */
long double LogLikelihood(Individuals *indiv, Global global, Incidence incidence, Parameters param){
  
  long double loglike = 0.;

  // for every individual case
  for(int j=0; j<indiv->Nindiv; j++){
    
    // for convenience
    double tj = indiv->infecTime[j]; // infection time of individual j
    
    // calculate instantaneous FOI at tau
    double lambda_tj = foi_at_t(j, (int)tj, indiv, global, incidence, param);
    
    // calculate exponential (cumulative hazard) up to day before infection
    double cum_haz = 0.0;
    for (int t = 1; t<(int)tj; t++){
      cum_haz += foi_at_t(j, t, indiv, global, incidence, param);
    }
    
    loglike += log((long double) lambda_tj) - (long double) cum_haz;
    
  } // close indiv case for loop
  
  return loglike;
}

/* Update unknown infection times using Metropolis-Hastings step*/
void updateInfecData(Individuals *indiv, Global global, Incidence incidence, Parameters param){
  
  // randomly sample 20% of infection times to explore uncertainty in infection time of observed cases
  double update_prob = 0.20;
  
  for(int i=0; i<indiv->Nindiv; i++){
    
    if (gsl_rng_uniform(RNG) > update_prob) {
      continue; // time not updated for person i this iteration
    }
    
    double old_time = indiv->infecTime[i]; // old time
    double old_loglike = LogLikelihood(indiv, global, incidence, param);
    
    // propose new infection time for individual i
    double incubation = 3.0 + gsl_rng_uniform(RNG) * (10.0 - 3.0); // Uniform(3, 10)
    double new_time = (double)indiv->symptonset[i] - incubation;
    indiv->infecTime[i] = new_time;
    
    // loglikelihood of proposed state
    double new_loglike = LogLikelihood(indiv, global, incidence, param);
    
    // Metropolis-Hastings acceptance test
    double log_uniform = log(gsl_rng_uniform(RNG));
    if(log_uniform >= (new_loglike - old_loglike)) {
      
      // reject and continue with old time
      indiv->infecTime[i] = old_time;
      
    } // else keep new value
    
  } // end for loop
} 

/* Draw from a reflected normal distribution about 0 */
double rnorm_reflect_pos(double mean, double sd){
  
  // draw raw value 
  double x = gsl_ran_gaussian(RNG, sd) + mean;
  
  // reflect around 0 if needed
  if (x < 0.0){
    x = -x;
  }
  
  // don't return 0
  if (x == 0.0){
    x += DBL_MIN; // in C, minimum normalized positive value...used as underflow here?
  }
  
  return x;
}

/* Compute vector for generation time density, evaluated from 0 to MAX_GENTIME with given parameters */
void GenDensity(Parameters *param, Global global){
  
  // loop over every week up to maximum generation time
  for(int i=0; i<=global.MAX_GENTIME; i++){
    param->genDensity[i] = 0.;
    
    // taking small slice of the 1-week interval to approximate the integral of the Gamma pdf over that week
    for(int j=0; j<steps; j++){
      // got same error as carla
      if (param->Gta <= 0) param->Gta = 1e-6;
      if (param->Gtb <= 0) param->Gtb = 1e-6;
      param->genDensity[i] += gsl_ran_gamma_pdf((double)(i+j/steps), (param->Gta), (param->Gtb))/steps;
    }
  }
  
}

/* Accept candidate */
int stageAccept(long double candidate, long double old){
  
  // log likelihood ratio of candidate and old states
  long double lr = candidate - old;
  
  // Draw a uniform number u in (0,1) and log‑transform (better for small u)
  long double u = log10(gsl_ran_flat(RNG, 0.0, 1.0));
  
  return (u < lr);
}

/* MCMC outpath file name */
void file_outpath(char *path, size_t path_sz, int chain_id) {
  
  snprintf(path, path_sz, "output/mcmc.%d.txt", chain_id + 1);
  
}

/* I just wanna see how long things take */
double what_time(void){
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Run MCMC chain */
void run_chain(int chain_id, int Nsim, int burnin, Individuals *indiv, 
               Global global, Incidence incidence, Parameters param,
               double init_alpha, double init_beta, double init_gamma, 
               double init_r, double init_Gta, double init_Gtb){
  
  // file names
  char address_output[100];
  file_outpath(address_output, sizeof address_output, chain_id);
  
  // write out
  FILE *mcmc;
  mcmc = fopen(address_output, "w");
  fprintf(mcmc, "%s %s %s %s %s %s %s %s %s\n", "nSim", "alpha", "beta", "gamma", "r", "Gta", "Gtb", "logL", "accepted");
  
  printf("\rIn progress 0/%d", Nsim); fflush(stdout);
  
  // Initialization
  long double current_alpha = init_alpha;
  long double current_beta = init_beta;
  long double current_gamma = init_gamma;
  long double current_r = init_r;
  long double current_Gta = init_Gta;
  long double current_Gtb = init_Gtb;
  long double currentLogLike = -999999.9;
  int accepted = 0;
  initInfecData(indiv);
  
  printf("\rChain %d – In progress 0/%d", chain_id + 1, Nsim);
  fflush(stdout);
  
  // MCMC iterations
  for (int it = 0; it < Nsim; it++) {
    double start_time = what_time();
    
    long double oldLogLike = currentLogLike;
    
    // update infection times
    GenDensity(&param, global); // for likelihood and FOI calculation
    updateInfecData(indiv, global, incidence, param);
    
    double after_infec = what_time();
    
    // go through parameters we want to update
    long double *current[] = {&current_alpha, &current_beta, &current_gamma, 
                              &current_r, &current_Gta, &current_Gtb};
    double *sigma[] = {&param.sigma_alpha, &param.sigma_beta, &param.sigma_gamma,
                       &param.sigma_r, &param.sigma_Gta, &param.sigma_Gtb};;
    // int *accept[] = {&accepted_alpha, &accepted_beta, &accepted_gamma,
    //                  &accepted_r, &accepted_Gta, &accepted_Gtb};
    void which_param(int idx, double value){
      switch (idx) {
      case 0: param.alpha = value; break;
      case 1: param.beta = value; break;
      case 2: param.gamma = value; break;
      case 3: param.r = value; break;
      case 4: param.Gta = value; break;
      case 5: param.Gtb = value; break;
      }
    }
    
    for (int p = 0; p<6; p++){
      double oldVal = *current[p];
      double proposed = rnorm_reflect_pos(oldVal, *sigma[p]);
      
      which_param(p, proposed);
      GenDensity(&param, global);
      long double candidateLogLike = LogLikelihood(indiv, global, incidence, param);
      
      if (stageAccept(candidateLogLike, oldLogLike)) {
        *current[p] = proposed;
        oldLogLike = currentLogLike;
        currentLogLike = candidateLogLike;
        //(*accept[p])++;
        accepted++;
      }
    }
    
    // Write to file after burn-in 
    if (it >= burnin) {
      fprintf(mcmc, "%d %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %.10Lf %d\n",
              it+1, current_alpha, current_beta, current_gamma, current_r,
              current_Gta, current_Gtb, currentLogLike, accepted);
      fflush(mcmc);
    } 
    
    printf("\rChain %d – In progress %d/%d", chain_id + 1, it + 1, Nsim);
    fflush(stdout);
    
    double end_time = what_time();
    double elapsed = end_time - start_time;
    printf("\tit %d: %.9f seconds", it, elapsed); fflush(stdout);
    printf("\tinfec times: %.9f seconds", after_infec - start_time); fflush(stdout);
    
  } // end MCMC iterations 
  
  fclose(mcmc);
  printf("\n");
}