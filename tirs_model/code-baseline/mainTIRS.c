#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "headerTIRS.h"

gsl_rng *RNG = NULL;

int main(int argc, char *argv[]){
  
  /* Command line and basics ------------------------------------------------------------------------------------------- */
  int Nsim, burnin;
  char *Nsim_, *burnin_;
  char *address_global, *address_indiv, *address_incidence;
  char help[] = "\nCOMMAND LINE ERROR\nusage: ./execute/executable -sm [Nsim] -bn [BurnIn] -ag [address_global] -ai [address_indiv] -ac [address_incidence]\n";
  
  if(argc < 11){fprintf(stderr,"%s\n",help); exit(1);}
  
  Nsim_ = argv[2]; if(Nsim_ == NULL){Nsim_ = "0";}
  burnin_ = argv[4]; if(burnin_ == NULL){burnin_ = "0";}
  Nsim = atoi(Nsim_); // string to int
  burnin = atoi(burnin_);

  address_global = argv[6]; if(address_global == NULL){address_global = "0";}
  address_indiv = argv[8]; if(address_indiv == NULL){address_indiv = "0";}
  address_incidence = argv[10]; if(address_incidence == NULL){address_incidence = "0";}
  
  /* Variables ---------------------------------------------------------------------------------------------------------- */
  Individuals indiv;
  Incidence incidence;
  Global global;
  Parameters param;
  
  /* Allocation and Reading Files ----------------------------------------------------------------------------------------- */
  ReadFiles1(address_global, &indiv, &global, &incidence, &param);
  
  param.genDensity = (double *)calloc(global.MAX_GENTIME+1, sizeof(double));
  
  indiv.IDperson = (int *)calloc(indiv.Nindiv, sizeof(int)); 
  indiv.IDblock = (int *)calloc(indiv.Nindiv, sizeof(int));
  indiv.IDcluster = (int *)calloc(indiv.Nindiv, sizeof(int));
  indiv.TIRS = (int *)calloc(indiv.Nindiv, sizeof(int));
  indiv.symptonset = (int *)calloc(indiv.Nindiv, sizeof(int));
  indiv.infecTime = (double *)calloc(indiv.Nindiv, sizeof(double));
  
  incidence.date = (int *)calloc(global.Nweeks, sizeof(int)); 
  incidence.incidence = (double *)calloc(global.Nweeks, sizeof(double));
  
  ReadFiles2(address_indiv, address_incidence, &indiv, global, &incidence); 
  
  gsl_rng_env_setup();
  RNG = gsl_rng_alloc(gsl_rng_default);
  
  /* MCMC --------------------------------------------------------------------------------------------------------------- */
  
  // double starters[3][6] = {
  //   {0.0020, 0.000222, 2, 0.04, 5, 1.5},
  //   {0.0024, 0.000178, 5, 0.03, 3, 0.5},
  //   {0.0028, 0.000267, 0.5, 0.05, 1, 0.333}};
  // 
  // int chain_id = atoi(getenv("SLURM_ARRAY_TASK_ID"));
  // double init_alpha  = starters[chain_id][0];
  // double init_beta   = starters[chain_id][1];
  // double init_gamma  = starters[chain_id][2];
  // double init_r      = starters[chain_id][3];
  // double init_Gta    = starters[chain_id][4];
  // double init_Gtb    = starters[chain_id][5];

  // testing bc what is going on
  int chain_id = 0;
  double init_alpha  = 0.0024;
  double init_beta   = 0.000222;
  double init_gamma  = 2;
  double init_r      = 0.039;
  double init_Gta    = 4;
  double init_Gtb    = 0.5;
  
  
  // unique seed?
  gsl_rng_set(RNG, 987654321UL + chain_id);
  
  run_chain(chain_id, Nsim, burnin, &indiv, global, incidence, param, init_alpha,
            init_beta, init_gamma, init_r, init_Gta, init_Gtb);

  /* Free --------------------------------------------------------------------------------------------------------------- */
  gsl_rng_free(RNG);  
  free(param.genDensity);
  
  free(indiv.IDperson);
  free(indiv.IDblock);
  free(indiv.IDcluster);
  free(indiv.TIRS);
  free(indiv.symptonset);
  free(indiv.infecTime);
  
  free(incidence.date);
  free(incidence.incidence);
} // end main




