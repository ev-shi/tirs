#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

extern gsl_rng *RNG;

// Constants: -------------------------------------------------------------------------------------------------
static const double steps = 500.0;

// static const double init_alpha = 0.002;
// static const double init_beta = 0.0002;
// static const double init_gamma = 0.5;
// static const double init_r = 0.05;
// static const double init_Gta = 5;
// static const double init_Gtb = 1;

// Define structures of interest 
typedef struct{
  int Nindiv;
  int *IDperson;   // ID of individual
  int *IDblock;    // ID of block
  int *IDcluster;  // ID of cluster
  int *TIRS;       // 0 control arm, 1 treatment arm
  int *symptonset; // date of symptom onset
  double *infecTime; // inferred date of infection
} Individuals;

typedef struct{
  int *date; // week of Jan 1st is Week 0, week of Dec 31 is Week 52
  double *incidence; // Incidence of each week of the year in the city of Mérida
  int POBTOT; // total population size
} Incidence;

typedef struct{
  long double alpha; // cumulative hazard in block
  long double beta;  // cumulative hazard in the cluster
  long double gamma; // transmissibility in Mérida
  long double r;     // efficacy of TIRS intervention
  long double Gta;    // shape parameter of generation time Gamma distribution
  long double Gtb;    // scale parameter of generation time Gamma distribution
  double sigma_alpha;  // how much to change alpha by
  double sigma_beta;
  double sigma_gamma;
  double sigma_r; 
  double sigma_Gta;
  double sigma_Gtb;
  double *genDensity; // vector with generation time density evaluated from 0 to MAX_GENTIME
} Parameters;

typedef struct{
  int Nclusters; // number of clusters in the trial
  int Nweeks;    // 52
  int Week0;     // 0
  int WeekL;     // 51
  int MAX_GENTIME; // maximum generation time allowed
  int MINT; // initial time of simulations -- min_infection_time - MAX_GENTIME
  int MAXT; // last time of simulations -- last symptom onset
} Global;

// typedef struct{
//   long double *probInf; // probability of infection of individual j
//   int *infector;        // Index i of infector of j, -1 if community infection
//   int *placeInf;        // location of infection, 0 block, 1 cluster, 2 city
// } Chain;

// Functions: -------------------------------------------------------------------------------------------------

/* Reads global variables, the number of individuals, and the sigmas for the MCMC parameters */
void ReadFiles1(char *address_global, Individuals *indiv, Global *global, Incidence *incidence, Parameters *param);

/* Reads individuals data and city incidence */
void ReadFiles2(char *address_indiv, char *address_incidence, Individuals *indiv,
                Global global, Incidence *incidence);

/* Initialize complete data for model with observed symptom onset times.*/
void initInfecData(Individuals *indiv);

/* Calculate instantaneous FOI for given individual j and week t*/
long double foi_at_t(int j, int t, Individuals *indiv, Global global, Incidence incidence, Parameters param);

/* Calculate Likelihood */
long double LogLikelihood(Individuals *indiv, Global global, Incidence incidence, Parameters param);

/* Update unknown infection times using Metropolis-Hastings step*/
void updateInfecData(Individuals *indiv, Global global, Incidence incidence, Parameters param);

/* Draw from a reflected normal distribution about 0 */
double rnorm_reflect_pos(double mean, double sd);

/* Compute vector for generation time density, evaluated from 0 to MAX_GENTIME with given parameters */
void GenDensity(Parameters *param, Global global);

/* Accept candidate? */
int stageAccept(long double candidate, long double old);

/* MCMC outpath file name */
void file_outpath(char *path, size_t path_sz, int chain_id);

/* I just wanna see how long things take */
double what_time(void);

/* Run MCMC chain */
void run_chain(int chain_id, int Nsim, int burnin, Individuals *indiv, 
               Global global, Incidence incidence, Parameters param,
               double init_alpha, double init_beta, double init_gamma, 
               double init_r, double init_Gta, double init_Gtb);

