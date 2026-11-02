#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <time.h>
#include "headerSchools.h"

gsl_rng *RNG;   // adding this here bc i was getting errors otherwise?

int main(int argc, char *argv[]){
    
    
    // Command line and basics ---------------------------------------------------------------------------------------------------//

    int NinfSets, NparamSets;
    char *seed, *NinfSets_, *NparamSets_; 
    char *address_glo, *address_ind, *address_neg, *address_schchar, *address_idxFam, *address_idxSch, *address_idxHNegSch, *address_idxCom, *address_idxComNeg, *address_out, *address_inf, *address_param;
    char help[] = "\nERROR COMMAND LINE\nusage: ./execute/executable -sd [seed] -ag [address_glo] -ai [address_ind] -an [address_neg] -axf [address_idxFam] -axs [address_idxSch] -axc [address_indexCom] -axcn [address_indexComNeg] -ao [address_outIncidence] -ainf [address_infectiontimes] -asc [address_schchar] -axns [address_idxHNegSch] -apar [address_parameters] -ninft [Ninft] -np [Np]\n";
    // -sd: seed; -ag: input global data; -ai: input individual data; -an: input negative household individuals' data
	// -axf: daily index of attendance to the household (for positive cases); -axs: daily index of attendance to the school (for positive cases); -axc: daily index of attendance to the community (for positive cases); -axcn: daily index of attendance to the community (for the negative household members)
	// -ao: daily contribution of community cases; -ainf: input infection times; -asc: input school characteristics; -axns: daily index of attendance to the school (for the negative school contacts)
	// -apar: selection of posterior parameters; -ninft: = np; -np: number of posterior parameters selected

    if(argc < 31){fprintf(stderr,"%s\n",help); exit(1);}

    seed = argv[2]; if(seed == NULL){seed = "0";}
    address_glo = argv[4]; if(address_glo == NULL){address_glo = "0";}
    address_ind = argv[6]; if(address_ind == NULL){address_ind = "0";}
    address_neg = argv[8]; if(address_neg == NULL){address_neg = "0";}
    address_idxFam = argv[10]; if(address_idxFam == NULL){address_idxFam = "0";}
    address_idxSch = argv[12]; if(address_idxSch == NULL){address_idxSch = "0";}
    address_idxCom = argv[14]; if(address_idxCom == NULL){address_idxCom = "0";}
    address_idxComNeg = argv[16]; if(address_idxComNeg == NULL){address_idxComNeg = "0";}
    address_out = argv[18]; if(address_out == NULL){address_out = "0";}
    address_inf = argv[20]; if(address_inf == NULL){address_inf = "0";}
    address_schchar = argv[22]; if(address_schchar == NULL){address_schchar = "0";}
    address_idxHNegSch = argv[24]; if(address_idxHNegSch == NULL){address_idxHNegSch = "0";}
    address_param = argv[26]; if(address_param == NULL){address_param = "0";}
    NinfSets_ = argv[28]; if(NinfSets_ == NULL){NinfSets_ = "0";}
    NinfSets = atoi(NinfSets_);
    NparamSets_ = argv[30]; if(NparamSets_ == NULL){NparamSets_ = "0";}
    NparamSets = atoi(NparamSets_);
    
    
    // Output address:
    char variant[10] = {};
    memcpy(variant,&address_ind[18],5);
    char address_output[100] = {};
    strcat(address_output, "output/chain.");
    strcat(address_output, variant);
    strcat(address_output, "-s1.txt");
    
    // Seed:
    setenv("GSL_RNG_SEED", seed, 1);
    gsl_rng_env_setup();
    RNG = gsl_rng_alloc(gsl_rng_default);


    // Variables: ----------------------------------------------------------------------------------------------------------------//
   
    int i, it, infector, place;
    Indiv indiv;
    Negs negs;
    Out out;
    Global global;
    SchChar schchar;
    Param param;
    Chain chain;
    int **infecTimes; 
    long double **setParams;  


    // Allocation and Reading Files: ---------------------------------------------------------------------------------------------//

    ReadFiles1(address_glo, &indiv, &negs, &global, &param);

    param.genDensity = (double *)calloc(global.MAX_GENTIME+1, sizeof(double));

    schchar.IDcluster = (int *)calloc(global.Nclusters, sizeof(int));
    schchar.typeSch = (int *)calloc(global.Nclusters, sizeof(int));
    schchar.Nexposures = (int *)calloc(global.Nclusters, sizeof(int));
    schchar.Npositives = (int *)calloc(global.Nclusters, sizeof(int));
    schchar.Nattending = (int *)calloc(global.Nclusters, sizeof(int));
    schchar.idxSch = (int **)calloc(global.Nclusters, sizeof(int *));
    for(i=0; i<global.Nclusters; i++){
        (schchar.idxSch)[i] = (int *)calloc(global.NDAYS, sizeof(int));
    }

    indiv.IDperson = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.IDfamily = (int *)calloc(indiv.Nindiv, sizeof(int));
    indiv.IDcluster = (int *)calloc(indiv.Nindiv, sizeof(int));
    indiv.symonset = (int *)calloc(indiv.Nindiv, sizeof(int));
    indiv.diagnose = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg1 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg2 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg3 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg4 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg5 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg6 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg7 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg8 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.neg9 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.vax1 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.vax2 = (int *)calloc(indiv.Nindiv, sizeof(int)); 
    indiv.vax3 = (int *)calloc(indiv.Nindiv, sizeof(int));
    indiv.attend = (int *)calloc(indiv.Nindiv, sizeof(int));
    
    indiv.idxCom = (int **)calloc(indiv.Nindiv, sizeof(int *));
    for(i=0; i<indiv.Nindiv; i++){
        (indiv.idxCom)[i] = (int *)calloc(global.NDAYS, sizeof(int));
    }
    indiv.idxFam = (int **)calloc(indiv.Nindiv, sizeof(int *));
    for(i=0; i<indiv.Nindiv; i++){
        (indiv.idxFam)[i] = (int *)calloc(global.NDAYS, sizeof(int));
    }
    indiv.idxSch = (int **)calloc(indiv.Nindiv, sizeof(int *));
    for(i=0; i<indiv.Nindiv; i++){
        (indiv.idxSch)[i] = (int *)calloc(global.NDAYS, sizeof(int));
    }
    indiv.xi = (long double **)calloc(indiv.Nindiv, sizeof(long double *));
    for(i=0; i<indiv.Nindiv; i++){
        (indiv.xi)[i] = (long double *)calloc(global.NDAYS, sizeof(long double));
    }
    indiv.rho = (double **)calloc(indiv.Nindiv, sizeof(double *));
    for(i=0; i<indiv.Nindiv; i++){
        (indiv.rho)[i] = (double *)calloc(global.NDAYS, sizeof(double));
    }

    negs.IDperson = (int *)calloc(negs.Nindiv, sizeof(int)); 
    negs.IDfamily = (int *)calloc(negs.Nindiv, sizeof(int));
    negs.vax1 = (int *)calloc(negs.Nindiv, sizeof(int)); 
    negs.vax2 = (int *)calloc(negs.Nindiv, sizeof(int)); 
    negs.vax3 = (int *)calloc(negs.Nindiv, sizeof(int));

    negs.idxCom = (int **)calloc(negs.Nindiv, sizeof(int *));
    for(i=0; i<negs.Nindiv; i++){
        (negs.idxCom)[i] = (int *)calloc(global.NDAYS, sizeof(int));
    }
    negs.xi = (long double **)calloc(negs.Nindiv, sizeof(long double *));
    for(i=0; i<negs.Nindiv; i++){
        (negs.xi)[i] = (long double *)calloc(global.NDAYS, sizeof(long double));
    }

    out.date = (int *)calloc(global.NDAYS, sizeof(int)); 
    out.incidence = (double *)calloc(global.NDAYS, sizeof(double));

    infecTimes = (int **)calloc(NinfSets, sizeof(int *));
    for(i=0; i<NinfSets; i++){
        infecTimes[i] = (int *)calloc(indiv.Nindiv, sizeof(int));
    }

    chain.probInf = (long double *)calloc(indiv.Nindiv, sizeof(long double));
    chain.infector = (int *)calloc(indiv.Nindiv, sizeof(int));
    chain.placeInf = (int *)calloc(indiv.Nindiv, sizeof(int));

    ReadFiles2(address_ind, address_neg, address_schchar, address_idxFam, address_idxSch, address_idxCom, 
               address_idxComNeg, address_idxHNegSch, address_out, address_inf, &indiv, &negs, &schchar, 
               global, &out, &infecTimes, NinfSets);

    Xi(&indiv, global);
    Rho(&indiv, global); 
    XiNeg(&negs, global);

    setParams = (long double **)calloc(NparamSets, sizeof(long double *));
    for(i=0; i<NparamSets; i++){
        setParams[i] = (long double *)calloc(4, sizeof(long double));
    }

    // Reading here the file with the parameters sets
    i = 0;
    FILE *paramFile;
    paramFile = fopen(address_param, "r");
    if(paramFile == NULL){
        fprintf(stderr, "\nERROR opening the file %s for reading\n", address_param);
        exit(2);
    }
    while(!feof(paramFile) && (i < NparamSets)){
        fscanf(paramFile, "%Lf %Lf %Lf %Lf",
               &(setParams[i][0]), // alpha
               &(setParams[i][1]), // gamma
               &(setParams[i][2]), // Gta
               &(setParams[i][3])); // Gtb
        i++;
    }
    fclose(paramFile);


    // Computations: -------------------------------------------------------------------------------------------------------------//

    FILE *chainFile;
    chainFile = fopen(address_output,"w");
    fprintf(chainFile, "%s %s %s %s %s\n", "idxIndiv", "IDindiv", "idxInfector", "IDinfector", "placeInf");
    
    for(it=0; it<NparamSets; it++){

        param.alpha = setParams[it][0];
        param.beta = 0;
        param.gamma = setParams[it][1];
        param.Gta = setParams[it][2];
        param.Gtb = setParams[it][3];

        GenDensity(&param, global);
        ChainRec(indiv, global, out, param, &chain, infecTimes, it);

        for(i=0; i<indiv.Nindiv; i++){

            infector = chain.infector[i];
            place = chain.placeInf[i];
        
            if(infector == -1){
                fprintf(chainFile, "%d %d %d %d %d\n", i, indiv.IDperson[i], infector, -1, place);
            } 
            else if(infector == -999){
                fprintf(chainFile, "%d %d %d %d %d\n", i, indiv.IDperson[i], infector, -999, place);
            }
            else{
                fprintf(chainFile, "%d %d %d %d %d\n", i, indiv.IDperson[i], infector, indiv.IDperson[(int)(infector)], place);
            }
        }

        printf("\rIn progress %d/%d", it+1, NparamSets);
        fflush(stdout);

    }

    
    fclose(chainFile); 
    printf("\n");
    

    // Free: ---------------------------------------------------------------------------------------------------------------------//
    
    free(param.genDensity);

    free(schchar.IDcluster);
    free(schchar.typeSch);
    free(schchar.Nexposures);
    free(schchar.Npositives);
    free(schchar.Nattending);
    for(i=0; i<global.Nclusters; i++){
        free(schchar.idxSch[i]);
    }
    free(schchar.idxSch);

    free(indiv.IDperson);
    free(indiv.IDfamily);
    free(indiv.IDcluster);
    free(indiv.symonset);
    free(indiv.diagnose);
    free(indiv.neg1);
    free(indiv.neg2);
    free(indiv.neg3);
    free(indiv.neg4);
    free(indiv.neg5);
    free(indiv.neg6);
    free(indiv.neg7);
    free(indiv.neg8); 
    free(indiv.neg9);
    free(indiv.vax1);
    free(indiv.vax2);
    free(indiv.vax3);
    free(indiv.attend); 
    for(i=0; i<indiv.Nindiv; i++){
        free(indiv.idxCom[i]);
        free(indiv.idxFam[i]);
        free(indiv.idxSch[i]);
        free(indiv.xi[i]);
        free(indiv.rho[i]);
    }
    free(indiv.idxCom); 
    free(indiv.idxFam); 
    free(indiv.idxSch); 
    free(indiv.xi); 
    free(indiv.rho);

    free(negs.IDperson);
    free(negs.IDfamily);
    free(negs.vax1);
    free(negs.vax2);
    free(negs.vax3);
    for(i=0; i<negs.Nindiv; i++){
        free(negs.idxCom[i]);
        free(negs.xi[i]);
    }
    free(negs.idxCom);
    free(negs.xi);

    free(out.date);
    free(out.incidence);

    for(i=0; i<NinfSets; i++){
        free(infecTimes[i]);
    }
    free(infecTimes);

    free(chain.probInf);
    free(chain.infector);
    free(chain.placeInf);

    for(i=0; i<NparamSets; i++){
        free(setParams[i]);
    }
    free(setParams);
   

    return 0;
}




