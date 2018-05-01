#include "Parameter_files/INIT_PARAMS.H"
#include "Parameter_files/ANAL_PARAMS.H"
#include "Parameter_files/Variables.h"

//int main(int argc, char ** argv){
float ReturnTau_e(char* arg1, char* arg2, int READ_FROM_FILE, int PRINT_FILES,
                  struct CosmoParamStruct CosmologicalParams,
                  float z_Hist[], float xHI_Hist[], int N_REDSHIFTS) {

    char filename[500];
    char dummy_string[500];
    FILE *F;

    int i;

    float taue;

    INDIVIDUAL_ID = atof(arg1);
    INDIVIDUAL_ID_2 = atof(arg2);


    /////////////////   Read in the cosmological parameter data     /////////////////

    if(READ_FROM_FILE) {
        double *PARAM_COSMOLOGY_VALS = calloc(TOTAL_COSMOLOGY_FILEPARAMS,sizeof(double));

        sprintf(filename,"WalkerCosmology_%1.6lf_%1.6lf.txt",INDIVIDUAL_ID,INDIVIDUAL_ID_2);
        F = fopen(filename,"rt");

        for(i=0;i<TOTAL_COSMOLOGY_FILEPARAMS;i++) {
            fscanf(F,"%s\t%lf\n",&dummy_string,&PARAM_COSMOLOGY_VALS[i]);
        }
        fclose(F);

        // Assign these values. Hard-coded, so order is important
        RANDOM_SEED = (unsigned long long)PARAM_COSMOLOGY_VALS[0];
        SIGMA8 = (float)PARAM_COSMOLOGY_VALS[1];
        hlittle = (float)PARAM_COSMOLOGY_VALS[2];
        OMm = (float)PARAM_COSMOLOGY_VALS[3];
        OMl = (float)PARAM_COSMOLOGY_VALS[4];
        OMb = (float)PARAM_COSMOLOGY_VALS[5];
        POWER_INDEX = (float)PARAM_COSMOLOGY_VALS[6]; //power law on the spectral index, ns
    }
    else {

        RANDOM_SEED = CosmologicalParams.RANDOM_SEED;
        SIGMA8 = CosmologicalParams.SIGMA8;
        hlittle = CosmologicalParams.hlittle;
        OMm = CosmologicalParams.OMm;
        OMl = CosmologicalParams.OMl;
        OMb = CosmologicalParams.OMb;
        POWER_INDEX = CosmologicalParams.POWER_INDEX;

    }

    float *Redshifts = calloc(N_REDSHIFTS,sizeof(float));
    float *xH = calloc(N_REDSHIFTS,sizeof(float));

    for(i=0;i<N_REDSHIFTS;i++) {
        Redshifts[i] = (float)(z_Hist[i]);
        xH[i] = (float)(xHI_Hist[i]);
    }

    taue = tau_e(0, Redshifts[N_REDSHIFTS-1], Redshifts, xH, N_REDSHIFTS);


    if(PRINT_FILES) {
        sprintf(filename, "Tau_e_%1.6lf_%1.6lf.txt",INDIVIDUAL_ID,INDIVIDUAL_ID_2);
        F=fopen(filename, "wt");
        fprintf(F, "%lf\n",taue);
        fclose(F);
    }

    free(Redshifts);
    free(xH);

    return taue;
}
