from cffi import FFI
import os

ffi = FFI()
LOCATION = os.path.dirname(os.path.abspath(__file__))
CLOC = LOCATION + '/src/py21cmmc/_21cmfast'

# This is the overall C library
ffi.set_source(
    "py21cmmc._wrapped_21cmfast",
    '''
    #include "drive_21cmMC_streamlined.c"
    
    // THIS IS A PART OF DRIVE_21cmMC_STREAMLINED, BUT NEEDS TO BE CALLED IN C BECAUSE OF THE FFTWF_MALLOC
    void init_inhomo_reco(){
        int ct, i, j, k;
        
        z_re = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS); // the redshift at which the cell is ionized
        Gamma12 = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);  // stores the ionizing backgroud
        N_rec_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // cumulative number of recombinations
        N_rec_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++) {
            z_re[ct] = -1.0;
        }

        // initialize N_rec
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k)) = 0.0;
                }
            }
        }
    }
    
    void destroy_inhomo_reco(){
        fftwf_free(N_rec_unfiltered);
        fftwf_free(N_rec_filtered);
        fftwf_free(z_re);
        fftwf_free(Gamma12);
    }
    ''',
    include_dirs = [CLOC],
    libraries=['m','gsl','gslcblas','fftw3f_omp', 'fftw3f'],
    #extra_compile_args = ['-fopenmp', '-Ofast', '-w']
    extra_compile_args = ['-fopenmp',  '-w', '-g', '-O0'] # For debugging only
)


# Now define the stuff we want access to (CAN'T DO #defines and #includes here!)
# First all the global variables
unallowed = [
    "#",
    "fftwf",
]
# for fl in ['Variables.h', "INIT_PARAMS.H", "ANAL_PARAMS.H"]:
#     with open(CLOC+"/Parameter_files/%s"%fl) as f:
#         lines = [ln for ln in f.readlines() if not any([ln.startswith(u) for u in unallowed])]
#         for i, ln in enumerate(lines):
#
#             if ln.startswith("#define"):
#                 lnst = ln.split(" ")
#                 if len(lnst)>2:
#                     if "(" not in lnst[1] and ")" not in lnst[1]:
#                         lines[i] = " ".join(lnst[:2])+ " ...\n"
#                     else:
#                         lines[i] = '\n'
#                 else:
#                     lines[i] = "\n"
#         with open("/home/steven/DERP%s"%fl,'w') as f:
#             f.write("".join(lines))
#         ffi.cdef("".join(lines))

for fl in ['Variables.h']:
    with open(CLOC+"/Parameter_files/%s"%fl) as f:
        lines = [ln for ln in f.readlines() if not any([ln.startswith(u) for u in unallowed])]
        ffi.cdef("".join(lines))
# Now the signatures of the various functions that we care about
ffi.cdef(
    """
    // Extra "defines" constants from INIT and ANAL.
    static const int MIDDLE;
    static const unsigned long long D;
    static const unsigned long long MID;
    static const unsigned long long TOT_NUM_PIXELS;
    static const unsigned long long TOT_FFT_NUM_PIXELS;
    static const unsigned long long KSPACE_NUM_PIXELS;

    static const unsigned long long HII_D;
    static const int HII_MIDDLE;
    static const unsigned long long HII_MID;
    
    static const unsigned long long HII_TOT_NUM_PIXELS;
    static const unsigned long long HII_TOT_FFT_NUM_PIXELS;
    static const unsigned long long HII_KSPACE_NUM_PIXELS;
    
    static const float LC_BOX_PADDING_IN_MPC;
    
    static const float DELTA_K;
    
    
    float REDSHIFT; // This should really be in Variables.h actually. 
    
    static const double CMperMPC;
    static const double NU_over_EV;
    
    double drdz(float z);
    double splined_erfc(double x);
    void init_MHR();
    void init_inhomo_reco();
    void destroy_inhomo_reco();
    
    // Stuff from driver
    struct ReturnData drive_21CMMC(char* arg1, char* arg2, struct BoxDimStruct BoxDim, struct FlagOptionStruct FlagOptions,
                               struct AstroParamStruct AstrophysicalParams, struct CosmoParamStruct CosmologicalParams);
                               
    
    void init_21cmMC_Ts_arrays();
    void init_21cmMC_HII_arrays();
    
    void ComputeTsBoxes();
    void ComputeIonisationBoxes(int sample_index, float REDSHIFT_SAMPLE, float PREV_REDSHIFT, int ICsInMemory);
    
    void ComputeInitialConditions();
    void ComputePerturbField(float REDSHIFT_SAMPLE);
    
    void GeneratePS(int CO_EVAL, double AverageTb);
    void ReadFcollTable();
    
    void destroy_21cmMC_Ts_arrays();
    void destroy_21cmMC_HII_arrays(int skip_deallocate);
    
    // Stuff from ps.c
    void init_ps();
    void free_ps();
    
    double tau_e(float zstart, float zend, float *zarry, float *xHarry, int len);
    """
)

if __name__=="__main__":
    ffi.compile()
