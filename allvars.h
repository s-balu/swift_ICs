#ifdef USE_FFTW3
#include <fftw3-mpi.h>
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
#elif USE_PFFT
#else
#ifdef SINGLEPRECISION_FFTW2
#include <srfftw_mpi.h>
#else
#include <rfftw_mpi.h>
#endif
#endif


#define ASSERT_ALLOC(cond) {                                                                                  \
   if(cond)                                                                                                   \
    {                                                                                                         \
      if(ThisTask == 0)                                                                                       \
	printf("\nallocated %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                     \
    }                                                                                                         \
  else                                                                                                        \
    {                                                                                                         \
      printf("failed to allocate %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                \
      printf("bailing out.\n");                                                                               \
      FatalError(1);                                                                                          \
    }                                                                                                         \
}

extern int OutputFormat;

#define  PI          3.14159265358979323846 
#define  GRAVITY     6.672e-8
#define  HUBBLE      3.2407789e-18   /* in h/sec */

double PowerSpec(double kmag);
double GrowthFactor(double astart, double aend);
double F_Omega(double a);
double F2_Omega(double a);
int    read_parameter_file(char *fname);
double PowerSpec_EH(double k);
double PowerSpec_Efstathiou(double k);

#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
typedef unsigned short int uint4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
typedef unsigned int uint4byte;
typedef long long int64byte;
typedef unsigned long long uint64byte;
#endif

#define N_GADGET_TYPE 6
extern struct io_header_1
{
#ifdef NO64BITID
  int npart[N_GADGET_TYPE];
#else
  unsigned long long npart[N_GADGET_TYPE];  
#endif
  double       mass[N_GADGET_TYPE];
  double       time;
  double       redshift;
  int4byte     flag_sfr;
  int4byte     flag_feedback;
#ifdef NO64BITID
  int npartTotal[N_GADGET_TYPE];
#else
  unsigned long long npartTotal[N_GADGET_TYPE];
#endif
  int          flag_cooling;
  int          num_files;
  double       BoxSize;
  double       Omega0;
  double       OmegaLambda;
  double       HubbleParam;
  int4byte     flag_stellarage;
  int4byte     flag_metals;
#ifdef NO64BITID
  unsigned int    n_all_high[N_GADGET_TYPE];
#else
  unsigned long long   n_all_high[N_GADGET_TYPE];  
#endif
  int4byte     flag_entropyICs;
  int4byte     flag_doubleprecision;
  char         unused[52];
}
  header, header1;


extern int      Nglass;
extern int      *Local_nx_table;
extern int      WhichSpectrum;


extern FILE     *FdTmp, *FdTmpInput;

extern int      Nmesh, Nsample;

extern int      SphereMode;

extern long long IDStart;


extern char     GlassFile[500]; 
extern char     FileWithInputSpectrum[500];

extern int      GlassTileFac; 

extern double   Box;
extern int Seed;

extern long long TotNumPart;

extern int      NumPart;
#ifdef SINGLE_HDF5
extern int      MaxNumPart;
#endif

extern int      NTaskWithN;

extern int      NumPartTask;


extern int      *Slab_to_task;


extern struct part_data 
{
  float Pos[3];
  float Vel[3];
#ifdef  MULTICOMPONENTGLASSFILE                      
  int   Type;
#endif
  long long ID;
} *P;


extern double InitTime;
extern double Redshift;
extern double MassTable[6];


extern char OutputDir[100], FileBase[100];
extern int  NumFilesWrittenInParallel;


extern int      ThisTask, NTask;

extern int      Local_nx, Local_x_start;

extern int  IdStart;

extern unsigned int TotalSizePlusAdditional;

#ifdef USE_FFTW3
extern fftw_plan Inverse_plan;
extern fftw_plan Forward_plan;
#elif USE_PFFT
#else
extern rfftwnd_mpi_plan Inverse_plan;
extern rfftwnd_mpi_plan Forward_plan;
extern fftw_real        *Workspace;
#endif


extern double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
extern double InputSpectrum_UnitLength_in_cm;
extern double G, Hubble;
extern double RhoCrit;

extern double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
extern double OmegaBaryon, HubbleParam;
extern double PrimordialIndex;
extern double ShapeGamma;

extern double Dplus; /* growth factor */


#ifdef DIFFERENT_TRANSFER_FUNC
extern int Type, MinType, MaxType;
#endif

extern int    WDM_On;
extern int    WDM_Vtherm_On;
extern double WDM_PartMass_in_kev;
