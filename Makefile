EXEC   = 2LPTic

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o compute_ffts.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile

FFT_TYPE = FFTW3
FORMAT = GADGET
GAS = NO

#OPT += -DWDM_GAUSSIAN_VELOCITIES

ifeq ($(FFT_TYPE),FFTW3)
OPT += -DUSE_FFTW3
else ifeq ($(FFT_TYPE),PFFT)
OPT += -DUSE_PFFT
else
OPT += -DSINGLEPRECISION_FFTW2
endif

ifeq ($(FORMAT),SWIFT)
OPT += -DWRITE_MASSES
OPT += -DSWIFT_ICS
OPT += -DSINGLE_HDF5
endif

ifeq ($(GAS),YES)
OPT += -DPRODUCEGAS
endif

#OPT += -DSINGLE_HDF5

OPT += -DUSE_CAMB    # Allow input tabulated transfer function to be in
                     # CAMB format

#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing

#OPT  += -DWRITE_MASSES  # write mass block

#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type

#OPT   +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles start from a glass (as opposed to grid)

#OPT += -DONLY_ZA # swith this on if you want ZA initial conditions (2LPT otherwise)

OPT += -DLONGIDS

OPTIONS =  $(OPT)

OPTIONS =  $(OPT)
CC      =   mpicc -m64 #-std=c11
OPTIMIZE =  -O3

FFTW2_BASE=${HOME}/Codes/
GSL_BASE=${GSL_ROOT}

ifeq (SINGLE_HDF5,$(findstring SINGLE_HDF5,$(OPT)))
HDF_BASE=${HDF5_BASE}
HDF_EXTRAS=#${HDF5_BASE}/lib/libhdf5_hl.a ${HDF5_BASE}/lib/libhdf5.a
else
HDF_BASE=/usr/local
HDF_EXTAS=
endif

FFTW_INCL = -I$(FFTW2_BASE)/include
FFTW_LIBS = -L$(FFTW2_BASE)/lib

GSL_LIBS =   -L$(GSL_BASE)/lib
GSL_INCL =  -I$(GSL_BASE)/include

HDF_LIBS = -L$(HDF_BASE)/lib $(HDF_EXTRAS) -lz -ldl -lm
HDF_INCL = -I$(HDF_BASE)/include #-DH5_USE_18_API

ifeq (USE_FFTW3,$(findstring USE_FFTW3,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -lfftw3_mpi -lfftw3 -lm
else ifeq (USE_PFFT,$(findstring USE_PFFT,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
else
ifeq (SINGLEPRECISION_FFTW2,$(findstring SINGLEPRECISION_FFTW2,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
endif
endif

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS) $(HDF_LIBS) -lgsl -lgslcblas -lhdf5_hl -lhdf5

#Note that the order in which the includes are set is important for ensuring
#that the correct HDF5 header is picked up
CFLAGS =   $(OPTIONS)  $(OPTIMIZE) $(HDF_INCL) $(FFTW_INCL) $(GSL_INCL) 

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



