#=======================================================================
#  Compiler options
#=======================================================================
# FC = h5pfc -fpp # ifort preprocessor
FC = h5pfc -cpp # gfortran preprocessor
## Laptop
# FC += -r8 -O3 # ifort options
FC += -fdefault-real-8 -fdefault-double-8 -O2 # gfortran options
# FC += -fdefault-real-8 -fdefault-double-8 -O0 -g -fbacktrace -fbounds-check
## Cartesius
# FC += -r8 -O3 -xAVX -axCORE-AVX2
## Irene
# FC += -r8 -O3 -mavx2 $(FFTW3_FFLAGS)
## MareNostrum
# FC += -r8 -O3 $(FFTW_FFLAGS)
## Traceback / Debug
# FC += -r8 -O0 -g -traceback -check bounds 
# FC += -DSHM -DSHM_DEBUG

#=======================================================================
#  Library
#======================================================================
# Common build flags
##  Laptop (with intel libraries)
# LDFLAGS = -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lhdf5_fortran -lhdf5 -lsz -lz -ldl -lm
## Laptop (with GNU libraries)
LDFLAGS = -lfftw3 -llapack -lblas -lpthread -lhdf5_fortran -lhdf5 -lsz -lz -ldl -lm

## Cartesius (Before compiling, load modules 2019 intel/2018b HDF5 FFTW)
# FFTW3_LIBS = -lfftw3
# BLAS_LIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
# HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz -ldl -lm
# LDFLAGS = $(FFTW3_LIBS) $(BLAS_LIBS) $(HDF5_LIBS)

## Irene (Before compiling, load modules flavor/hdf5/parallel hdf5 fftw3/gnu)
# FFTW3_LIBS = $(FFTW3_LDFLAGS)
# BLAS_LIBS = $(MKL_LDFLAGS)
# HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz -ldl -lm
# LDFLAGS = $(FFTW3_LIBS) $(BLAS_LIBS) $(HDF5_LIBS)

## MareNostrum (Before compiling, load modules fftw hdf5)
# FFTW3_LIBS = $(FFTW_LIBS)
# BLAS_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
# HDF5_LIBS = -lhdf5_fortran -lhdf5  -lsz -lz -ldl -lm
# LDFLAGS = $(FFTW3_LIBS) $(BLAS_LIBS) $(HDF5_LIBS)

## SuperMUC-NG (Before compiling, load modules fftw hdf5 szip)
# BLAS_LIBS = $(MKL_LIB)
# FFTW3_LIBS = $(FFTW_LIB)
# HDF5_LIBS = $(HDF5_F90_SHLIB) $(HDF5_SHLIB) -L$(SZIP_LIBDIR) -lz -ldl -lm
# LDFLAGS = $(BLAS_LIBS) $(FFTW3_LIBS) $(HDF5_LIBS)

#=======================================================================
#  Non-module Fortran files to be compiled:
#=======================================================================
EXTRA_DIST = transpose_z_to_x.F90 transpose_x_to_z.F90 transpose_x_to_y.F90\
	     transpose_y_to_x.F90 transpose_y_to_z.F90 transpose_z_to_y.F90\
	     factor.F90 halo.F90 fft_common.F90 alloc.F90 halo_common.F90

# Object files associated with standard flow solver
OBJS = obj/main.o obj/CalcLocalDivergence.o obj/CalcMaxCFL.o \
	obj/CalcMeanProfiles.o obj/CheckDivergence.o obj/CorrectPressure.o \
	obj/CorrectVelocity.o obj/CreateGrid.o obj/CreateInitialConditions.o \
	obj/DeallocateVariables.o obj/DebugRoutines.o obj/ExplicitTermsTemp.o \
	obj/ExplicitTermsVX.o obj/ExplicitTermsVY.o obj/ExplicitTermsVZ.o \
	obj/factorize.o obj/HdfReadContinua.o obj/HdfRoutines.o \
	obj/ImplicitAndUpdateTemp.o obj/ImplicitAndUpdateVX.o obj/ImplicitAndUpdateVY.o \
	obj/ImplicitAndUpdateVZ.o obj/InitPressureSolver.o obj/InitTimeMarchScheme.o \
	obj/InitVariables.o obj/LocateLargeDivergence.o obj/MakeMovieXCut.o \
	obj/MakeMovieYCut.o obj/MakeMovieZCut.o obj/MpiAuxRoutines.o \
	obj/QuitRoutine.o obj/ReadInputFile.o obj/ResetLogs.o \
	obj/SetTempBCs.o obj/SolveImpEqnUpdate_Temp.o obj/SolveImpEqnUpdate_X.o \
	obj/SolveImpEqnUpdate_YZ.o obj/SolvePressureCorrection.o obj/SpecRoutines.o \
	obj/TimeMarcher.o obj/WriteFlowField.o obj/WriteGridInfo.o \
	obj/CalcWriteQ.o obj/GlobalQuantities.o obj/ReadFlowInterp.o

# Object files associated with multiple resolution grids
OBJS += obj/CreateMgrdGrid.o obj/CreateMgrdStencil.o obj/InitMgrdVariables.o \
	obj/DeallocateMgrdVariables.o

# Object files associated with initial condition interpolation
OBJS += obj/CreateInputStencil.o obj/CreateOldGrid.o obj/CreateSalStencil.o \
	obj/InterpInputSal.o obj/InterpInputVel.o obj/InterpSalMgrd.o \
	obj/InterpVelMgrd.o obj/InitInputVars.o obj/DeallocateInputVars.o \
	obj/InterpInputPhi.o

# Object files associated with the salinity field
OBJS += obj/ExplicitTermsSal.o obj/ImplicitAndUpdateSal.o obj/SolveImpEqnUpdate_Sal.o \
	obj/UpdateScalarBCs.o obj/CreateICSal.o obj/InitSalVariables.o \
	obj/DeallocateSalVariables.o obj/SetSalBCs.o

# Object files associated with the phase-field method
OBJS += obj/AddLatentHeat.o obj/DeallocatePFVariables.o obj/ExplicitTermsPhi.o \
	obj/ImplicitAndUpdatePhi.o obj/InitPFVariables.o obj/InterpPhiMgrd.o \
	obj/InterpTempMgrd.o obj/SolveImpEqnUpdate_Phi.o obj/CreateICPF.o \
	obj/ImmersedBoundary.o

# Module object files
MOBJS = obj/param.o obj/decomp_2d.o obj/AuxiliaryRoutines.o obj/decomp_2d_fft.o

#=======================================================================
#  Files that create modules:
#=======================================================================
MFILES = param.F90 decomp_2d.F90 AuxiliaryRoutines.F90 decomp_2d_fft.F90

# Object and module directory:
OBJDIR=obj
OUTDIR=outputdir outputdir/stst3 outputdir/stst outputdir/flowmov outputdir/flowmov3d
# OBJS  := $(FFILES:%.F90=$(OBJDIR)/%.o)
# MOBJS := $(MFILES:%.F90=$(OBJDIR)/%.o)

# when using ifort compiler:
# FC += -module $(OBJDIR) 
# when using gfortran:
FC += -J $(OBJDIR) 

#============================================================================ 
#  make PROGRAM   
#============================================================================
PROGRAM = afid 

#Compiling 
all: objdir outdir $(PROGRAM) 
$(PROGRAM): $(MOBJS) $(OBJS) 
	$(FC) -o $@ $^ $(LDFLAGS) 

#============================================================================
#  Dependencies 
#============================================================================
$(OBJDIR)/param.o: source/flow_solver/param.F90
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/AuxiliaryRoutines.o: source/flow_solver/AuxiliaryRoutines.F90 
	$(FC) -c -o $@ $< $(LDFLAGS) 
$(OBJDIR)/decomp_2d.o: source/flow_solver/2decomp/decomp_2d.F90
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/decomp_2d_fft.o: source/flow_solver/2decomp/decomp_2d_fft.F90
	$(FC) -c -o $@ $< $(LDFLAGS) 
$(OBJDIR)/%.o: source/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/%.o: source/flow_solver/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/%.o: source/multires/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/%.o: source/multires/IC_interpolation/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/%.o: source/multires/phase-field/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/%.o: source/multires/salinity/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)

#============================================================================
#  Clean up 
#============================================================================
clean: 
	/bin/rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*genmod* $(OBJDIR)/*.o obj\

.PHONY: objdir outdir
objdir: $(OBJDIR) 
$(OBJDIR): 
	mkdir -p ${OBJDIR} 
#outdir: $(OUTDIR) 
# $(OUTDIR): 
# 	mkdir -p ${OUTDIR} 
