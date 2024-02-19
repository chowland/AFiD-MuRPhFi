# Choose the machine being used
# Options: PC, SNELLIUS, IRENE, MARENOSTRUM, SUPERMUC, DISCOVERER
MACHINE=PC
FLAVOUR=GNU
# Modules required for each HPC system as follows:
# SNELLIUS:
#	GNU: 2022 foss/2022a HDF5/1.12.2-gompi-2022a
# 	Intel: 2022 intel/2022a FFTW/3.3.10-GCC-11.3.0 HDF5/1.12.2-iimpi-2021a
# IRENE (Intel): flavor/hdf5/parallel hdf5 fftw3/gnu
# MARENOSTRUM (Intel): fabric intel mkl impi hdf5 fftw szip
# SUPERMUC (Intel): fftw hdf5
# DISCOVERER:
#	GNU: hdf5/1/1.14/latest-gcc-openmpi fftw/3/latest-gcc-openmpi lapack
#	Intel: hdf5/1/1.14/latest-intel-openmpi fftw/3/latest-gcc-openmpi mkl

#=======================================================================
#  Compiler options
#=======================================================================

# Object and module directory:
OBJDIR=obj

ifeq ($(FLAVOUR),GNU)
	FC = h5pfc -cpp -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch
else
	FC = h5pfc -fpp -r8
endif

ifeq ($(MACHINE),PC)
# GNU Debug Flags
	# FC += -O0 -g -fbacktrace -Wall -Wextra
	# FC += -Wpedantic
	# FC += -Warray-temporaries
	# FC += -fcheck=all -finit-real=snan -ffpe-trap=invalid #-std=f2018
# FC += -O0 -pg -fbacktrace -fbounds-check
# Intel Debug Flags
# FC += -O0 -g -traceback -check bounds
	ifeq ($(FLAVOUR),GNU)
		LDFLAGS = -lfftw3 -llapack -ldl
	else
		LDFLAGS = -lfftw3 -qmkl=sequential
	endif
endif
ifeq ($(MACHINE),DISCOVERER)
	ifeq ($(FLAVOUR),GNU)
		FC = h5pfc -cpp -fdefault-real-8 -fdefault-double-8
		LDFLAGS += -lfftw3 -llapack -ldl
	else
		LDFLAGS += -lfftw3 -qmkl=sequential
	endif
endif
ifeq ($(MACHINE),SNELLIUS)
	ifeq ($(FLAVOUR),GNU)
		LDFLAGS = -lfftw3 -lopenblas -ldl
	else
		LDFLAGS = -lfftw3 -qmkl=sequential
	endif
endif
ifeq ($(MACHINE),IRENE)
	FC += -mtune=skylake -xCORE-AVX512 -m64 -fPIC $(FFTW3_FFLAGS)
	LDFLAGS = $(FFTW3_LDFLAGS) $(MKL_LDFLAGS) -ldl
endif
ifeq ($(MACHINE),MARENOSTRUM)
	FC += -mtune=skylake -xCORE-AVX512 -m64 -fPIC $(FFTW_FFLAGS)
	LDFLAGS = $(FFTW_LIBS) -mkl=sequential
endif
ifeq ($(MACHINE),SUPERMUC)
	FC = mpif90 -fpp -r8 -O3 $(HDF5_INC)
	LDFLAGS = $(FFTW_LIB) $(HDF5_F90_SHLIB) $(HDF5_SHLIB) -qmkl=sequential
endif

ifeq ($(FLAVOUR),GNU)
	FC += -J $(OBJDIR)
else
	FC += -module $(OBJDIR)
endif

#=======================================================================
#  Non-module Fortran files to be compiled:
#=======================================================================
EXTRA_DIST = transpose_z_to_x.F90 transpose_x_to_z.F90 transpose_x_to_y.F90\
	     transpose_y_to_x.F90 transpose_y_to_z.F90 transpose_z_to_y.F90\
	     factor.F90 halo.F90 fft_common.F90 alloc.F90 halo_common.F90

# Object files associated with standard flow solver
OBJS = obj/main.o obj/CalcMaxCFL.o \
	obj/CalcMeanProfiles.o obj/CheckDivergence.o \
	obj/CorrectVelocity.o obj/CreateGrid.o obj/CreateInitialConditions.o \
	obj/DeallocateVariables.o obj/DebugRoutines.o obj/ExplicitTermsTemp.o \
	obj/ExplicitTermsVX.o obj/ExplicitTermsVY.o obj/ExplicitTermsVZ.o \
	obj/factorize.o obj/HdfReadContinua.o obj/HdfRoutines.o \
	obj/ImplicitAndUpdateTemp.o obj/ImplicitAndUpdateVX.o obj/ImplicitAndUpdateVY.o \
	obj/ImplicitAndUpdateVZ.o obj/InitTimeMarchScheme.o \
	obj/InitVariables.o obj/LocateLargeDivergence.o obj/MakeMovieXCut.o \
	obj/MakeMovieYCut.o obj/MakeMovieZCut.o obj/MpiAuxRoutines.o \
	obj/QuitRoutine.o obj/ReadInputFile.o obj/ResetLogs.o \
	obj/SetTempBCs.o obj/SolveImpEqnUpdate_Temp.o obj/SolveImpEqnUpdate_X.o \
	obj/SolveImpEqnUpdate_YZ.o \
	obj/TimeMarcher.o obj/WriteFlowField.o obj/WriteGridInfo.o \
	obj/CalcWriteQ.o obj/GlobalQuantities.o obj/ReadFlowInterp.o

# Object files associated with multiple resolution grids
OBJS += obj/CreateMgrdGrid.o obj/InitMgrdVariables.o \
	obj/DeallocateMgrdVariables.o obj/CreateMgrdStencil.o

# Object files associated with initial condition interpolation
OBJS += obj/CreateNewInputStencil.o obj/CreateOldGrid.o obj/CreateNewSalStencil.o \
	obj/InterpInputSal.o obj/InterpInputVel.o \
	obj/InterpVelMgrd.o obj/InitInputVars.o obj/DeallocateInputVars.o \
	obj/InterpInputPhi.o

# # Object files associated with the immersed boundary method
OBJS += obj/SolveImpEqnUpdate_Temp_ibm.o obj/SolveImpEqnUpdate_X_ibm.o \
	obj/SolveImpEqnUpdate_YZ_ibm.o obj/topogr_ibm.o obj/SolveImpEqnUpdate_Sal_ibm.o \
	obj/DeallocateIBMVars.o

# Object files for plane writing
OBJS += obj/mean_zplane.o

# Module object files
MOBJS = obj/param.o obj/decomp_2d.o obj/AuxiliaryRoutines.o obj/decomp_2d_fft.o \
	obj/pressure.o obj/HermiteInterpolations.o obj/grid.o obj/h5_tools.o obj/means.o \
	obj/ibm_param.o obj/IBMTools.o obj/moisture.o obj/salinity.o obj/phasefield.o \
	obj/time_averaging.o obj/spectra.o obj/potential_energy.o

#=======================================================================
#  Files that create modules:
#=======================================================================
MFILES = param.F90 decomp_2d.F90 AuxiliaryRoutines.F90 decomp_2d_fft.F90 \
	pressure.F90 HermiteInterpolations.F90 grid.F90 ibm_param.F90 IBMTools.F90 \
	moisture.F90 salinity.F90 phasefield.F90 time_averaging.F90 spectra.F90 \
	potential_energy.F90

#============================================================================ 
#  make PROGRAM   
#============================================================================
PROGRAM = afid 

#Compiling 
all: objdir $(PROGRAM) 
$(PROGRAM): $(MOBJS) $(OBJS) 
	$(FC) -o $@ $^ $(LDFLAGS) 

#============================================================================
#  Dependencies 
#============================================================================
$(OBJDIR)/param.o: src/flow_solver/param.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/AuxiliaryRoutines.o: src/flow_solver/AuxiliaryRoutines.F90 
	$(FC) -c -o $@ $< 
$(OBJDIR)/decomp_2d.o: src/flow_solver/2decomp/decomp_2d.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/decomp_2d_fft.o: src/flow_solver/2decomp/decomp_2d_fft.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/ibm_param.o: src/ibm/ibm_param.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/grid.o: src/grid.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/pressure.o: src/pressure.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/HermiteInterpolations.o: src/multires/HermiteInterpolations.F90 obj/ibm_param.o
	$(FC) -c -o $@ $<
$(OBJDIR)/h5_tools.o: src/h5tools/h5_tools.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/means.o: src/h5tools/means.F90 obj/ibm_param.o
	$(FC) -c -o $@ $<
$(OBJDIR)/IBMTools.o: src/ibm/IBMTools.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/salinity.o: src/salinity.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/phasefield.o: src/phasefield.F90 obj/salinity.o
	$(FC) -c -o $@ $<
$(OBJDIR)/moisture.o: src/moisture.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/time_averaging.o: src/time_averaging.F90
	$(FC) -c -o $@ $<
$(OBJDIR)/spectra.o: src/spectra.F90 obj/time_averaging.o obj/pressure.o
	$(FC) -c -o $@ $<
$(OBJDIR)/%.o: src/%.F90 $(MOBJS)
	$(FC) -c -o $@ $<
$(OBJDIR)/%.o: src/flow_solver/%.F90 $(MOBJS)
	$(FC) -c -o $@ $<
$(OBJDIR)/%.o: src/h5tools/%.F90 $(MOBJS)
	$(FC) -c -o $@ $<
$(OBJDIR)/%.o: src/multires/%.F90 $(MOBJS)
	$(FC) -c -o $@ $<
$(OBJDIR)/%.o: src/multires/IC_interpolation/%.F90 $(MOBJS)
	$(FC) -c -o $@ $<
$(OBJDIR)/%.o: src/ibm/%.F90 $(MOBJS)
	$(FC) -c -o $@ $<

#============================================================================
#  Clean up 
#============================================================================
clean: 
	/bin/rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*genmod* $(OBJDIR)/*.o obj\

.PHONY: objdir
objdir: $(OBJDIR) 
$(OBJDIR): 
	mkdir -p ${OBJDIR}
