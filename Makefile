MACHINE=PC
FLAVOUR=GNU

# Object and module directory:
OBJDIR=obj


FC =mpiifort -fpp
FC += -r8 -O3 -march=core-avx2 -I/home/software/hdf5-1.12.0/include
FC += -module $(OBJDIR)
LDFLAGS = -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lhdf5_fortran -L/home/software/hdf5-1.12.0/lib -lhdf5  -lz -ldl -lm

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
	obj/time_averaging.o obj/spectra.o

#=======================================================================
#  Files that create modules:
#=======================================================================
MFILES = param.F90 decomp_2d.F90 AuxiliaryRoutines.F90 decomp_2d_fft.F90 \
	pressure.F90 HermiteInterpolations.F90 grid.F90 ibm_param.F90 IBMTools.F90 \
	moisture.F90 salinity.F90 phasefield.F90 time_averaging.F90 spectra.F90

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
