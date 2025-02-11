.SUFFIXES: .F .o

all: dummy physics_monan

dummy:
	echo "****** compiling physics_monan ******"

OBJS = \
	modVector.o                    \
	modConstants.o                 \
	modGate.o                      \
	modHenrysLawCts.o              \
	lib_module_cu_gf.monan.o       \
	module_cu_gf.monan.o           \
	module_mp_cldfra_cb.o          \
	module_bl_mixingscalars.o 

physics_monan: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_cu_gf.monan.o:   \
	modVector.o             \
	modConstants.o          \
	modGate.o               \
	modHenrysLawCts.o       \
	lib_module_cu_gf.monan.o


clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
endif
