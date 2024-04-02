#!/bin/bash
#Usage: make target CORE=[core] [options]
#Example targets:
#    ifort
#    gfortran
#    xlf
#    pgi
#Availabe Cores:
#    atmosphere
#    init_atmosphere
#    landice
#    ocean
#    seaice
#    sw
#    test
#Available Options:
#    DEBUG=true    - builds debug version. Default is optimized version.
#    USE_PAPI=true - builds version using PAPI for timers. Default is off.
#    TAU=true      - builds version using TAU hooks for profiling. Default is off.
#    AUTOCLEAN=true    - forces a clean of infrastructure prior to build new core.
#    GEN_F90=true  - Generates intermediate .f90 files through CPP, and builds with them.
#    TIMER_LIB=opt - Selects the timer library interface to be used for profiling the model. Options are:
#                    TIMER_LIB=native - Uses native built-in timers in MPAS
#                    TIMER_LIB=gptl - Uses gptl for the timer interface instead of the native interface
#                    TIMER_LIB=tau - Uses TAU for the timer interface instead of the native interface
#    OPENMP=true   - builds and links with OpenMP flags. Default is to not use OpenMP.
#    OPENACC=true  - builds and links with OpenACC flags. Default is to not use OpenACC.
#    USE_PIO2=true - links with the PIO 2 library. Default is to use the PIO 1.x library.
#    PRECISION=single - builds with default single-precision real kind. Default is to use double-precision.
#    SHAREDLIB=true - generate position-independent code suitable for use in a shared library. Default is false.


export NETCDF=/mnt/beegfs/monan/libs/netcdf
export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF
# PIO is not necessary for version 8.* If PIO is empty, MPAS Will use SMIOL
export PIO=


make clean CORE=atmosphere
make -j 8 gfortran CORE=atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee make-all.output

#CR: TODO: put verify here if executable was created ok
mv /mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/sources/MONAN-Model_v0.1.0/atmosphere_model /mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/execs
mv /mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/sources/MONAN-Model_v0.1.0/build_tables /mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/execs
make clean CORE=atmosphere

make clean CORE=init_atmosphere
make -j 8 gfortran CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee make-all.output

mv /mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/sources/MONAN-Model_v0.1.0/init_atmosphere_model /mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/execs
make clean CORE=init_atmosphere


if [ -s "/mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/execs/init_atmosphere_model" ] && [ -e "/mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/execs/atmosphere_model" ]; then
    echo ""
    echo -e "\033[1;32m==>\033[0m Files init_atmosphere_model and atmosphere_model generated Successfully in /mnt/beegfs/carlos.souza/repo_Monan/scripts_CD-CT/scripts/../../MONAN/execs !"
    echo
else
    echo -e "\033[1;31m==>\033[0m !!! An error occurred during build. Check output"
    exit -1
fi



