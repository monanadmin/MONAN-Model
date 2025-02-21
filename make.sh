#!/bin/bash
#Usage: make target CORE=[core] [options]
#Example targets:
#   gnu             - GNU Fortran, C, and C++ compilers
#   xlf             - IBM XL compilers
#   xlf-summit-omp-offload - IBM XL compilers w/OpenMP offloading on ORNL Summit
#   ftn             - Cray compilers
#   titan-cray      - (deprecated) Cray compilers with options for ORNL Titan
#   nvhpc           - NVIDIA HPC SDK
#   pgi             - PGI compiler suite
#   pgi-summit      - PGI compiler suite w/OpenACC options for ORNL Summit
#   pgi-nersc       - (deprecated) PGI compilers on NERSC machines
#   pgi-llnl        - (deprecated) PGI compilers on LLNL machines
#   ifort           - Intel Fortran, C, and C++ compiler suite
#   ifort-scorep    - Intel compiler suite with ScoreP profiling library
#   ifort-gcc       - Intel Fortran compiler and GNU C/C++ compilers
#   intel-mpi       - Intel compiler suite with Intel MPI library
#   gfortran        - GNU Fortran, C, and C++ compilers
#   gfortran-clang  - GNU Fortran compiler with LLVM clang/clang++ compilers
#   g95             - (deprecated) G95 Fortran compiler with GNU C/C++ compilers
#   pathscale-nersc - (deprecated) Pathscale compilers on NERSC machines
#   cray-nersc      - (deprecated) Cray compilers on NERSC machines
#   gnu-nersc       - (deprecated) GNU compilers on NERSC machines
#   intel-nersc     - (deprecated) Intel compilers on NERSC machines
#   bluegene        - (deprecated) IBM XL compilers on BlueGene/Q systems
#   llvm            - LLVM flang, clang, and clang++ compilers
#   nag             - NAG Fortran compiler and GNU C/C++ compilers
#   cray            - Cray Programming Environment
#   intel           - Intel oneAPI Fortran, C, and C++ compiler suite
#
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
#    PRECISION=double - builds with default double-precision real kind. Default is to use single-precision.
#    SHAREDLIB=true - generate position-independent code suitable for use in a shared library. Default is false.
#. ./load_monan_app_modules.sh
export DIRhome=/home/paulo.kubota/monan_linux/monan/model
export DIRroot=/mnt/beegfs/paulo.kubota/monan_linux/monan/model

export MONAN_SRC_DIR=${DIRhome}/sources/MONAN-Model_GPU
export MONAN_EXE_MODEL=${DIRroot}/exec/MONAN-Model_GPU
export MONAN_EXE_PRE=${DIRroot}/../pre/exec/MONAN-Model_GPU 

#export NETCDF=/mnt/beegfs/monan/libs/netcdf
#export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF
export NETCDF='/mnt/beegfs/paulo.kubota/lib/lib_gnu/netcdf'
export PNETCDF='/mnt/beegfs/paulo.kubota/lib/lib_gnu/pnetcdf/pnetcdf-1.12.3/PnetCDF'
# PIO is not necessary for version 8.* If PIO is empty, MPAS Will use SMIOL
export PIO=
mkdir -p ${MONAN_SRC_DIR}/bin
mkdir -p ${MONAN_EXE_MODEL}/exec/
#######################################################################################
#
#
#
#######################################################################################
rm  ${MONAN_SRC_DIR}/bin/atmosphere_model      
rm  ${MONAN_EXE_MODEL}/exec/atmosphere_model   
  
make clean CORE=atmosphere
#make -j 8 gfortran CORE=atmosphere DEBUG=true OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee make.output
make -j 8 gfortran CORE=atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee make.output
mv ${MONAN_SRC_DIR}/atmosphere_model ${MONAN_SRC_DIR}/bin/
mv ${MONAN_SRC_DIR}/build_tables ${MONAN_SRC_DIR}/bin/
make clean CORE=atmosphere
cp -f ${MONAN_SRC_DIR}/bin/atmosphere_model      ${MONAN_EXE_MODEL}/exec/
cp -f ${MONAN_SRC_DIR}/bin/build_tables          ${MONAN_EXE_MODEL}/exec/
if [ -s "${MONAN_EXE_MODEL}/exec/atmosphere_model" ]; then
    echo ""
    echo -e "\033[1;32m==>\033[0m Files atmosphere_model generated Successfully in ${MONAN_SRC_DIR}/bin and copied to ${MONAN_EXE_MODEL}/exec !"
    echo
else
    echo -e "\033[1;31m==>\033[0m !!! An error occurred during build. Check output"
    exit -1
fi
#######################################################################################
#
#
#
#######################################################################################
rm  ${MONAN_SRC_DIR}/bin/init_atmosphere      
rm  ${MONAN_EXE_PRE}/exec/init_atmosphere   
make clean CORE=init_atmosphere
#make -j 8 gfortran CORE=init_atmosphere  DEBUG=true  OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee make.output
make -j 8 gfortran CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee make.output
mv ${MONAN_SRC_DIR}/init_atmosphere_model ${MONAN_SRC_DIR}/bin/
make clean CORE=init_atmosphere

mkdir -p ${MONAN_EXE_PRE}/exec/
cp -f ${MONAN_SRC_DIR}/bin/init_atmosphere_model ${MONAN_EXE_PRE}/exec/

if [ -s "${MONAN_EXE_PRE}/exec/init_atmosphere_model" ]; then
    echo ""
    echo -e "\033[1;32m==>\033[0m Files init_atmosphere_model generated Successfully in ${MONAN_SRC_DIR}/bin and copied to ${MONAN_EXE_PRE}/exec !"
    echo
else
    echo -e "\033[1;31m==>\033[0m !!! An error occurred during build. Check output"
    exit -1
fi

