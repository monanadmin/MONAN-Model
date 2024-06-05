1) clone
    git clone https://github.com/saulorfreitas/MONAN_Phys_SRF
    git checkout develop

2) change Makefile: 
    For GFORTRAN and GNU, add " -fallow-argument-mismatch " 
    FFLAGS_OPT =  -fallow-argument-mismatch 
    FFLAGS_DEBUG = -fallow-argument-mismatch 

3)--coloque no ./bashrc estas

    #....... use the exports below for using MONAN with gfortran/gcc 11
    LIBBASE=/home/sfreitas/local_bin/monan_gnu11/
    #....... use the exports below for using MONAN with PGI
    #LIBBASE=/home/sfreitas/local_bin/monan_pgi

    export PATH=${LIBBASE}/bin:$PATH
    export LD_LIBRARY_PATH=${LIBBASE}/lib:$LD_LIBRARY_PATH
    export NETCDF=${LIBBASE}
    export PNETCDF=${LIBBASE}
    export PIO=$LIBBASE
    export MPAS_EXTERNAL_LIBS="-L${LIBBASE}/lib -lhdf5_hl -lhdf5 -ldl -lz"
    export MPAS_EXTERNAL_INCLUDES="-I${LIBBASE}/include"

    #export HDF5_DISABLE_VERSION_CHECK=1
    alias cdo='/home/sfreitas/local_bin/cdo/bin/cdo'


4) Compilar
      General MPAS build command:
      make target CORE=core <options>
      - target := clean or gnu, gfortran, nvhpc, pgi
      - core   := atmosphere or init_atmosphere
      - options:= DEBUG=true, AUTOCLEAN=true, PRECISION=single, OPENMP=true
     
     for GNU/GFORTRAN :
     make gnu CORE=init_atmosphere PRECISION=single
     make clean CORE=atmosphere
     make gnu CORE=atmosphere PRECISION=single

5) resultado da compilação:
    *******************************************************************************
    MPAS was built with default single-precision reals.
    Debugging is off.
    Parallel version is on.
    Papi libraries are off.
    TAU Hooks are off.
    MPAS was built without OpenMP support.
    MPAS was built without OpenMP-offload GPU support.
    MPAS was built without OpenACC accelerator support.
    Position-dependent code was generated.
    MPAS was built with .F files.
    The native timer interface is being used
    Using the SMIOL library.
    *******************************************************************************
    make[1]: Leaving directory '/home/sfreitas/models/MONAN_Phys_SRF'

6) Executando uma simulação

  a) gerar IC atm e static fields (topo, veg, soil) 
     mpirun -np 12 init_atmosphere_model 

  b) executar o modelo  
     mpirun -np 12 atmosphere_model 
