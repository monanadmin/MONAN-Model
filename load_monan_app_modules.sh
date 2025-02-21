module purge
module load ohpc
module unload openmpi4
module load phdf5
module load netcdf
module load netcdf-fortran
module load mpich-4.0.2-gcc-9.4.0-gpof2pv
module load hwloc
module list

export OMP_NUM_THREADS=1
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"
export PMIX_MCA_gds=hash

export NETCDF=/mnt/beegfs/monan/libs/netcdf
export PNETCDF=/mnt/beegfs/monan/libs/PnetCDF

MPI_PARAMS="-iface ib0 -bind-to core -map-by core"
export MKL_NUM_THREADS=1
export I_MPI_DEBUG=5
export MKL_DEBUG_CPU_TYPE=5
export I_MPI_ADJUST_BCAST=12 ## NUMA aware SHM-Based (AVX512)

