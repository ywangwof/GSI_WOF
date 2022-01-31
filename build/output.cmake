-- The C compiler identification is Intel 18.0.0.20180823
-- The CXX compiler identification is Intel 18.0.0.20180823
-- Check for working C compiler: /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/icc
-- Check for working C compiler: /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/icc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working CXX compiler: /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/icpc
-- Check for working CXX compiler: /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/icpc -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- The Fortran compiler identification is Intel
-- Check for working Fortran compiler: /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/ifort
-- Check for working Fortran compiler: /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/ifort  -- works
-- Detecting Fortran compiler ABI info
-- Detecting Fortran compiler ABI info - done
-- Checking whether /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/ifort supports Fortran 90
-- Checking whether /apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/ifort supports Fortran 90 -- yes
Build the EnKF with FV3reg module
Control path is 
-- Try OpenMP C flag = [-openmp]
-- Performing Test OpenMP_FLAG_DETECTED
-- Performing Test OpenMP_FLAG_DETECTED - Failed
-- Try OpenMP C flag = [ ]
-- Performing Test OpenMP_FLAG_DETECTED
-- Performing Test OpenMP_FLAG_DETECTED - Failed
-- Try OpenMP C flag = [-fopenmp]
-- Performing Test OpenMP_FLAG_DETECTED
-- Performing Test OpenMP_FLAG_DETECTED - Success
-- Try OpenMP CXX flag = [-openmp]
-- Performing Test OpenMP_FLAG_DETECTED
-- Performing Test OpenMP_FLAG_DETECTED - Failed
-- Try OpenMP CXX flag = [ ]
-- Performing Test OpenMP_FLAG_DETECTED
-- Performing Test OpenMP_FLAG_DETECTED - Failed
-- Try OpenMP CXX flag = [-fopenmp]
-- Performing Test OpenMP_FLAG_DETECTED
-- Performing Test OpenMP_FLAG_DETECTED - Success
-- Found OpenMP: -fopenmp  
found openmp with flag 
The hostname is  fe1
done figuring out host--fe1
-- BUILD_CORELIBS manually-specified as ON
Setting paths for Generic System
setting values for corelibs
Host is set to GENERIC
/lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF
Setting Intel flags
-- Found MPI_C: /apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/libmpifort.so;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/release_mt/libmpi.so;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/libmpigi.a;/usr/lib64/libdl.so;/usr/lib64/librt.so;/usr/lib64/libpthread.so  
-- Found MPI_CXX: /apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/libmpicxx.so;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/libmpifort.so;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/release_mt/libmpi.so;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/libmpigi.a;/usr/lib64/libdl.so;/usr/lib64/librt.so;/usr/lib64/libpthread.so  
-- Found MPI_Fortran: /apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/libmpifort.so;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/release_mt/libmpi.so;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/lib/libmpigi.a;/usr/lib64/libdl.so;/usr/lib64/librt.so;/usr/lib64/libpthread.so  
MPI version is 
MPI f90 version is 
MPI f08 version is 
Enviroment NETCDF is 
netcdf_libs is /apps/netcdf/4.2.1.1-intel/lib/libnetcdf.so
-- Found NetCDF: /apps/netcdf/4.2.1.1-intel/lib/libnetcdff.so;/apps/netcdf/4.2.1.1-intel/lib/libnetcdf.so  
-- Found CURL: /usr/lib64/libcurl.so (found version "7.29.0") 
 trying to find lapack, GENERIC
-- Looking for include file pthread.h
-- Looking for include file pthread.h - found
-- Looking for pthread_create
-- Looking for pthread_create - not found
-- Looking for pthread_create in pthreads
-- Looking for pthread_create in pthreads - not found
-- Looking for pthread_create in pthread
-- Looking for pthread_create in pthread - found
-- Found Threads: TRUE  
-- Looking for Fortran sgemm
-- Looking for Fortran sgemm - found
-- A library with BLAS API found.
-- Looking for Fortran cheev
-- Looking for Fortran cheev - found
-- A library with LAPACK API found.
setting values for corelibs
BUFR library /lfs4/BMC/wrfruc/gge/precompiled/GSILIBS/b_intel18.0.5.274_impi2018.4.274/lib/libbufr_v.a set via Environment variable
SIGIO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsigio_4.a set via Environment variable
NEMSIO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libnemsio.a set via Environment variable
CRTM library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libcrtm.a set via Environment variable
SP library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsp_d.a set via Environment variable
SP library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsp_4.a set via Environment variable
SFCIO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsfcio_4.a set via Environment variable
Setting W3EMC library via environment variable /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libw3emc_d.a
W3NCO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libw3nco_d.a set via Environment variable
W3NCO_4 library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libw3nco_4.a set via Environment variable
IP library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libip_d.a set via Environment variable
IP 4 library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libip_4.a set via Environment variable
HEY!!! ncdiag flags are -free -assume byterecl -convert big_endian
BUFR library /lfs4/BMC/wrfruc/gge/precompiled/GSILIBS/b_intel18.0.5.274_impi2018.4.274/lib/libbufr_v.a set via Environment variable
SIGIO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsigio_4.a set via Environment variable
NEMSIO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libnemsio.a set via Environment variable
CRTM library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libcrtm.a set via Environment variable
SP library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsp_d.a set via Environment variable
SP library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsp_4.a set via Environment variable
SFCIO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libsfcio_4.a set via Environment variable
Setting W3EMC library via environment variable /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libw3emc_d.a
W3NCO library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libw3nco_d.a set via Environment variable
W3NCO_4 library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libw3nco_4.a set via Environment variable
IP library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libip_d.a set via Environment variable
IP 4 library /lfs4/BMC/wrfruc/gge/precompiled/NCEPLIBS/b_intel18.0.5.274_impi2018.4.274/install/lib/libip_4.a set via Environment variable
MPI include PATH  /apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/include/gfortran;/apps/intel/compilers_and_libraries_2018.5.274/linux/mpi/intel64/include
-- Configuring done
-- Generating done
-- Build files have been written to: /lfs4/NAGAPE/hpc-wof1/ywang/EPIC2/oumap/GSI_WoF/build
