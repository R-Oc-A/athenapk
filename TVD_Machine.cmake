#========================================================================================
#========================================================================================

message(STATUS "Loading machine configuration for Prof. Tom VD machine\n"  "Last tested: 2025-03-12\n\n")

# 1.CUDA Configuration
set(Kokkos_ENABLE_CUDA ON)
set(Kokkos_ARCH_ADA89 ON)
# 2. OpenMPI Configuration
set(OpenMPI_ROOT_DIR "$OPENMPI_HOME") 

# 3. Kokkos MPI configuration
set(Kokkos_ENABLE_MPI ON)# CACHE BOOL "enable mpi")
set(Kokkos_ENABLE_SERIAL OFF)# CACHE BOOL "disable serial")
set(Kokkos_ARCH_ZEN3 ON) #CACHE BOOL "CPU architecture")

# 4. compiler settings (Using nvcc_wrapper)
set(CMAKE_CXX_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/external/Kokkos/bin/nvcc_wrapper")

# 5. HDF5
set (PARTHENON_DISABLE_HDF5 OFF )#CACHE BOOL "Enable HDF5")
set( CMAKE_PREFIX_PATH "$OPENMPI_HOME;$HDF5_HOME;$CMAKE_PREFIX_PATH")
set( HDF5_DIR "$HDF5_HOME")


# common options
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default release build")

#if(${MACHINE_VARIANT} MATCHES "cuda")
#endif()

# Setting launcher options independent of parallel or serial test as the launcher always
# needs to be called from the batch node (so that the tests are actually run on the
# compute nodes.
set(TEST_MPIEXEC mpiexec CACHE STRING "Command to launch MPI applications")
set(TEST_NUMPROC_FLAG "-np" CACHE STRING "Flag to set number of processes")
