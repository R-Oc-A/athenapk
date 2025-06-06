#========================================================================================
# Parthenon performance portable AMR framework
# Copyright(C) 2020 The Parthenon collaboration
# Licensed under the 3-clause BSD License, see LICENSE file for details
#========================================================================================
# (C) (or copyright) 2020. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001 for Los
# Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
# for the U.S. Department of Energy/National Nuclear Security Administration. All rights
# in the program are reserved by Triad National Security, LLC, and the U.S. Department
# of Energy/National Nuclear Security Administration. The Government is granted for
# itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works, distribute copies to
# the public, perform publicly and display publicly, and to permit others to do so.
#========================================================================================

message(STATUS "Loading machine configuration for Ricardo's laptop. Last tested: 2025-05-26\n\n")

# No CUDA
set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Disable cuda for my machine")

# OpenMPI config
# before doing anything I'd like to see if this thing works with my previous athena env
# set(OpenMPI_ROOT_DIR "$OPENMPI_HOME")

# Kokkos MPI configuration
set(Kokkos_ENABLE_MPI ON)
set(Kokkos_ENABLE_SERIAL OFF)
set(Kokkos_ARCH_ZEN3 ON)

# compiler settings
#set(CMAKE_CXX_COMPILER "${CMAKE_CURRENT_SOURCE_DIR}/external/Kokkos/bin/nvcc_wrapper")

# HDF5
set(PARTHENON_DISABLE_HDF5 OFF)
set(CMAKE_PREFIX_PATH "$HDF5_HOME;$CMAKE_PREFIX_PATH")
set(HDF5_DIR "$HDF5_HOME")

# Common options
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default release build")



