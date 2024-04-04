# Copyright (c) 2014 Thomas Heller
# Copyright (c) 2015 Andreas Schaefer
#
# SPDX-License-Identifier: BSL-1.0
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(HPX_FIND_PACKAGE ON)
cmake_policy(VERSION 3.18)

include(CMakeFindDependencyMacro)

# Forward HPX_* cache variables
include("${CMAKE_CURRENT_LIST_DIR}/HPXCacheVariables.cmake")

# include HPX cmake utilities
include("${CMAKE_CURRENT_LIST_DIR}/HPXMacros.cmake")

# include exported targets

# Asio can be installed by HPX or externally installed. In the first case we use
# exported targets, in the second we find Asio again using find_package.
if(HPX_WITH_FETCH_ASIO)
  include("${CMAKE_CURRENT_LIST_DIR}/HPXAsioTarget.cmake")
else()
  set(HPX_ASIO_ROOT "/work/chuanqiu/2Dsystem/my_NewMark/workhpx/hpx/mybuild/_deps/asio-src")
  include(HPX_SetupAsio)
endif()

# LCI
if(HPX_WITH_NETWORKING AND HPX_WITH_PARCELPORT_LCI)
  if(HPX_WITH_FETCH_LCI)
    find_dependency(Threads)
    find_dependency()
    include("${CMAKE_CURRENT_LIST_DIR}/HPXLCITarget.cmake")
  else()
    set(LCI_ROOT "")
    include(HPX_SetupLCI)
    hpx_setup_lci()
  endif()
endif()

# Eve can be installed by HPX or externally installed. In the first case we use
# exported targets, in the second we find Eve again using find_package.
if("${HPX_WITH_DATAPAR_BACKEND}" STREQUAL "EVE")
  if(HPX_WITH_FETCH_EVE)
    include("${CMAKE_CURRENT_LIST_DIR}/HPXEveTarget.cmake")
  else()
    set(EVE_ROOT "")
    include(HPX_SetupEve)
  endif()
endif()

if("${HPX_WITH_DATAPAR_BACKEND}" STREQUAL "SVE")
  if(HPX_WITH_FETCH_SVE)
    include("${CMAKE_CURRENT_LIST_DIR}/HPXSVETarget.cmake")
  else()
    set(SVE_ROOT "")
    include(HPX_SetupSVE)
  endif()
endif()

include("${CMAKE_CURRENT_LIST_DIR}/HPXInternalTargets.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/HPXTargets.cmake")

get_filename_component(
  _hpx_root_dir "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE
)

set(HPX_VERSION_STRING "1.9.0")
set(HPX_VERSION_MAJOR 1)
set(HPX_VERSION_MINOR 9)
set(HPX_VERSION_SUBMINOR 0)

set(HPX_PREFIX "${_hpx_root_dir}")
set(HPX_DEBUG_POSTFIX "d")
set(HPX_BUILD_TYPE "Release")
set(HPX_CXX_STANDARD "17")
# We explicitly set the default to 98 to force CMake to emit a -std=c++XX flag.
# Some compilers (clang) have a different default standard for cpp and cu files,
# but CMake does not know about this difference. If the standard is set to the
# .cpp default in CMake, CMake will omit the flag, resulting in the wrong
# standard for .cu files.
set(CMAKE_CXX_STANDARD_DEFAULT 98)

set(HPX_GIT_COMMIT
    "ff8ad23d200a76b1b4343da20c982c41df29ffa2"
    CACHE STRING "Revision of HPX from Git" FORCE
)

set(HPX_CXX_COMPILER
    "/opt/apps/gcc/12.2.0/bin/g++"
    CACHE STRING "CXX compiler for HPX" FORCE
)
set(HPX_CXX_COMPILER_ID
    "GNU"
    CACHE STRING "CXX compiler id for HPX" FORCE
)
set(HPX_CXX_COMPILER_VERSION
    "12.2.0"
    CACHE STRING "CXX compiler version for HPX" FORCE
)

# ##############################################################################
# Setup the imported libraries (publicly linked) #

# Allocator
set(HPX_JEMALLOC_ROOT "")
set(HPX_TCMALLOC_ROOT "")
set(HPX_TBBMALLOC_ROOT "")
# Special handle for mimalloc cause we can't specify HPX_MIMALLOC_ROOT as a HINT
# to find_package
set(HPX_MIMALLOC_ROOT "")
if(NOT MIMALLOC_ROOT AND NOT "$ENV{MIMALLOC_ROOT}")
  set(MIMALLOC_ROOT ${HPX_MIMALLOC_ROOT})
endif()
include(HPX_SetupAllocator)

include(HPX_SetupThreads)

# Boost Separate boost targets to be unarily linked to some modules
set(HPX_BOOST_ROOT "/opt/apps/gcc12/boost/1.80.0/release")
# By default BOOST_ROOT is set to HPX_BOOST_ROOT (not necessary for PAPI or
# HWLOC cause we are specifying HPX_<lib>_ROOT as an HINT to find_package)
if(NOT BOOST_ROOT AND NOT "$ENV{BOOST_ROOT}")
  set(BOOST_ROOT ${HPX_BOOST_ROOT})
endif()
include(HPX_SetupBoost)
include(HPX_SetupBoostFilesystem)
include(HPX_SetupBoostIostreams)

# HIP
include(HPX_SetupHIP)

# Hwloc
set(HPX_HWLOC_ROOT "/opt/apps/hwloc/2.7.1")
include(HPX_SetupHwloc)

# Papi
set(HPX_PAPI_ROOT "")
include(HPX_SetupPapi)

# CUDA
include(HPX_SetupCUDA)

include(HPX_SetupMPI)
if((HPX_WITH_NETWORKING AND HPX_WITH_PARCELPORT_MPI) OR HPX_WITH_ASYNC_MPI)
  hpx_setup_mpi()
endif()

# APEX
set(APEX_WITH_MSR "")
set(MSR_ROOT "")
set(APEX_WITH_ACTIVEHARMONY "")
set(ACTIVEHARMONY_ROOT "")
set(APEX_WITH_OTF2 "")
set(OTF2_ROOT "")
include(HPX_SetupApex)
# ##############################################################################

set(HPX_WITH_MALLOC_DEFAULT system)

if(("${HPX_WITH_DATAPAR_BACKEND}" STREQUAL "VC") AND NOT Vc_DIR)
  set(Vc_DIR "")
endif()
include(HPX_SetupVc)

if(NOT HPX_CMAKE_LOGLEVEL)
  set(HPX_CMAKE_LOGLEVEL "WARN")
endif()

hpx_check_compiler_compatibility()
hpx_check_boost_compatibility()
hpx_check_allocator_compatibility()

if(NOT DEFINED ${CMAKE_FIND_PACKAGE_NAME}_FOUND)
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND true)
endif()

# Set legacy variables for linking and include directories
set(HPX_LIBRARIES "HPX::hpx;")
# All properties are propagated from HPX::hpx, so the following can be empty
set(HPX_INCLUDE_DIRS "")
set(HPX_LIBRARY_DIR "")
set(HPX_LINK_LIBRARIES "")
set(HPX_LINKER_FLAGS "")
