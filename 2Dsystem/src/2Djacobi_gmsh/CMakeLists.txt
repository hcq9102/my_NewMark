# for 2DJacobi 
cmake_minimum_required(VERSION 3.19)
project(2Dsystem LANGUAGES CXX C)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

# Add own cmake scripts
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_FLAGS "-lgmp") #link gmp lib

find_package(blaze REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
include("${blaze_DIR}/blaze-config-version.cmake")
link_libraries(blaze::blaze)

# BLAS and LAPACK
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

find_package(HPX REQUIRED)

# Find the external library
find_package(Gmsh REQUIRED)
# Add the include directory for the external library
include_directories(${Gmsh_INCLUDE_DIRS})

# Enable HPX inside blaze
add_definitions(-DBLAZE_USE_HPX_THREADS)
add_hpx_executable(2Djaco SOURCES 2Djacobi.cpp          
        DEPENDENCIES ${LAPACK_LIBRARIES} ${GMSH_LIB})


# include important directories
include_directories(../preprocess_gmsh)
set_property(TARGET 2Djaco PROPERTY CXX_STANDARD 17)



