# for 2DJacobi 
cmake_minimum_required(VERSION 3.13)
project(2Dsystem LANGUAGES CXX C)
set(CMAKE_BUILD_TYPE Release)

find_package(blaze REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
include("${blaze_DIR}/blaze-config-version.cmake")
link_libraries(blaze::blaze)


# BLAS and LAPACK
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

find_package(HPX REQUIRED)
#if (HPX_FOUND)
#  include_directories(${HPX_INCLUDE_DIRS})
#endif()

# Enable HPX inside blaze
add_definitions(-DBLAZE_USE_HPX_THREADS)
add_hpx_executable(2Djaco SOURCES 2Djacobi.cpp DEPENDENCIES ${LAPACK_LIBRARIES})
add_hpx_executable(2Dgaus SOURCES 2Dgausidel.cpp DEPENDENCIES ${LAPACK_LIBRARIES})
add_hpx_executable(2Dnewmark SOURCES 2Dnewmark.cpp DEPENDENCIES ${LAPACK_LIBRARIES})
add_hpx_executable(2Dnewmark_finer SOURCES 2Dfiner_Newmark_refer_sol.cpp DEPENDENCIES ${LAPACK_LIBRARIES})

#target_link_libraries(2Djaco ${BLAS_LIBRARIES} ${HPX_LIBRARIES})


# include important directories
include_directories(sourceflie)



