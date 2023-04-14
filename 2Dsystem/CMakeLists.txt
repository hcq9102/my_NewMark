cmake_minimum_required(VERSION 3.20.3)
project(2Dsystem LANGUAGES CXX)
set(CMAKE_BUILD_TYPE Release)

if(HPX_DIR)
    set(HPX_DIR ${HPX_DIR} CACHE PATH "HPX_DIR")
else()
    error("HPX_DIR not set")
endif()

if(blaze_DIR)
    set(blaze_DIR ${HPX_DIR} CACHE PATH "blaze_DIR")
else()
    error("blaze_DIR not set")
endif()



include_directories(src/preprocess)

add_subdirectory(src/2Djacobi)
add_subdirectory(src/2Dgausidel)
add_subdirectory(src/2Dnewmark)
add_subdirectory(src/2DfinerNewmark)
