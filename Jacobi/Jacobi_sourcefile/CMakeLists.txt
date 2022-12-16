cmake_minimum_required(VERSION 3.13)
project(FirstExample CXX)

#find_package(blaze REQUIRED)
#add_executable(FirstExample FirstExample.cpp)
#target_link_libraries( FirstExample INTERFACE blaze::blaze )

find_package(blaze)

if( blaze_FOUND )
        add_library( two_mass INTERFACE )
        target_link_libraries( two_mass INTERFACE blaze::blaze )
endif()
