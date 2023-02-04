cmake_minimum_required(VERSION 3.13)
project(2Dproject LANGUAGES CXX C)


#find_package(blaze REQUIRED)
#add_executable(FirstExample FirstExample.cpp)
#target_link_libraries( FirstExample INTERFACE blaze::blaze )

find_package(blaze REQUIRED)

if( blaze_FOUND )
        add_library( 2Djacobi 2Djacobi.cpp )
        add_executable(2Djacobi_exe 2Djacobi.cpp)
        target_link_libraries(2Djacobi_exe INTERFACE blaze::blaze )
endif()


#preprocess.h getRes.hpp feglqd2.hpp feasmbl1_1m.hpp federiv2.hpp feeldof_m.hpp feglqd1xy.hpp feisoq4.hpp fejacobi2.hpp fematiso.hpp