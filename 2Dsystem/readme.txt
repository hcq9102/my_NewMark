




//chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$  cmake -DCMAKE_DIR=home/chuanqiu/Blaze/blaze .
chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ g++ -g -I/home/chuanqiu/Blaze/blaze -llapack 2Djacobi.cpp -o 2d_jacobi

2.23 faster:
   
   chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ g++ -g -I/home/chuanqiu/Blaze/blaze -llapack -fopenmp 2Djacobi.cpp -o jaco3
   chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$export OMP_NUM_THREADS=4
   chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ ./jaco3 
   
   
   3.29 add release type: -O2 -DNDEBUG
   g++ -O2 -DNDEBUG -I/home/chuanqiu/Blaze/blaze -llapack -fopenmp 2Djacobi.cpp -o jaco3
   
   result:
   it will run in 80 seconds on 1 thread.
   
   //////////////////////////next stage//////////////////////////
   1. see the performance benefit using hpx threads
   
   https://github.com/PeriHPX/PeriHPX/blob/main/CMakeLists.txt#L65
   
   ...
   https://github.com/PeriHPX/PeriHPX/blob/main/CMakeLists.txt#L121-L122
   
   
   
   


