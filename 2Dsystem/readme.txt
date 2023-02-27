


jacobi_compiled folder:

//chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$  cmake -DCMAKE_DIR=home/chuanqiu/Blaze/blaze .
chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ g++ -g -I/home/chuanqiu/Blaze/blaze -llapack 2Djacobi.cpp -o 2d_jacobi

2.23 faster:
   
   chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ g++ -g -I/home/chuanqiu/Blaze/blaze -llapack -fopenmp 2Djacobi.cpp -o jaco3
   chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$export OMP_NUM_THREADS=4
   chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ ./jaco3 
   
   
   
   


