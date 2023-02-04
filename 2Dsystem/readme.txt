


////jacobi_compiled folder:

chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$  cmake -DCMAKE_DIR=/Users/chuanqiuhe/Blaze/blaze ..
chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ g++ -g -I/home/chuanqiu/Blaze/blaze -llapack 2Djacobi.cpp -o 2d_jacobi
2.4 compiled. but got issue:
    chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ ./2d_jacobi 
    terminate called after throwing an instance of 'std::runtime_error'
      what():  Inversion of singular matrix failed
    Aborted (core dumped)


