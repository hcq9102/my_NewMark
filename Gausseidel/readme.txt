/********** on Rostam: Gausseidel method
               
               chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/gaussseidel$ cmake -DCMAKE_DIR=/home/chuanqiu/Blaze/blaze .
               chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/gaussseidel$ make
               chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/gaussseidel$ g++ -I/home/chuanqiu/Blaze/blaze -llapack gausseidel.cpp -o gaus
               
               ### Get GDB exe. file (twomassGdb)
               g++ -g -I/home/chuanqiu/Blaze/blaze -llapack gausseidel.cpp -o Gdb
               
               
               
 // exit(1) 用于停止程序运行至此，用于debug
 // matrix * vector = vector
 // Jacobi_GS__WR_dt.ipynb
      里面只是修改了输入的文件名，与Jacobi_WR_dt.ipynb 几乎一样。





1.25
I am done with implementing the one dimensional jacobi and gauss-seidel methods, and got the correct results on velocity, displacement, acceleration and plot graphs, seperately. Now, we can create a combined graph with these methods.  
//
after that, to implement two-dimensional plate with a single material.
then we can do parallelization afterwards.
