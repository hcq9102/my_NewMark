
There are 4 source files of different methods: Jacobi, Gauss seidel, Newmark, Newmark_reference_solution.

#How to run Jacobi method (g++):
1).Source files:
      2Djacobi.cpp
      preprocess.hpp
      assembly.hpp
      
2) Run: g++ -O2 -DNDEBUG -g -I/home/chuanqiu/Blaze/blaze -llapack -fopenmp 2Djacobi.cpp -o 2Djacobi


// cmake

cmake -DCMAKE_BUILD_TYPE=Release -DHPX_DIR=/work/chuanqiu/2Dsystem/my_NewMark/workhpx/installhpx/lib64/cmake/HPX -Dblaze_DIR=/home/chuanqiu/Blaze/blaze/share/blaze/cmake ..

-------------------
two way to compile all files.

1. ~/sourcefile/
    four cpp files & one CMakeLists.txt file

2. ~/CMakeLists.txt
   ~/src/
      --  cpp file1 & CMakeLists.txt file1
      --  cpp file2 & CMakeLists.txt file2
      --  cpp file3 & CMakeLists.txt file3
      --  cpp file4 & CMakeLists.txt file4
-------------------------------------------------      
How to Build
// cmake -DHPX_DIR=/work/chuanqiu/2Dsystem/my_NewMark/workhpx/installhpx/lib64/cmake/HPX -Dblaze_DIR=/home/chuanqiu/Blaze/blaze/share/blaze/cmake ..
      