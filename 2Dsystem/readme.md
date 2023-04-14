

### 1. g++ compile

chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ g++ -g -I/home/chuanqiu/Blaze/blaze -llapack 2Djacobi.cpp -o 2d_jacobi

   2.23 faster:
   
      chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ g++ -g -I/home/chuanqiu/Blaze/blaze -llapack -fopenmp 2Djacobi.cpp -o jaco3

      chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$export OMP_NUM_THREADS=4

      chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ ./jaco3 
   
   
   3.29 add Release type: -O2 -DNDEBUG
   
      g++ -O2 -DNDEBUG -I/home/chuanqiu/Blaze/blaze -llapack -fopenmp 2Djacobi.cpp -o jaco3
   
      result:
      it will run in 80 seconds on 1 thread.
   
   ////////////////////////////////////////////////////
### 2. see the performance benefit using hpx threads
     
   Write CMakeLists.txt file.
   
   ref: 
   
      https://github.com/PeriHPX/PeriHPX/blob/main/CMakeLists.txt#L65
       ...
      https://github.com/PeriHPX/PeriHPX/blob/main/CMakeLists.txt#L121-L122
   
   -------------------cmake files update------------------------------------
   
   #### 2.1 build HPX: 
   
      cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/work/chuanqiu/2Dsystem/my_NewMark/workhpx/installhpx -DHPX_WITH_FETCH_ASIO=ON -DHPX_WITH_MALLOC=system ..
      
   #### 2.2 Write cmake files for four methods: 


         -----------------------------------------------
         two ways to compile all files.

         1. ~/sourcefile/
             four cpp files & one CMakeLists.txt file

         2. ~/CMakeLists.txt
            ~/src/
               --  cpp file1 & CMakeLists.txt file1
               --  cpp file2 & CMakeLists.txt file2
               --  cpp file3 & CMakeLists.txt file3
               --  cpp file4 & CMakeLists.txt file4
         ------------------------------------------------- 
         
         
#### How to Build

      cmake -DHPX_DIR=/work/chuanqiu/2Dsystem/my_NewMark/workhpx/installhpx/lib64/cmake/HPX -Dblaze_DIR=/home/chuanqiu/Blaze/blaze/share/blaze/cmake ..
   
   
   
   


