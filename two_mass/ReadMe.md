## 9.4 There is no error report on code anymore. But the compile message shows as follow:


### Local : Visual studio code: 

        Starting build...
      /usr/bin/clang++ -g /Users/chuanqiuhe/Blaze/two_mass_sys/two_mass.cpp -o /Users/chuanqiuhe/Blaze/two_mass_sys/two_mass -std=c++17
      Undefined symbols for architecture x86_64:
        "_dgesv_", referenced from:
            blaze::gesv(int, int, double*, int, int*, double*, int, int*) in two_mass-41c3f2.o
        "_dgetrf_", referenced from:
            blaze::getrf(int, int, double*, int, int*, int*) in two_mass-41c3f2.o
      ld: symbol(s) not found for architecture x86_64
      clang: error: linker command failed with exit code 1 (use -v to see invocation)

      Build finished with error(s).

       *  The terminal process terminated with exit code: -1. 
       *  Terminal will be reused by tasks, press any key to close it. 
       
       
       

### Rostam: The similar message shows:

          chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/test825$ g++ -I/home/chuanqiu/Blaze/blaze two_mass.cpp
      /tmp/cc6NPGrP.o: In function `blaze::getrf(int, int, double*, int, int*, int*)':
      two_mass.cpp:(.text._ZN5blaze5getrfEiiPdiPiS1_[_ZN5blaze5getrfEiiPdiPiS1_]+0x3f): undefined reference to `dgetrf_'
      /tmp/cc6NPGrP.o: In function `blaze::gesv(int, int, double*, int, int*, double*, int, int*)':
      two_mass.cpp:(.text._ZN5blaze4gesvEiiPdiPiS0_iS1_[_ZN5blaze4gesvEiiPdiPiS0_iS1_]+0x48): undefined reference to `dgesv_'
      collect2: error: ld returned 1 exit status
      
      
### 9.6 compile issue solved: 
                1. On rostam(lapack lib already exist), just add another lib link: -llapack

                g++ -I/home/chuanqiu/Blaze/blaze -llapack two_mass.cpp
        
                2. local,  On Mac, need homebrew for that: https://formulae.brew.sh/formula/lapack
                
                follow the second link, After that, install lapack
                
# 922 core_dump issue solved:

        The definition of D_m , D_k using DynamicMatrix<double>D_m,D_k  , it wont allocate memory for D_m, D_k automatically.  so when calculate it,
        it becomes empty.
        so use blaze::ZeroMatrix<double> D_m(3UL, 3UL);
               blaze::ZeroMatrix<double> D_k(3UL, 3UL);  will be fine.
               
               
               
               
               
               

                
                
