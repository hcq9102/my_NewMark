# Visual studio code: There is no error report on code anymore. But the compile message shows as follow:

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
       
       
       

# The similar rostam message shows:

          chuanqiu@rostam1:/work/chuanqiu/Dynamic_project/test825$ g++ -I/home/chuanqiu/Blaze/blaze two_mass.cpp
      /tmp/cc6NPGrP.o: In function `blaze::getrf(int, int, double*, int, int*, int*)':
      two_mass.cpp:(.text._ZN5blaze5getrfEiiPdiPiS1_[_ZN5blaze5getrfEiiPdiPiS1_]+0x3f): undefined reference to `dgetrf_'
      /tmp/cc6NPGrP.o: In function `blaze::gesv(int, int, double*, int, int*, double*, int, int*)':
      two_mass.cpp:(.text._ZN5blaze4gesvEiiPdiPiS0_iS1_[_ZN5blaze4gesvEiiPdiPiS0_iS1_]+0x48): undefined reference to `dgesv_'
      collect2: error: ld returned 1 exit status
