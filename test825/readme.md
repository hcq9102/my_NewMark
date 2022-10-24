// g++ -I/home/chuanqiu/Blaze/blaze -llapack two_mass.cpp 



// Get GBD exe file: g++ -g -I/home/chuanqiu/Blaze/blaze -llapack two_mass.cpp -o twomassGdb

/////////////

debug log:

1. e,e2
check L324 &340 

WR_STOR(STOR_line,j) = WR(2,j);


2.  using following conditions: go 35 iterations

    size_t e=0;
    size_t e2=0;
    size_t i = 0;
    
    while(e<pow(10.0,-14) || e2<pow(10.0,-14)){


