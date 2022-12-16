#include <iostream>
#include <cmath> 
#include <blaze/Blaze.h>
#include <blaze/math/Submatrix.h>
#include <blaze/Forward.h>

#include <fstream>
#include <string>
using namespace std;
using namespace blaze;
using blaze::unaligned;
using blaze::DynamicVector ;

// ASSEMBLY  func of element matrices
// assembly matrix : project m,k to M,K according to the element unit LM((1,2), (2,3)) ---> index ((0,1),(1,2))
DynamicMatrix<double> getGlobalMatrix(const size_t &a, DynamicMatrix<double> &X, DynamicMatrix<double> &x ){ // nel, M0, m0
    for(size_t i = 0UL; i <a; i++){
        for(size_t j = 0UL; j<a; j++){
            X(i,j) = X(i,j) + x(i,j);
        }
    }

    for(size_t i = 1UL; i <= a; i++){
        for(size_t j = 1UL; j<= a; j++){
            X(i,j) = X(i,j) + x(i-1,j-1);
        }
    }
    return X;

}

// diagnoal matrice
DynamicMatrix<double> diagMatrix(DynamicMatrix<double> Y){
    DynamicMatrix<double> D(3UL, 3UL); 
    for(size_t i = 0UL; i <Y.rows(); i++){
        for(size_t j = 0UL; j< Y.columns(); j++){
            if(i == j){
                D(i,j) = Y(i,j);
            }else{
                D(i,j) = 0;
            }
        }
    }
    return D;
}


int main()
{

    // 1. (given)initial parameters for two spring obj.
    // *** ...INPUT DATA...
    const double m1 = 6.0, m2 = 3.0;
    double k1 = 1.0/6, k2 = 1.0/6;

    double test = 4.23456789123412;
    cout << "test prec " << k1<< endl;
    cout << "test prec " <<setprecision(12) << k1 << endl;



    const size_t nsd = 1; // Number of space dimensions
    const size_t ndof = 1; // Number of degrees of freedom per node (1D)
    
    const size_t nel = 2;  // number of element
    const size_t nnp = nel+1; // total number of node 3
    const int nen = 2; // nodes on each element

    // Initialize vector of  displacement,velocity , accelaration: d0, d_v, d_a : use vector???
    DynamicVector<double,columnVector > d{ 0.0, 6.0, 12.0 }; //displacement0
    DynamicVector<double> v(3UL); // Initialize velocity vector
    const int nd = 1; // Number of prescribed displacement degrees-of-freedom

    // ZeroVector<double> d_v(3UL); // column vec // velocity0
    // ZeroVector<double> d_a(3UL); // accelaration0

    // define element unit |--|--| i.e. 0--1--2 (unit1: node1-> node2; unit2: node2->node3)
    StaticMatrix<int,2UL,2UL> LM{{0,1}, {1,2}};  // did not use LM array to calculate the M, K

    // *** ...TIME ITERATION...
    const double dt = 0.01; // define the time step size
    const int nt = 8000; // the number of time steps

    // NEWMARK PARAMETERS

    double beta_b = 0.2500; 
    double gamma_b = 0.5000; 

    //////////////
    // 2.  discretization, object -> elements
    // m, k vector for each element
    // note: skip calculate process
    
    // generate element mass matrix
    DynamicMatrix<double> mel( 2UL , 2UL ); // mel from calculating(FEM), here just assign it a matrix
    mel = {{3.0, 0.0}, 
            {0.0, 3.0}};

    // generate element stiffness matrix
    DynamicMatrix<double> k( 2UL , 2UL ); // mel from calculating(FEM), here just assign it a matrix
    k = { { 1.0/6, -1.0/6},
          { -1.0/6, 1.0/6 }};
    
    // define APPLIED FORCE
    //fd=zeros(nnp,nt+1);
    // ZeroVector<double> fd(nnp, nt+1);
    DynamicMatrix <double> fd(nnp, nt+1);
    fd(2, 0) = 1.0;  // fd(3,1) = 1 --> project to c++ index rule : fd(2,0) = 1
    
    // fd.set( 2, 0, 1 );
     // ?????????????err1: why take this assignment as function operation.
    

    //////////////////////////////////////////////////////////////
    // 3.  assembly of global  matrix
    // [m] --> [M]; [k] -->[K]

    // initialize global M vector ---(use dynamic matrix in blaze)
    DynamicMatrix<double> M(nel+1 , nel+1);
    // initialize global K vector ---(use dynamic matrix in blaze)
    DynamicMatrix<double> K(nel+1 , nel+1);


    // assembly of global matrix function
    M = getGlobalMatrix(nel, M, mel);
    K = getGlobalMatrix(nel, K, k);

    std::cout << "K " << K<<endl;
    std::cout << setprecision(14) << K <<endl;
    

    // 4. split matrix: Partitioning
    //  Partition M
    DynamicMatrix<double> L_m, U_m, P_m;
    DynamicMatrix<double> D_m(3UL, 3UL);
    //ZeroMatrix<double> D_m(3UL, 3UL);
    lu( M, U_m, L_m,  P_m); // decomposition M, get L_m, U_m
    
    // D_m = diagMatrix(M, D_m);
    D_m = diagMatrix(M);  

    // std::cout << "L_m:"<< L_m <<std::endl;
    // std::cout << "U_m:"<< U_m <<std::endl; 
    // std::cout << "D_m:"<< D_m <<std::endl;  

    //  Partition K
    DynamicMatrix<double> L_k, U_k, P_k; //DynamicMatrix<double> means DynamicMatrix<double,rowMajor> ==> so lu(k,u,l)
    lu(K, U_k, L_k, P_k);

    DynamicMatrix<double> D_k(3UL, 3UL);
    D_k = diagMatrix(K);

    // **... get K, M _plus and minus matrices
    // for jacobi:
    DynamicMatrix<double> M_plus_j, M_minus_j, K_plus_j, K_minus_j;
    M_plus_j = D_m;
    M_minus_j = -(L_m + U_m);

    K_plus_j = D_k;
    K_minus_j = K_plus_j - K;
    //K_minus_j = -(L_k + U_k);
    std::cout << "K_plus_j:"<< K_plus_j <<std::endl;
    std::cout << "K_minus_j:"<< K_minus_j <<std::endl; 
    

    // for Gaus_seidel:
    DynamicMatrix<double> M_plus_g, M_minus_g, K_plus_g, K_minus_g;
   // cout << "M size   -" << size(M) << endl;
   //cout << "D_m size " << size(D_m) << endl;
   // cout << "L_m      " << size(L_m) << endl;
    M_plus_g = D_m + L_m;

    M_minus_g = -(U_m);

    K_plus_g = D_k + L_k;
    K_minus_g = -(U_k);

    
    //  Starting solving
    //****************************************************************************//
    // Jacobi
    // 1. get initial d_0, v_0, a_0;
    DynamicVector<double,columnVector > d_0, v_0, a_0, a0;
    // error0 here is just for corection
    
    d_0=d;
    v_0=v; 
    //take the initial displacement and the associated initial istantaneous acceleration as I.C. (not for the analytical solution
    // a0=M\(-K*d);  
    DynamicVector<double,columnVector> temp = -K * d;
    a0 = inv(M) * temp; // -K * d is a vector;
    cout << "a0 :" << a0 <<endl;
    a0[0] = 0.0;
    a0 = a0 ;
    a_0=a0;

    std::cout << "a0" << a0<<endl;
    std::cout << "a_0" << a_0<<endl;

    // 2. set initial d,v,a for the whole process
    DynamicMatrix <double> U_d(nnp, nt+1); // use zeroMatrix<>(), got complain
    DynamicMatrix <double> U_v(nnp, nt+1);
    DynamicMatrix <double> U_a(nnp, nt+1);
   // std::cout << "U_d:"<< U_d <<std::endl;
    
    /*
    U_d(:,1)=d;
    U_v(:,1)=v;
    U_a(:,1)=a0;
    */
    // initialize first column
    for( size_t i=0UL; i<nnp; ++i ) { 
        std::cout << "d[i]:"<< d[i] <<std::endl;    
        U_d(i,0) = d[i];
        U_v(i,0) = v[i];
        U_a(i,0) = a0[i];
        
    }
    

    //Iniziatise WR time-discrete matrix (initial WRs=initial conditions)
    DynamicMatrix <double> WR(3, nt+1);
    DynamicMatrix <double> WR2(3, nt+1);

    DynamicMatrix <double> WRv(3, nt+1);
    DynamicMatrix <double> WRa(3, nt+1);

    for( int i=0; i<3; i++ ) {
        for( size_t j=0UL; j<nt+1; j++ ) {
            WR(i,j)= d[i];
            WRv(i,j) = v_0[i];
            WRa(i,j) = a_0[i];
        }
    }

    // Initialize first row of WR_STOR=ITERATION 0, THEN COMPARED WITH ITERATION %1!
    /*
    WR_STOR(1,:)=WR(nnp,:);
    WR_STOR2(1,:)=WR2(nnp,:);
    */
    // how to make wr_store dynamic ? 随着迭代增加，一行行填充
    DynamicMatrix<double> WR_STOR(35, nt+1); 
    DynamicMatrix<double> WR_STOR2(35, nt+1); 
    int STOR_line = 0;
    
   
    for( size_t j=0UL; j<nt+1; ++j ) {
        WR_STOR(0,j) = WR(2,j);
        WR_STOR2(0,j) = WR2(2,j);
    } 
    
     
    //std::cout<<"1col "<< WR(2,0)  <<" ;2col " <<WR(2,1)<<" ;3col "<<WR(2,2)<<std::endl; 

    // let's initialize e (error) and the counting i
    double e=1.0;
    double e2=1.0;
    size_t i = 0;
    DynamicVector<double, columnVector> d1p(nnp, 1UL);
    DynamicVector<double> v1p(nnp, 1UL);
    DynamicVector<double,columnVector> b(nnp);
    DynamicMatrix<double,columnMajor> A;
    DynamicMatrix<double,columnMajor> L, U, P;
    DynamicMatrix <double> fd1(3, nt+1);
    DynamicVector<double,columnVector > a(3);
    DynamicVector<double,columnVector > z(nnp);
    DynamicVector<double,columnVector > a1(nnp);
    // corrector error1
    DynamicVector<double,columnVector > error1(nnp);
    DynamicVector<double,columnVector > d1(nnp);
    DynamicVector<double,columnVector > v1(nnp);
    DynamicVector<double,columnVector > error2(nnp);
    DynamicMatrix <double> cc(1UL, nt+1);
    DynamicMatrix <double> dd(1UL, nt+1);
    DynamicMatrix <double> cc2(1, nt+1);
    DynamicMatrix <double> dd2(1, nt+1);

    while((e>pow(10.0,-14)) || (e2>pow(10.0,-14) )){
        
        // why again initialize???
        //std::cout << "e0  " << e << std::endl;
        //std::cout << "e2_0  " << e2 << std::endl;
        d = d_0;
        v = v_0;
        a = a0;

        std::cout <<"a " << a <<endl;
       
        // SOLUTION OF THE MATRIX SYSTEM: vector force
       
        //DynamicMatrix<double> temp(nnp, 1);
        fd1 =  K_minus_j * WR;
        // std::cout <<"K_minus_j " << K_minus_j << endl;
        // std::cout << "WR0 " << WR(0,0) << endl;
        // std::cout << "WR1 " << WR(1,0) << endl;
        // std::cout << "WR2 " << WR(2,0) << endl;

        // std::cout << "fd1 " << fd1(0, 2) << endl;
        
        // for( size_t jj=0UL; jj<nt+1; ++jj ) {
        //     // fd1(ii,jj) = K_minus_j * WR(ii, jj);  per ogni colonna di WR (spostamento a un t moltiplico per Kminus che � come moltiplicare per l'1/6 di prima (coeff riduttivo/moltiplicativo))
        //     submatrix(fd1, 0UL, jj, 3UL, 1UL) = K_minus_j * submatrix(WR, 0UL, jj, 3UL, 1UL);
            
        // }
        
        //SOLUTION
        
        A = M + beta_b *(dt*dt)*K_plus_j;
        lu(A, L, U, P);
       
        
        int n = 0;
        while(n<nt){
            
            // PREDICTOR PHASE // ?????? d1p, v1p didnt change as n increment, so move it outside of while
            d1p = d + dt*v + ((dt*dt)/2) * (1-2*beta_b)*a;
            v1p = v + (1-gamma_b)*dt*a;

            // std::cout << "d1p "<< d1p << std::endl;
            // std::cout << "v1p" << v1p << std::endl;
            // std::cout << "d1p " << setprecision(16) << d1p <<endl;
            // std::cout << "v1p " << setprecision(16) << v1p <<endl;

        // for(size_t i=0UL; i<nnp; ++i){
        //     d1p[i] = d[i] + dt*v[i] + ((dt*dt)/2) * (1-2*beta_b)*a[i];
        //     v1p[i] = v[i] + (1-gamma_b)*dt*a[i];       
        // }
            
            // std::cout << "column(fd1, n+1)  " <<column(fd1, n+1) << std::endl;
            // std::cout << "K_plus_j " << K_plus_j << std::endl;
            // std::cout << "d1p  " <<d1p << std::endl;
            b = column(fd1, n+1) - K_plus_j * d1p ;
            //std::cout << "b " << setprecision(16) << b <<endl;
            // b[0] = (int)b[0] ;
            // b[1] = (int)b[1] ;
            // b[2] = (int)b[2] ;
            
            //std::cout << "b " << b << std::endl;
           
            // LU_decomposition
            // error1 = {0,0,0.000033};
            z = inv(L) * b;
            a1 = inv(U) * z ;
            a1[0] = 0;

            // std::cout << "z " << z << std::endl;
            // std::cout << "a1 " << a1 << std::endl;
            
            // CORRECTOR PHASE
            d1 = d1p + beta_b * (dt*dt)*a1;
            v1 = v1p+(1-gamma_b)*dt*a1;
            // std::cout << "d1 " << d1 << std::endl;
            // std::cout << "v1 " << v1 << std::endl;
            /**
             * SUBSTITUTING
                U_d(:,n+1)=d1;
                U_v(:,n+1)=v1;
                U_a(:,n+1)=a1;
             */
            for( size_t i=0UL; i<nnp; ++i ) {     
                U_d(i,n+1) = d1[i];
                U_v(i,n+1) = v1[i];
                U_a(i,n+1) = a1[i];

            } 
            

            d = d1;
            v = v1;
            a = a1;
            // std::cout << "d1 " << setprecision(16) << d1 <<endl;
            // std::cout << "v1 " << setprecision(16) << v1 <<endl;
            // std::cout << "a1 " << setprecision(16) << a1 <<endl;

            for( size_t i=0UL; i<nnp; ++i ) {     
                WR(i,n+1) = d1[i];
                WRv(i,n+1) = v1[i];
                WRa(i,n+1) = a1[i];
            }
            n++;
            
        }       

        i = i+1;
        cout <<"J current i: " << i << endl;
        
        STOR_line = STOR_line + 1;
        cout <<" STOR_line  : " <<  STOR_line  << endl;
        // WR_STOR(STOR_line,:) = WR(nnp,:);
        for( size_t j=0UL; j<nt+1; ++j ) {
            WR_STOR(STOR_line,j) = WR(2,j);
        }
        // std::cout << "WR_STOR " << setprecision(16) << WR_STOR(STOR_line,7996) <<endl;
        // std::cout << "WR_STOR " << setprecision(16) << WR_STOR(STOR_line,7997) <<endl;
        // std::cout << "WR_STOR " << setprecision(16) << WR_STOR(STOR_line,7998) <<endl;
        // std::cout << "WR_STOR " << setprecision(16) << WR_STOR(STOR_line,7999) <<endl;
        // std::cout << "WR_STOR " << setprecision(16) << WR_STOR(STOR_line,8000) <<endl;

        

        // std::cout << "WR_STOR_1 " << setprecision(16) << WR_STOR(STOR_line-1,7996) <<endl;
        // std::cout << "WR_STOR_1 " << setprecision(16) << WR_STOR(STOR_line-1,7997) <<endl;
        // std::cout << "WR_STOR_1 " << setprecision(16) << WR_STOR(STOR_line-1,7998) <<endl;
        // std::cout << "WR_STOR_1 " << setprecision(16) << WR_STOR(STOR_line-1,7999) <<endl;
        // std::cout << "WR_STOR_1 " << setprecision(16) << WR_STOR(STOR_line-1,8000) <<endl;
        
        for( size_t j=0UL; j<nt+1; ++j ) {
            cc(0,j) = WR_STOR(STOR_line,j) - WR_STOR(STOR_line-1,j);
        }
        //std::cout << "cc " << setprecision(16) << cc(0,7998) <<endl;
        //std::cout << "cc " << setprecision(16) << cc(0,7999) <<endl;
        //std::cout << "cc " << setprecision(16) << cc(0,8000) <<endl;
        
        // cc=WR_STOR(STOR_line,:)-WR_STOR(STOR_line-1,:);
        
        dd = abs(cc); 
        //std::cout << "dd " << setprecision(16) << dd(0,7998) <<endl;
        //std::cout << "dd " << setprecision(16) << dd(0,7999) <<endl;
        //std::cout << "dd " << setprecision(16) << dd(0,8000) <<endl;
        e = max(dd);
        std::cout << "e  " << setprecision(16)<< e << std::endl;

        // WR_STOR2
        for( size_t j=0UL; j<nt+1; j++ ) {
            WR_STOR2(STOR_line,j) = WR(1,j);
        }

        
        for( size_t j=0UL; j<nt+1; j++ ) {
            cc2(0,j) = WR_STOR2(STOR_line,j) - WR_STOR2(STOR_line-1,j);
        }
        // cc=WR_STOR(STOR_line,:)-WR_STOR(STOR_line-1,:);
        
        dd2 = abs(cc2);
        e2 = max(dd2);
        std::cout << "e2  "<< setprecision(16) << e2 << std::endl;
        
    }
        auto rs = rows( WR, { 2UL } );
        std::cout << "WR(2,.) data:  \n";
        std::ofstream fout("WR_plot_res.csv");
        fout << "dt,WR\n";
        for (std::size_t step = 0; step <= nt; step++){
            fout << dt*step<< ","
                <<WR(2,step) << "\n";
        }
        fout.close();

        // test results in the end
        std::cout << "U_d" << setprecision(16) << U_d(2,8000) <<endl;
        std::cout << "U_d" << setprecision(16) << U_d(2,7999) <<endl;
        std::cout << "U_d" << setprecision(16) << U_d(2,7998) <<endl;
        std::cout << "U_d" << setprecision(16) << U_d(2,7997) <<endl;
        std::cout << "U_d" << setprecision(16) << U_d(2,7996) <<endl;

        std::cout << "U_v" << setprecision(16) << U_v(2,8000) <<endl;
        std::cout << "U_v" << setprecision(16) << U_v(2,7999) <<endl;
        std::cout << "U_v" << setprecision(16) << U_v(2,7998) <<endl;
        std::cout << "U_v" << setprecision(16) << U_v(2,7997) <<endl;
        std::cout << "U_v" << setprecision(16) << U_v(2,7996) <<endl;

        std::cout << "U_a" << setprecision(16) << U_a(2,8000) <<endl;
        std::cout << "U_a" << setprecision(16) << U_a(2,7999) <<endl;
        std::cout << "U_a" << setprecision(16) << U_a(2,7998) <<endl;
        std::cout << "U_a" << setprecision(16) << U_a(2,7997) <<endl;
        std::cout << "U_a" << setprecision(16) << U_a(2,7996) <<endl;


    
}


