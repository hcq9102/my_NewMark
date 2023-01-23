#include <iostream>
#include <cmath>
#include <blaze/Blaze.h>
#include <blaze/math/Submatrix.h>
#include <blaze/Forward.h>
#include <blaze/math/Columns.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Subvector.h>

//#include "assembly.cpp"

#include <fstream>
#include <string>
using namespace std;
using namespace blaze;
using blaze::unaligned;
using blaze::DynamicVector;

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

// decompose
// diagnoal matrice
DynamicMatrix<double> lowerMatrix(DynamicMatrix<double> Y){
    DynamicMatrix<double> X(3UL, 3UL);
    for(size_t i = 0UL; i <Y.rows(); i++){
        for(size_t j = 0UL; j< Y.columns(); j++){
            if(i >= j){
                X(i,j) = Y(i,j);
            }
        }
    }
    return X;
}

int main()
{
    // 1. (given)initial parameters for two spring obj.
    // *** ...INPUT DATA...
    const double m1 = 6.0, m2 = 3.0;
    double k1 = 1.0/6, k2 = 1.0/6;

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
    std::cout << "M " << M << endl;
    K = getGlobalMatrix(nel, K, k);
    std::cout << "K \n" << K << endl;
    // 4. split matrix: Partitioning
    //  Partition M
    DynamicMatrix<double> L_m, U_m, P_m;
    DynamicMatrix<double> D_m(3UL, 3UL);
    //ZeroMatrix<double> D_m(3UL, 3UL);
    lu( M, L_m,U_m,P_m); // decomposition M, get L_m, U_m
    
    // D_m = diagMatrix(M, D_m);
    D_m = diagMatrix(M);

    // std::cout << "L_m:"<< L_m <<std::endl;
    // std::cout << "U_m:"<< U_m <<std::endl;
    // std::cout << "D_m:"<< D_m <<std::endl;

    //  Partition K
    DynamicMatrix<double> L_k, U_k, P_k; //DynamicMatrix<double> means DynamicMatrix<double,rowMajor> ==> so lu(k,u,l)
    lu(K,L_k,U_k,P_k);   // llu(K, L);

    std::cout << "P_k_g \n" << P_k << endl;
    std::cout << "L_k_g \n" << L_k << endl;
    std::cout << "U_k_g \n" << U_k << endl;
    

    DynamicMatrix<double> D_k(3UL, 3UL);
    D_k = lowerMatrix(K);
    std::cout << "D_k_g \n" << D_k << endl;


    // **... get K, M _plus and minus matrices
    // for jacobi:

    // for Gaus_seidel:
    DynamicMatrix<double> M_plus_g, M_minus_g, K_plus_g, K_minus_g;

    M_plus_g = D_m + L_m;
    M_minus_g = -(U_m);

    K_plus_g = D_k;
    std::cout << "K_plus_g " << K_plus_g << endl;
    K_minus_g = K_plus_g - K;

    //  Starting solving
    //**********************************************//
    // GAUSS-SEIDEL
    // 1. get initial d_0, v_0, a_0;
    DynamicVector<double,columnVector > d_0, v_0, a_0, a0;
    // error0 here is just for corection

    d_0=d;
    v_0=v;
    //take the initial displacement and the associated initial istantaneous acceleration as I.C. (not for the analytical solution
    // a0=M\(-K*d);
    DynamicVector<double,columnVector> temp = -K * d;
    a0 = inv(M) * temp; // -K * d is a vector;
    //cout << "a0 :" << a0 <<endl;
    a0[0] = 0.0;
    a0 = a0 ;
    a_0=a0;

    std::cout << "I am here2" <<endl;
   // std::cout << "a0" << a0<<endl;
   // std::cout << "a_0" << a_0<<endl;

    // 2. set initial d,v,a for the whole process
    DynamicMatrix <double> U_d(nnp, nt+1); // use zeroMatrix<>(), got complain
    DynamicMatrix <double> U_v(nnp, nt+1);
    DynamicMatrix <double> U_a(nnp, nt+1);
    // std::cout << "U_d:"<< U_d <<std::endl;


    // initialize first column
    for( size_t i=0UL; i<nnp; ++i ) {
        //std::cout << "d[i]:"<< d[i] <<std::endl;
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
    DynamicMatrix<double> WR_STOR(20, nt+1);
    DynamicMatrix<double> WR_STOR2(20, nt+1);
    int STOR_line = 0;


    for( size_t j=0UL; j<nt+1; ++j ) {
        WR_STOR(0,j) = WR(2,j);
        WR_STOR2(0,j) = WR(1,j);
    }


    //std::cout<<"1col "<< WR(2,0)  <<" ;2col " <<WR(2,1)<<" ;3col "<<WR(2,2)<<std::endl;

    // let's initialize e (error) and the counting i
    double e=1.0;
    double e2=1.0;
    size_t i = 0;
    DynamicVector<double, columnVector> d1p(nnp, 1UL);
    DynamicVector<double> v1p(nnp, 1UL);
    DynamicVector<double,columnVector> b(2);
    DynamicMatrix<double,columnMajor> A;
    DynamicMatrix<double,columnMajor> L, U, P;
    DynamicMatrix <double> fd1(3, nt+1);
    DynamicVector<double,columnVector > a(3);
    DynamicVector<double,columnVector > z(nnp);
    DynamicVector<double,columnVector > a1(nnp);


    DynamicVector<double,columnVector > d1(nnp);
    DynamicVector<double,columnVector > v1(nnp);
    DynamicMatrix <double> cc(1UL, nt+1);
    DynamicMatrix <double> dd(1UL, nt+1);
    DynamicMatrix <double> cc2(1, nt+1);
    DynamicMatrix <double> dd2(1, nt+1);

    while((e>pow(10.0,-16)) || (e2>pow(10.0,-16) )){
        // why again initialize???
        //std::cout << "e0  " << e << std::endl;
        //std::cout << "e2_0  " << e2 << std::endl;
        d = d_0;
        v = v_0;
        a = a0;

        // SOLUTION OF THE MATRIX SYSTEM: vector force

        //DynamicMatrix<double> temp(nnp, 1);
        fd1 =  K_minus_g * WR;
        // std::cout <<"K_minus_g" << K_minus_g << endl;
        // std::cout << "WR0 " << WR(0,0) << endl;
        // std::cout << "WR1 " << WR(1,0) << endl;
        // std::cout << "WR2 " << WR(2,0) << endl;

        //SOLUTION

        A = M + beta_b *(dt*dt)*K_plus_g;
        // A=A(2:3,2:3);
        // Creating a dense submatrix of size 2*2, starting in row 2 and column 2 (runtime arguments)
        A = submatrix(A, 1UL, 1UL, 2UL, 2UL);
        lu(A, L, U, P);
        // std::cout << "A_sub " << A << endl;
        // std::cout << "L \n " << L << endl;
        // std::cout << "U \n " << U << endl;

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

            //b = column(submatrix(fd1,1UL, 1UL, 2UL, n+1), n+1) - submatrix( K_plus_g,1UL, 1UL, 2UL, 2UL)  * submatrix(d1p,1UL, 1UL, 2UL, 2UL) ;
            // b=fd1(2:3,n+1)-Kplus(2:3,2:3)*d1p(2:3);
            // std::cout << "test subm \n "<<submatrix(fd1,1UL,n,2UL,1UL) << std::endl;
            // std::cout << "test tem1\n "<<submatrix(K_plus_g,1UL,1UL,2UL,2UL) << std::endl;
            // std::cout << "test tem2\n "<<subvector( d1p, 1UL, 2UL) << std::endl;
            // std::cout << "mult "<< submatrix(K_plus_g,1UL,1UL,2UL,2UL)*subvector( d1p, 1UL, 2UL) <<std::endl;
            DynamicVector<double> temp = (submatrix(K_plus_g,1UL,1UL,2UL,2UL)) * subvector( d1p, 1UL, 2UL);
            //DynamicVector<double> temp2 = subvector(fd1,1UL,n,2UL,1UL);
            //std::cout << "test fd1 \n "<<{fd1(1,n),fd1(2,n)} << std::endl;
            DynamicVector<double> temp3 = {fd1(1,n+1),fd1(2,n+1)};
            //std::cout << "test tem \n "<<temp << std::endl;
            

            b = temp3 -temp;
            
            // LU_decomposition
            z = inv(L) * b;
            a1 = inv(U) * z ;
            // a1=[0;a1];
            a1 = {0.0,a1[0],a1[1]}; 
            
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
             // submatrix?
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
        std::cout << "e 2 " << setprecision(16)<< e2 << std::endl;

    }
    auto rs = rows(WR, { 2UL } );
    std::cout << "WR(2,.) data:  \n";
    std::ofstream fout("GS_WR_plot_res.csv");
    fout << "dt,WR\n";
    for (std::size_t step = 0; step <= nt; step++){
        fout << dt*step<< ","
             <<WR(2,step) << "\n";
    }
    fout.close();

    // test results in the end
    //std::cout << "U_d" << setprecision(16) << U_d(2,8001) <<endl;
    std::cout << "U_d" << setprecision(16) << U_d(2,8000) <<endl;
    std::cout << "U_d" << setprecision(16) << U_d(2,7999) <<endl;
    std::cout << "U_d" << setprecision(16) << U_d(2,7998) <<endl;
    std::cout << "U_d" << setprecision(16) << U_d(2,7997) <<endl;
    std::cout << "U_d" << setprecision(16) << U_d(2,7996) <<endl;

    //std::cout << "U_d" << setprecision(16) << U_v(2,8001) <<endl;
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

    std::cout << "A " << A <<endl; // good
    std::cout << "a \n " << a <<endl;



}
