#include <iostream>
#include <cmath> 
#include <blaze/Blaze.h>

using namespace std;
using namespace blaze;
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
DynamicMatrix<double> diagMatrix(DynamicMatrix<double> Y, DynamicMatrix<double> D){ 
    for(size_t i = 0UL; i <Y.rows(); i++){
        for(size_t j = 0UL; j< Y.columns(); j++){
            if(i == j){
                D(i,j) = Y(i,j);
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
    const double k1 = 0.1667, k2 = 0.1667;

    const size_t nsd = 1; // Number of space dimensions
    const size_t ndof = 1; // Number of degrees of freedom per node (1D)
    
    const size_t nel = 2;  // number of element
    const size_t nnp = nel+1; // total number of node 3
    const int nen = 2; // nodes on each element

    // Initialize vector of  displacement,velocity , accelaration: d0, d_v, d_a : use vector???
    DynamicVector<double,columnVector > d{ 0, 6, 12 }; //displacement0
    ZeroVector<double> v(3UL); // Initialize velocity vector
    const int nd = 1; // Number of prescribed displacement degrees-of-freedom

    // ZeroVector<double> d_v(3UL); // column vec // velocity0
    // ZeroVector<double> d_a(3UL); // accelaration0

    // define element unit |--|--| i.e. 0--1--2 (unit1: node1-> node2; unit2: node2->node3)
    StaticMatrix<int,2UL,2UL> LM{{0,1}, {1,2}};  // did not use LM array to calculate the M, K

    // *** ...TIME ITERATION...
    const double dt = 0.01; // define the time step size
    const int nt = 8000; // the number of time steps

    // NEWMARK PARAMETERS

    double beta_b = 0.25; 
    double gamma_b = 0.5; 

    //////////////
    // 2.  discretization, object -> elements
    // m, k vector for each element
    // note: skip calculate process
    
    // generate element mass matrix
    DynamicMatrix<double> mel( 2UL , 2UL ); // mel from calculating(FEM), here just assign it a matrix
    mel = {{3.0, 0}, 
            {0, 3.0}};

    // generate element stiffness matrix
    DynamicMatrix<double> k( 2UL , 2UL ); // mel from calculating(FEM), here just assign it a matrix
    k = { { 0.1667, -0.1667},
          { -0.1667, 0.1667 }};
    
    // define APPLIED FORCE
    //fd=zeros(nnp,nt+1);
    // ZeroVector<double> fd(nnp, nt+1);
    ZeroMatrix <double> fd(nnp, nt+1);
    fd(2,0) = 1.0;  // fd(3,1) = 1 --> project to c++ index rule : fd(2,0) = 1

    

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

    // 4. split matrix: Partitioning
    //  Partition M
    DynamicMatrix<double,blaze::columnMajor> L_m, U_m, P_m, D_m;
    lu( M, L_m, U_m, P_m);
    D_m = diagMatrix(M, D_m);

    //  Partition K
    DynamicMatrix<double,blaze::columnMajor> L_k, U_k, P_k, D_k;
    lu(K, L_k, U_k, P_k);
    D_k = diagMatrix(K, D_k);

    // **... get K, M _plus and minus matrices
    // for jacobi:
    DynamicMatrix<double> M_plus_j, M_minus_j, K_plus_j, K_minus_j;
    M_plus_j = D_m;
    M_minus_j = -(L_m + U_m);

    K_plus_j = D_k;
    K_minus_j = -(L_k + U_k);

    // for Gaus_seidel:
    DynamicMatrix<double> M_plus_g, M_minus_g, K_plus_g, K_minus_g;
    M_plus_g = D_m + L_m;
    M_minus_g = -(U_m);

    K_plus_g = D_k + L_k;
    K_minus_g = -(U_k);

    
    //  Starting solving
    //****************************************************************************//
    // Jacobi
    // 1. get initial d_0, v_0, a_0;
    DynamicVector<double,columnVector > d_0, v_0, a_0, a0;
    d_0=d;
    v_0=v;
    //take the initial displacement and the associated initial istantaneous acceleration as I.C. (not for the analytical solution
    // a0=M\(-K*d);
    a0 = inv(M) * (-K * d);
    a0[0] = 0;
    a_0=a0;

    // 2. set initial d,v,a for the whole process
    ZeroMatrix <double> U_d(nnp, nt+1);
    ZeroMatrix <double> U_v(nnp, nt+1);
    ZeroMatrix <double> U_a(nnp, nt+1);
    /*
    U_d(:,1)=d;
    U_v(:,1)=v;
    U_a(:,1)=a0;
    */
    for( size_t i=0UL; i<nnp; ++i ) {     
        U_d(i,0) = d[i];
        U_v(i,0) = v[i];
        U_a(i,0) = a0[i];
    }

    //Iniziatise WR time-discrete matrix (initial WRs=initial conditions)
    ZeroMatrix <double> WR(nnp, nt+1);
    ZeroMatrix <double> WR2(nnp, nt+1);

    ZeroMatrix <double> WRv(nnp, nt+1);
    ZeroMatrix <double> WRa(nnp, nt+1);

    for( size_t i=0UL; i<nnp; ++i ) {
        for( size_t j=0UL; j<nt+1; ++j ) {
            WR(i,j) = d[i];
            WRv(i,j) = v[i];
            WRa(i,j) = a0[i];
        }
    }

    // Initialize first row of WR_STOR=ITERATION 0, THEN COMPARED WITH ITERATION %1!
    /*
    WR_STOR(1,:)=WR(nnp,:);
    WR_STOR2(1,:)=WR2(nnp,:);
    */
    blaze::DynamicMatrix<double> WR_STOR, WR_STOR2;
    int STOR_line = 0;
    
    for( size_t j=0UL; j<nt+1; ++j ) {
            WR_STOR(0,j) = WR(2,j);
            WR_STOR2(0,j) = WR2(2,j);
    }   

    // let's initialize e (error) and the counting i
    size_t e=1;
    size_t e2=1;
    size_t i = 0;
    DynamicVector<double,columnVector > a(nnp);
    while(e>pow(10,-14) || e2>pow(10,-14)){
        // why again initialize???
        d = d_0;
        v = v_0;
        a = a0;
        // SOLUTION OF THE MATRIX SYSTEM: vector force
        ZeroMatrix <double> fd1(nnp, nt+1);

        DynamicMatrix<double> temp(nnp, 1);
        // how to select one row/ column data in matrix?????????????
        for( size_t ii=0UL; ii<nnp; ++ii ) {
            for( size_t jj=0UL; jj<nt+1; ++jj ) {
                fd1(ii,jj) = K_minus_j * WR(ii, jj); // per ogni colonna di WR (spostamento a un t moltiplico per Kminus che � come moltiplicare per l'1/6 di prima (coeff riduttivo/moltiplicativo))
            }
        }

        // for( size_t jj=0UL; jj<nt+1; ++jj ) {
        //     for( size_t ii=0UL; ii<nnp; ++ii ) {
        //         temp(ii, jj) = WR(ii, jj);
        //     }
        // }
        // fd1(ii,jj) = K_minus_j * WR(ii, jj); // per ogni colonna di WR (spostamento a un t moltiplico per Kminus che � come moltiplicare per l'1/6 di prima (coeff riduttivo/moltiplicativo))

        //SOLUTION
        DynamicMatrix<double> A;
        A = M + beta_b *(dt*dt)*K_plus_j;
        DynamicMatrix<double,blaze::columnMajor> L, U, P;
        lu(A, L, U, P);

        DynamicVector<double,columnVector > d1p(nnp);
        DynamicVector<double,columnVector > v1p(nnp);
        // PREDICTOR PHASE // ?????? d1p, v1p didnt change as n increment, so move it outside of while
        for(size_t i=0UL; i<nnp; ++i){
            d1p[i] = d[i] + dt*v[i] + ((dt*dt)/2) * (1-2*beta_b)*a[i];
            v1p[i] = v[i] + (1-gamma_b)*dt*a[i];
        }

        int n = 0;
        while(n<nt+1){   
            //SOLUTION : disp('b')
            DynamicVector<double,columnVector > b(nnp);
            for(size_t i=0UL; i<nnp; ++i){
                b[i] = fd1(i, n+1) - (K_plus_j * d1p)(i,0); // K_plus_j * d1p : matrix * vector
            }
            // LU_decomposition
            DynamicVector<double,columnVector > z(nnp);
            DynamicVector<double,columnVector > a1(nnp);
            z = inv(L) * b;
            a1 = inv(U) * z;
            a1[0] = 0;

            // CORRECTOR PHASE
            DynamicVector<double,columnVector > d1(nnp);
            DynamicVector<double,columnVector > v1(nnp);

            d1 = d1p + beta_b * (dt*dt)*a1;
            v1 = v1p+(1-gamma_b)*dt*a1;
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

            for( size_t i=0UL; i<nnp; ++i ) {     
                WR(i,n+1) = d1[i];
                WRv(i,n+1) = v1[i];
                WRa(i,n+1) = a1[i];
            }
        }

        i = i+1;
        cout <<"J current i: " << i << endl;
        cout << WR << endl;

        STOR_line = STOR_line + 1;
        // WR_STOR(STOR_line,:) = WR(nnp,:);
        for( size_t j=0UL; j<nt+1; ++j ) {
            WR_STOR(STOR_line,j) = WR(2,j);
        }

        ZeroMatrix <double> cc(1UL, nt+1);
        for( size_t j=0UL; j<nt+1; ++j ) {
            cc(0,j) = WR_STOR(STOR_line,j) - WR_STOR(STOR_line-1,j);
        }

        // cc=WR_STOR(STOR_line,:)-WR_STOR(STOR_line-1,:);
        ZeroMatrix <double> dd(1, nt+1);
        ZeroMatrix <double> e(1, nt+1);
        dd = abs(cc); 
        e = max(dd);

        // WR_STOR2
        for( size_t j=0UL; j<nt+1; ++j ) {
            WR_STOR2(STOR_line,j) = WR(1,j);
        }

        ZeroMatrix <double> cc2(1, nt+1);
        for( size_t j=0UL; j<nt+1; ++j ) {
            cc2(0,j) = WR_STOR2(STOR_line,j) - WR_STOR2(STOR_line-1,j);
        }
        // cc=WR_STOR(STOR_line,:)-WR_STOR(STOR_line-1,:);
        ZeroMatrix <double> dd2(1, nt+1);
        ZeroMatrix <double> e2(1, nt+1);
        dd2 = abs(cc2);
        e2 = max(dd2);
    }
    
}

