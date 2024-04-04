#include <iostream>
#include <hpx/hpx_main.hpp>  // Need main source file
#include <cmath>
#include <blaze/Blaze.h>
#include <blaze/Math.h>
#include <blaze/math/Submatrix.h>
#include <blaze/Forward.h>
#include <blaze/math/Columns.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Subvector.h>
#include <blaze/math/Elements.h>

#include "preprocess.hpp"
#include "assembly.hpp"

#include <fstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace blaze;

using blaze::DynamicVector;


int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();
    //compute ELEMENT MATRICES and vectors and ASSEMBLY
    /*
    -----------------------------------------
     STIFFNESS MATRIX  && MASS MATRIX
    -----------------------------------------
    */
    geGlobal_k(nglx, ngly, iopt, emodule, poisson, nel, nnel, ndof);
    geGlobal_m(lumped, nel, nnel, ndof);

   // -------------------- gauss solution -----------------------------// 

   //Initial acceleration
    blaze::DynamicVector<double> a0(sdof,0.0);
    
    DynamicVector<double,columnVector> temp = -K * d;
    a = inv(M) * temp;
    
    int iii = 1;
    for (int jj = 0; jj < gcoord.rows();++jj) {
        if (gcoord(jj, 0) == 0.0) {
            a[iii-1] = 0.0;
            a[iii] = 0.0;
        }
        iii += 2;
    }
    
    a0 = a;
    
    //fd(:,1);

    // Compute Kplus as the lower triangular part of K
    blaze::DynamicMatrix<double> Kplus_g,Kminus_g;
    Kplus_g  = lowerMatrix(K);
    Kminus_g = Kplus_g - K;

    //2.Initialization
    DynamicVector<double,columnVector> d_0, v_0, a_0;  
    d_0=d;        // vector of the initial imposed displacement = [0;0;1.5;0;3.0;0;4.5;0;6.0;0;7.5;0;9.0;0;10.5;0;12;0...]
    v_0=v;        // vector of the initial imposed velocity = [0;0;...0;0]
    a_0=a0;       // vector of the initial imposed acceleration = [0;0;...0;0]

    // set initial d,v,a for the whole process
    DynamicMatrix <double> U_d(sdof, nt+1, 0.0); // matrix where displacement solutions at each time step are stored
    DynamicMatrix <double> U_v(sdof, nt+1, 0.0); //matrix where velocity solutions at each time step are stored
    DynamicMatrix <double> U_a(sdof, nt+1, 0.0); //matrix where acceleration solutions at each time step are stored
    // initialize first column 
    column(U_d, 0) = d;
    column(U_v, 0) = v;
    column(U_a, 0) = a0; 
    /*
    ---------------------------------------------------------------
    Iniziatise WR time-discrete matrix (initial WRs=initial conditions)
    ---------------------------------------------------------------
    */ 
    DynamicMatrix <double> WR(sdof, nt+1,0.0);
    for (int qq = 0; qq < sdof; ++qq)
    {
        row(WR,qq) = d[qq];
    }

     /*
    ---------------------------------------------------------------
    initialize first row of WR_STOR=ITERATION 0 for all dofs
    ---------------------------------------------------------------
    */ 
    DynamicMatrix <double> WR_STOR_1(WR);
    DynamicMatrix<double> WR_STOR_2(sdof, nt+1, 0.0);

    /*
    ---------------------------------------------------------------
    Initialize error 'e' and the counting 'i' of the number of iterations performed
    ---------------------------------------------------------------
    */ 
    double e = 1.0;
    int i = 0;

    /*
    ---------------------------------------------------------------
   SOLUTION OF THE WAVEFORM RELAXATION STEPS
    ---------------------------------------------------------------
    */ 
    int rowsize = to_use.size();
   
    DynamicVector<double> d1p(sdof,0.0);
    DynamicVector<double,columnVector> v1p(sdof,0.0);
    //DynamicVector<double,columnVector> b;
    DynamicMatrix<double,columnMajor> A;
    DynamicMatrix<double,columnMajor> L, U, P;

    //DynamicVector<double,columnVector > z(sdof);
    DynamicVector<double,columnVector> a1(sdof,0.0);
    //DynamicVector<double,columnVector > a_int(96);
   
    DynamicVector<double> d1(sdof,0.0);
    DynamicVector<double> v1(sdof,0.0);
    DynamicVector<double> v0(sdof,0.0);

    // blaze::DynamicVector<double> d1p_touse(rowsize);
    // blaze::DynamicMatrix<double> K_plus_g_touse(rowsize,rowsize);
    // blaze::DynamicMatrix<double> fd1_touse(rowsize, nt+1);
  
    DynamicMatrix<double,columnMajor> A_touse(rowsize, rowsize);
    DynamicMatrix <double> fd1;
    DynamicVector <double> dd(nt+1,0.0);
    DynamicVector<double> e_t(nt+1,0.0);

    double corrector_factor1 = beta_b * (dt*dt);
    double corrector_factor2 = (1-gamma_b)*dt;
    double corrector_factor3 = ((dt*dt)/2) * (1-2*beta_b);

    // Solution (time steps with Newmark)
    A = M + corrector_factor1*Kplus_g;  
    A_touse = rows(columns(A, to_use),to_use); 
    
    lu(A_touse, L, U, P);


    blaze::invert(L);
    blaze::invert(U);

    while((e > pow(10.0,-14))){
        
        //Initial conditions (reimposed at the beginning of every cycle)
        d = d_0;
        v = v0; 
        a = a0;

        //Force vector
        //blaze::DynamicMatrix <double> fd1(fd1_0);     //forces matrices (over dofs and time)
        clear(fd1); 
        fd1 = Kminus_g * WR;      
        
        /*
        -----------------------------------
        Solution (time steps with Newmark)
        ----------------------------------- 
        */
        

        int n = 0;
        while(n<nt){
            // PREDICTOR PHASE
            d1p = d + dt*v + corrector_factor3 * a;

            v1p = v + corrector_factor2 * a;
            
            
            //SOLUTION
            // K_plus_g_touse = rows(columns(Kplus_g, to_use),to_use);  // rows, columns, num_rows, num_cols
            
            // //reset(d1p_touse);
            // d1p_touse = elements(d1p, to_use);

            // //clear(fd1_touse);
            // fd1_touse = rows(fd1,to_use);
            
            //b = column(fd1_touse, n+1) - K_plus_g_touse * d1p_touse ;

            //LU_decomposition
            // z = inv(L) * (column(rows(fd1,to_use), n+1) - K_plus_g_touse * d1p_touse); // 211 iterations
            // z = L      * (column(rows(fd1,to_use), n+1) - K_plus_g_touse * d1p_touse);// 213iterations: inv(L) -->L
            
            /*
            // in this way, its correct as well(211 iterations), but cannt combine these two lines in one. so its slow.
            z = blaze::solve( L, (column(rows(fd1,to_use), n+1) - rows(columns(Kplus_g, to_use),to_use) * elements(d1p, to_use)));
            a_int = blaze::solve(U,z) ;
            
            */
           // reconstruct the global acceleration
            reset(a1);
            elements(a1,to_use) = U*L * (column(rows(fd1,to_use), n+1) - rows(columns(Kplus_g, to_use),to_use) * elements(d1p, to_use)); 
            elements(a1,to_cancel) = 0.0;

            // CORRECTOR PHASE
            d1 = d1p + corrector_factor1 * a1;
            v1 = v1p + corrector_factor2 * a1;
            
            column(U_d, n+1) = d1;
            column(U_v, n+1) = v1; 
            column(U_a, n+1) = a1;

            d = d1;
            v = v1;
            a = a1;
         
            column(WR, n+1) = d1;
            n++;     
        }  
        
        i = i+1;

        if( i % 2 == 0){ //iteration number is even
            e = 1;
        }else{ //iteration number is odd 
            WR_STOR_2 = WR;
            auto diff = WR_STOR_2 - WR_STOR_1;
            for (int n1 = 0; n1 < nt + 1; ++n1) {
                e_t[n1] = blaze::l2Norm(blaze::column(diff, n1));
            }
            
            dd = abs(e_t);
            e = max(dd);
            std::cout<<"i : "<< i << std::endl;
            std::cout<<"e : "<< e << std::endl;
            WR_STOR_1=WR;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Gaus_runtime " << elapsed.count() << " s\n";


    // for plotting
    std::ofstream fout1("2D_Gaus_WR_plot_res_52_faster.csv");
    fout1 << "dt,WR\n";
    for (std::size_t step = 0; step <= nt; step++){
        fout1 << dt*step<< ","
             <<WR(52,step) << "\n";
    }
    fout1.close();

    return 0; 


}