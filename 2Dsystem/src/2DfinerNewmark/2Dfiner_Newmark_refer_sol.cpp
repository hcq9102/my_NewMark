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

using namespace std;
using namespace blaze;

using blaze::DynamicVector;


int main(int argc, char *argv[])
{
    //compute ELEMENT MATRICES and vectors and ASSEMBLY
    /*
    -----------------------------------------
     STIFFNESS MATRIX  && MASS MATRIX
    -----------------------------------------
    */
    geGlobal_k(nglx, ngly, iopt, emodule, poisson, nel, nnel, ndof);
    geGlobal_m(lumped, nel, nnel, ndof);
     
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
 
    // -------------------- finer NEWMARK PREDICTOR-CORRECTOR SOLUTION (reference solution)-----------------------------

    // set initial d,v,a for the whole process
    DynamicMatrix <double> U_dN(sdof, ntN+1, 0.0); // matrix where displacement solutions at each time step are stored
    DynamicMatrix <double> U_vN(sdof, ntN+1, 0.0); //matrix where velocity solutions at each time step are stored
    DynamicMatrix <double> U_aN(sdof, ntN+1, 0.0); //matrix where acceleration solutions at each time step are stored
    // initialize first column 
    column(U_dN, 0) = d;
    column(U_vN, 0) = v;
    column(U_aN, 0) = a; 

    //----- Boundary Conditions-----
    //imposing the B.C on the displacement
    d[0] = 0.0;   // Node 1 horizontal dof
    d[1] = 0.0;   // Node 1 vertical dof
    d[18] = 0.0;   // Node 1 horizontal dof
    d[19] = 0.0;   // Node 1 vertical dof
    d[36] = 0.0;   // Node 1 horizontal dof
    d[37] = 0.0;   // Node 1 vertical dof
    d[54] = 0.0;   // Node 1 horizontal dof
    d[55] = 0.0;   // Node 1 vertical dof
    d[72] = 0.0;   // Node 1 horizontal dof
    d[73] = 0.0;   // Node 1 vertical dof
    d[90] = 0.0;   // Node 1 horizontal dof
    d[91] = 0.0;   // Node 1 vertical dof

    // imposing the B.C on the acceleration
    v[0] = 0.0;   // Node 1 horizontal dof
    v[1] = 0.0;   // Node 1 vertical dof
    v[18] = 0.0;   // Node 1 horizontal dof
    v[19] = 0.0;   // Node 1 vertical dof
    v[36] = 0.0;   // Node 1 horizontal dof
    v[37] = 0.0;   // Node 1 vertical dof
    v[54] = 0.0;   // Node 1 horizontal dof
    v[55] = 0.0;   // Node 1 vertical dof
    v[72] = 0.0;   // Node 1 horizontal dof
    v[73] = 0.0;   // Node 1 vertical dof
    v[90] = 0.0;   // Node 1 horizontal dof
    v[91] = 0.0;   // Node 1 vertical dof

    // imposing the B.C on the acceleration
    a[0] = 0.0;   // Node 1 horizontal dof
    a[1] = 0.0;   // Node 1 vertical dof
    a[18] = 0.0;   // Node 1 horizontal dof
    a[19] = 0.0;   // Node 1 vertical dof
    a[36] = 0.0;   // Node 1 horizontal dof
    a[37] = 0.0;   // Node 1 vertical dof
    a[54] = 0.0;   // Node 1 horizontal dof
    a[55] = 0.0;   // Node 1 vertical dof
    a[72] = 0.0;   // Node 1 horizontal dof
    a[73] = 0.0;   // Node 1 vertical dof
    a[90] = 0.0;   // Node 1 horizontal dof
    a[91] = 0.0;   // Node 1 vertical dof

    int rowsize = to_use.size();
   
    blaze::DynamicVector<double> d1N(sdof,0.0);
    blaze::DynamicVector<double,columnVector> v1N(sdof,0.0);
    //blaze::DynamicVector<double,columnVector> b;
    blaze::DynamicMatrix<double,columnMajor> A;
    blaze::DynamicMatrix<double,columnMajor> L, U, P;

    // blaze::DynamicVector<double,columnVector > z(sdof);
    blaze::DynamicVector<double,columnVector > a1N(sdof,1);
    // blaze::DynamicVector<double,columnVector > a1N_int(96);
   

    // blaze::DynamicVector<double> d1N_touse(rowsize);
    // blaze::DynamicMatrix<double> K_touse(rowsize,rowsize);
    // blaze::DynamicMatrix<double> fdN_touse(rowsize, ntN+1);
  
    blaze::DynamicMatrix<double,columnMajor> A_touse(rowsize, rowsize);
    
    double corrector_factor1 = beta_b * (dtN*dtN);
    double corrector_factor2 = (1-gamma_b)*dtN;
    double corrector_factor3 = ((dtN*dtN)/2) * (1-2*beta_b);

    // calculate
    A = M + corrector_factor1*K;
        
    A_touse = rows(columns(A, to_use),to_use); 
    lu(A_touse, L, U, P);

    int n = 0;
    
    while(n<ntN){
        // PREDICTOR PHASE
        d1N = d + dtN*v + corrector_factor3*a;
        v1N = v + corrector_factor2*a;
            
           
        //SOLUTION
        // K_touse = rows(columns(K, to_use),to_use);  // rows, columns, num_rows, num_cols
            
        // //reset(d1p_touse);
        // d1N_touse = elements(d1N, to_use);

        // //clear(fd1_touse);
        // fdN_touse = rows(fdN,to_use);
            
        // b = column(fdN_touse, n+1) - K_touse * d1N_touse ;

         
        // //LU_decomposition
        // z = inv(L) * b; 
        // a1N_int = inv(U) * z ;

        // reconstruct the global acceleration a1N   
        reset(a1N);
        elements(a1N,to_use) = U*L * (column(rows(fdN,to_use), n+1) - rows(columns(K, to_use),to_use) * elements(d1N, to_use)); 
        elements(a1N,to_cancel) = 0.0;
        

        // CORRECTOR PHASE
        d1N = d1N + corrector_factor1 * a1N;
        v1N = v1N + corrector_factor2 * a1N;
        
        
        column(U_dN, n+1) = d1N;
        column(U_vN, n+1) = v1N; 
        column(U_aN, n+1) = a1N;

        d = d1N;
        v = v1N;
        a = a1N;
      
        n++;  
        cout << n << endl;
           
    } 

    // for plotting
    std::ofstream fout0("2D_finer_Newmark_U_dN_plot_res_106_faster.csv");
    fout0 << "dt,U_dN\n";
    for (std::size_t step = 0; step <= ntN; step++){
        fout0 << dtN*step<< ","
            <<U_dN(106,step) << "\n";
    }
    fout0.close();

    std::ofstream fout2("2D_finer_Newmark_U_dN_plot_res_52_faster.csv");
    fout2 << "dt,U_dN\n";
    for (std::size_t step = 0; step <= ntN; step++){
        fout2 << dtN*step<< ","
            <<U_dN(52,step) << "\n";
    }
    fout2.close();
   
    return 0;  
}
