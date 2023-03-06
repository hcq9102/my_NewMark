#include <iostream>
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
//#include "assembly.hpp"

#include <fstream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace blaze;

using blaze::DynamicVector;


int main()
{
    //compute ELEMENT MATRICES and vectors and ASSEMBLY

    /*
    -----------------------------------------
     STIFFNESS MATRIX
    -------------------------------------------
    */
    feglqd2(nglx,ngly);                       // [point2,weight2]=feglqd2(nglx,ngly);  
    fematiso(iopt,emodule,poisson);           //constitutive matrix matmtx=fematiso(iopt,emodule,poisson);

    double x = 0;
    double y = 0;
    double wtx =0.0;
    double wty =0.0;
    double detjacob=0.0;
    //int count = 0;

    for (int iel = 0; iel < nel; ++iel) { // loop over the total number of elements
        DynamicMatrix <double> kel = kel0;      //initialization of element (stiffness) matrix as all zeros
        for (int i = 0; i < nnel; ++i) {
            nd[i] = nodes(iel, i); // extract nodes for the (iel)-th element
            xcoord[i] = gcoord(nd[i]-1, 0); // extract x value of the node
            ycoord[i] = gcoord(nd[i]-1, 1); // extract y value of the node
        } 
         /*
         -----------------------
            numerical integration
          -----------------------
         */   
        // sampling point(it's in r-axis)/weight(it's ins-axis)  in x-axis; y-axis
        for (int intx = 0; intx < nglx; intx++) {
            x = point2(intx,0);
            wtx = weight2(intx,0);

            for (int inty = 0; inty < ngly; inty++) {
                y = point2(inty,1);
                wty = weight2(intx,1);


                // compute shape functions and derivatives at sampling points
                //blaze::DynamicVector<double> dhdx(nnel), dhdy(nnel), dhdr(nnel), dhds(nnel);
                feisoq4(x,y);
                // compute Jacobian

                jacobi2 = fejacobi2(nnel,dhdr,dhds,xcoord,ycoord);  

                // determinant of Jacobian
                detjacob= det(jacobi2);
                //std::cout << "detjacob " <<detjacob << std::endl;
                //inverse of Jacobian matrix
                blaze::DynamicMatrix<double> invjacob = inv(jacobi2);

                //derivatives w.r.t. physical coordinate
                federiv2(nnel,dhdr,dhds,invjacob);

                //compute kinematic matrix
                fekine2D(nnel,dhdx,dhdy);

                /*
                -----------------------------------
                compute ELEMENT (stiffness) matrix
                -----------------------------------
                */
                kel = kel + (trans(kinmtx2) * matmtx) * kinmtx2 * (wtx*wty *detjacob);
            }
        }
     
        //extract system dofs for the element
        indexk = feeldofk(nd,nnel,ndof);                     
        //assemble element matrices in the global one
        K=feasmbl1_m(K,kel,indexk); //K=feasmbl1(K,kel,index);  
    }

    /*
    -----------------------------------------
     MASS MATRIX
    -------------------------------------------
    */
    /*
    -----------------------------------
    ELEMENT (mass) matrix
    -----------------------------------
    */
    DynamicMatrix <double> mel = mel0;

    if (lumped == 0) {
        mel = ((rho * elArea) / 36) * consistent_mel;
    } else if (lumped == 1) {
        mel = ((rho * elArea) / 4) * blaze::IdentityMatrix<double>(8);
    }
   //loop for/over the total number of elements (considering one element at the time)
   for (int ielm = 0; ielm < nel; ++ielm) {
    
        for (int im = 0; im < nnel; ++im){ //extract nodes for the (iel)-th element
            ndm[im] = nodes(ielm,im);
        }
 
    indexm = feeldofm(ndm,nnel,ndof);
    M=feasmbl1_m(M,mel,indexm); //assemble element matrices in the global one  
   }


   
   //Initial acceleration
    blaze::DynamicVector<double> a0(sdof);
    
    DynamicVector<double,columnVector> temp = -K * d;
    a = inv(M) * temp;
    
    //a = tem * blaze::DynamicVector<double>(-K*d);   // compute a using matrix inversion and matrix-vector multiplication
    //std::cout << "a \n" << a << std::endl; ????? 小精度的都不太对
    
    int iii = 1;
    for (int jj = 0; jj < gcoord.rows();++jj) {
        if (gcoord(jj, 0) == 0.0) {
            a[iii-1] = 0.0;
            a[iii] = 0.0;
        }
        iii += 2;
    }
 
    // -------------------- NEWMARK PREDICTOR-CORRECTOR SOLUTION -----------------------------

    // set initial d,v,a for the whole process
    DynamicMatrix <double> U_dN(sdof, nt+1, 0.0); // matrix where displacement solutions at each time step are stored
    DynamicMatrix <double> U_vN(sdof, nt+1, 0.0); //matrix where velocity solutions at each time step are stored
    DynamicMatrix <double> U_aN(sdof, nt+1, 0.0); //matrix where acceleration solutions at each time step are stored
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
   
    blaze::DynamicVector<double> d1N(sdof);
    blaze::DynamicVector<double,columnVector> v1N(sdof);
    blaze::DynamicVector<double,columnVector> b;
    blaze::DynamicMatrix<double,columnMajor> A;
    blaze::DynamicMatrix<double,columnMajor> L, U, P;

    blaze::DynamicVector<double,columnVector > z(sdof);
    blaze::DynamicVector<double,columnVector > a1N(sdof,1);
    blaze::DynamicVector<double,columnVector > a1N_int(96);
   

    blaze::DynamicVector<double> d1N_touse(rowsize);
    blaze::DynamicMatrix<double> K_touse(rowsize,rowsize);
    blaze::DynamicMatrix<double> fd_touse(rowsize, nt+1);
  
    blaze::DynamicMatrix<double,columnMajor> A_touse(rowsize, rowsize);
    

    // calculate
    A = M + beta_b *(dt*dt)*K;
        
    A_touse = rows(columns(A, to_use),to_use); 
    lu(A_touse, L, U, P);

    int n = 0;
    while(n<nt){
        // PREDICTOR PHASE
        d1N = d + dt*v + ((dt*dt)/2) * (1-2*beta_b)*a;
        v1N = v + (1-gamma_b)*dt*a;
            
            
        //SOLUTION
        K_touse = rows(columns(K, to_use),to_use);  // rows, columns, num_rows, num_cols
            
        //reset(d1p_touse);
        d1N_touse = elements(d1N, to_use);

        //clear(fd1_touse);
        fd_touse = rows(fd,to_use);
            
        b = column(fd_touse, n+1) - K_touse * d1N_touse ;


        //LU_decomposition
        z = inv(L) * b; 
        a1N_int = inv(U) * z ;

        // reconstruct the global acceleration a1N   
        reset(a1N);
        for (int i = 0; i < to_use.size(); ++i) {
            a1N[to_use[i]] = a1N_int[i];
        }
        for (int j = 0; j< to_cancel.size();++j){
            a1N[to_cancel[j]] = 0.0;
        }
        //elements(a1,to_cancel) = 0.0;

        // CORRECTOR PHASE
        d1N = d1N + beta_b * (dt*dt)*a1N;
        v1N = v1N+(1-gamma_b)*dt*a1N;
            
            /**
             * SUBSTITUTING
                U_d(:,n+1)=d1;
                U_v(:,n+1)=v1;
                U_a(:,n+1)=a1;
             */
        for( size_t i=0UL; i<sdof; i++ ) {     
            U_dN(i,n+1) = d1N[i];
            U_vN(i,n+1) = v1N[i];
            U_aN(i,n+1) = a1N[i];

        }
        // column(U_d, n+1) = d1;
        // column(U_v, n+1) = v1; 
        // column(U_a, n+1) = a1;

        d = d1N;
        v = v1N;
        a = a1N;
        
        n++;     
    } 

    // for plotting
    std::ofstream fout0("2D_Newmark_U_dN_plot_res_106.csv");
    fout0 << "dt,U_dN\n";
    for (std::size_t step = 0; step <= nt; step++){
        fout0 << dt*step<< ","
            <<U_dN(106,step) << "\n";
    }
    fout0.close();

    std::ofstream fout2("2D_Newmark_U_dN_plot_res_52.csv");
    fout2 << "dt,U_dN\n";
    for (std::size_t step = 0; step <= nt; step++){
        fout2 << dt*step<< ","
            <<U_dN(52,step) << "\n";
    }
    fout2.close();


    return 0;  
}
