/*
print all values you want use ./exec >>log.txt
here : /huanqiu@rostam1:/work/chuanqiu/Dynamic_project/2Dsystem/2Dj_g_n$ ./jacogdb >>log.txt
*/

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
    
    
    a0 = a;
    
    //fd(:,1);
    /*
    -------------------------------------------------------------------------
    NEWMARK PREDICTOR-CORRECTOR SOLUTION
    -------------------------------------------------------------------------
    */


     /*
     -------------------------------------------------------------------------
     WAVEFORM RELAXATION SOLUTION (JACOBI)
     -------------------------------------------------------------------------
     */   

    //1. Partitioning
    // for jacobi:... get K, M _plus and minus matrices


    DynamicMatrix<double> M_plus_j, K_plus_j, K_minus_j;
    M_plus_j = diagMatrix(M);
    K_plus_j = diagMatrix(K);
    K_minus_j = K_plus_j - K;


    //2.Initialization
    DynamicVector<double,columnVector > d_0, v_0, a_0;  
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
    DynamicMatrix <double> WR(sdof, nt+1);
    // for( int qq=0; qq<sdof; qq++ ) {
    //     for( int j=0; j<nt+1; j++ ) {
    //         WR(qq,j)= d[qq];
    //     }
    // }
    for (int qq = 0; qq < sdof; ++qq)
    {
        blaze::submatrix(WR, qq, 0, 1, nt+1) = d[qq];
    }
  
    /*
    ---------------------------------------------------------------
    initialize first row of WR_STOR=ITERATION 0 for all dofs
    ---------------------------------------------------------------
    */ 
    DynamicMatrix <double> WR_STOR_1(WR);
    DynamicMatrix<double> WR_STOR_2(sdof, nt+1);
     

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
   
    DynamicVector<double> d1p(sdof);
    DynamicVector<double,columnVector> v1p(sdof);
    DynamicVector<double,columnVector> b;
    DynamicMatrix<double,columnMajor> A;
    DynamicMatrix<double,columnMajor> L, U, P;

    DynamicVector<double,columnVector > z(sdof);
    DynamicVector<double,columnVector > a1(sdof,1);
    DynamicVector<double,columnVector > a_int(96);
   
    DynamicVector<double,columnVector > d1(sdof);
    DynamicVector<double,columnVector > v1(sdof);
    DynamicVector<double> v0(sdof);
    
    DynamicMatrix <double> dd( nt+1,1UL);

    blaze::DynamicVector<double> d1p_touse(rowsize);
    blaze::DynamicMatrix<double> K_plus_j_touse(rowsize,rowsize);
    blaze::DynamicMatrix<double> fd1_touse(rowsize, nt+1);
  
    blaze::DynamicMatrix<double,columnMajor> A_touse(rowsize, rowsize);
    blaze::DynamicMatrix <double> fd1(fd1_0);
    blaze::DynamicMatrix<double> e_t(nt + 1, 1);

    while((e > pow(10.0,-14))){
        
        //Initial conditions (reimposed at the beginning of every cycle)
        d = d_0;
        v = v0; 
        a = a0;

        //Force vector
        //blaze::DynamicMatrix <double> fd1(fd1_0);     //forces matrices (over dofs and time)
        clear(fd1); 
        fd1 =  K_minus_j * WR;      
        
        /*
        -----------------------------------
        Solution (time steps with Newmark)
        ----------------------------------- 
        */
       // move out of while??????????????????????????
        A = M + beta_b *(dt*dt)*K_plus_j;
        
        A_touse = rows(columns(A, to_use),to_use); // correct here
        //clear(A_touse);
        //print_X_value(A_touse);
        lu(A_touse, L, U, P);
        

        int n = 0;
        //blaze::DynamicMatrix<double> K_plus_j_touse(to_use.size(), to_use.size());
        while(n<nt){
            // PREDICTOR PHASE
            d1p = d + dt*v + ((dt*dt)/2) * (1-2*beta_b)*a;
            
            v1p = v + (1-gamma_b)*dt*a;
            
            
            //SOLUTION
            //clear(K_plus_j_touse);
            K_plus_j_touse = rows(columns(K_plus_j, to_use),to_use);  // rows, columns, num_rows, num_cols
            
            //reset(d1p_touse);
            d1p_touse = elements(d1p, to_use);

            //clear(fd1_touse);
            fd1_touse = rows(fd1,to_use);
            
            b = column(fd1_touse, n+1) - K_plus_j_touse * d1p_touse ;


            //LU_decomposition
            z = inv(L) * b;
            
            a_int = inv(U) * z ;
           
            reset(a1);
            for (int i = 0; i < to_use.size(); ++i) {
                a1[to_use[i]] = a_int[i];
            }
            for (int j = 0; j< to_cancel.size();++j){
                a1[to_cancel[j]] = 0.0;
            }
            //elements(a1,to_cancel) = 0.0;

            // CORRECTOR PHASE
            d1 = d1p + beta_b * (dt*dt)*a1;
            v1 = v1p+(1-gamma_b)*dt*a1;
            
            /**
             * SUBSTITUTING
                U_d(:,n+1)=d1;
                U_v(:,n+1)=v1;
                U_a(:,n+1)=a1;
             */
            for( size_t i=0UL; i<sdof; i++ ) {     
                U_d(i,n+1) = d1[i];
                U_v(i,n+1) = v1[i];
                U_a(i,n+1) = a1[i];

            }
            // column(U_d, n+1) = d1;
            // column(U_v, n+1) = v1; 
            // column(U_a, n+1) = a1;

            d = d1;
            v = v1;
            a = a1;

            for( size_t i=0UL; i<sdof; i++ ) {     
                WR(i,n+1) = d1[i];
            }         
            //column(WR, n+1) = d1;
            n++;     
        }  
        
        i = i+1;
        // cout <<"J current i: " << i << endl;
        

        if( i % 2 == 0){ //iteration number is even
            e = 1;
        }else{ //iteration number is odd 
            WR_STOR_2 = WR;
            //clear(e_t);
            for (int n1 = 0; n1 < nt + 1; ++n1) {
                e_t(n1, 0) = blaze::l2Norm(blaze::column(WR_STOR_2, n1) - blaze::column(WR_STOR_1, n1));
            }
            
            dd = abs(e_t);
            e = max(dd);
            std::cout<<"i="<< i << std::endl;
            std::cout<<"e \n" << e<<std::endl;
            WR_STOR_1=WR;
        }
    }
    // for plotting
    std::ofstream fout0("2D_Jaco_WR_plot_res_52.csv");
    fout0 << "dt,WR\n";
    for (std::size_t step = 0; step <= nt; step++){
        fout0 << dt*step<< ","
             <<WR(52,step) << "\n";
    }
    fout0.close();

    return 0;   
}

