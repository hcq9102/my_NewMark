#include <iostream>
#include <cmath>
#include <blaze/Blaze.h>
#include <blaze/Math.h>
#include <blaze/math/Submatrix.h>
#include <blaze/Forward.h>
#include <blaze/math/Columns.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Subvector.h>

#include "preprocess.hpp"

#include <fstream>
#include <string>
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


    for (int iel = 0; iel < nel; ++iel) { // loop over the total number of elements
        
        for (int i = 0; i < nnel; ++i) {
            nd[i] = nodes(iel, i); // extract nodes for the (iel)-th element
            xcoord[i] = gcoord(nd[i]-1, 0); // extract x value of the node
            ycoord[i] = gcoord(nd[i]-1, 1); // extract y value of the node
        }

        DynamicMatrix <double> kel(edof,edof);      //initialization of element (stiffness) matrix
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
                //std::cout << "kel before "<< kel <<endl;
                //std::cout << "trans(kinmtx2) * matmtx * kinmtx2 * (wtx*wty *detjacob) : " << trans(kinmtx2) * matmtx * kinmtx2 * (wtx*wty *detjacob) << std::endl;
                //kel +=  trans(kinmtx2) * matmtx * kinmtx2 * eval(wtx * wty * detjacob);   //element stiffness matrix
                kel += (trans(kinmtx2) * matmtx) * kinmtx2 * (wtx*wty *detjacob);
                std::cout << "kel after "<< kel <<endl;
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
   
   //loop for/over the total number of elements (considering one element at the time)
   for (int ielm = 0; ielm < nel; ++ielm) {

       for (int im = 0; im < nnel; ++im){ //extract nodes for the (iel)-th element
            ndm[im] = nodes(ielm,im);
       }

    /*-----------------------------------
    ELEMENT (mass) matrix
    -----------------------------------
    */
       if (lumped == 0) {
            mel = ((rho * elArea) / 36) * consistent_mel;
        } else if (lumped == 1) {
            mel = ((rho * elArea) / 4) * unitmatrix;
        }

        indexm = feeldofm(ndm,nnel,ndof); 
        M=feasmbl1_m(M,mel,indexm); //assemble element matrices in the global one
    
   }

   
    
    // for (int i =0; i< M.rows(); i++){
    //     for (int j =0; j< M.columns(); j++){
    //         if(M(i,j) != 0){
    //             std::cout <<M(i,j) <<std::endl;
    //         }    
    //     }
    // }

   //Initial acceleration
    blaze::DynamicVector<double> a0(sdof);
    DynamicVector<double,columnVector> temp = -K * d;
    a = inv(M) * temp;

    int iii = 0;
    for (int jj = 0; jj < gcoord.rows();++jj) {
        if (gcoord(jj, 0) == 0.0) {
            a[iii] = 0.0;
            a[iii + 1] = 0.0;
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
    DynamicMatrix <double> U_d(sdof, nt+1); // matrix where displacement solutions at each time step are stored
    DynamicMatrix <double> U_v(sdof, nt+1); //matrix where velocity solutions at each time step are stored
    DynamicMatrix <double> U_a(sdof, nt+1); //matrix where acceleration solutions at each time step are stored
    // initialize first column
    for( size_t i=0UL; i< sdof; ++i ) {   
        U_d(i,0) = d[i];
        U_v(i,0) = v[i];
        U_a(i,0) = a0[i];
        
    }   
    /*
    ---------------------------------------------------------------
    Iniziatise WR time-discrete matrix (initial WRs=initial conditions)
    ---------------------------------------------------------------
    */ 
    DynamicMatrix <double> WR(sdof, nt+1);
    for( int qq=0; qq<sdof; qq++ ) {
        for( int j=0; j<nt+1; j++ ) {
            WR(qq,j)= d[qq];
        }
    }
    /*
    ---------------------------------------------------------------
    initialize first row of WR_STOR=ITERATION 0 for all dofs
    ---------------------------------------------------------------
    */ 
    DynamicMatrix <double> WR_STOR_1(WR);
    DynamicMatrix<double> WR_STOR_2; 

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
    DynamicVector<double> d1p(sdof);
    DynamicVector<double,columnVector> v1p(sdof);
    DynamicVector<double,columnVector> b;
    DynamicMatrix<double,columnMajor> A;
    DynamicMatrix<double,columnMajor> L, U, P;

    DynamicVector<double,columnVector > z(sdof);
    DynamicVector<double,columnVector > a1(sdof);
    DynamicVector<double,columnVector > a_int(96);
   
    DynamicVector<double,columnVector > d1(sdof);
    DynamicVector<double,columnVector > v1(sdof);
    DynamicVector<double> v0(sdof);
    DynamicMatrix <double> cc(1UL, nt+1);
    DynamicMatrix <double> dd(1UL, nt+1);

    while((e>pow(10.0,-14))){
        
        //Initial conditions (reimposed at the beginning of every cycle)
        d = d_0;
        v = v0; 
        a = a0;

        //Force vector
        blaze::DynamicMatrix <double> fd1(sdof,nt+1);     //forces matrices (over dofs and time)
        fd1 =  K_minus_j * WR;
        /*
        -----------------------------------
        Solution (time steps with Newmark)
        ----------------------------------- 
        */
        A = M + beta_b *(dt*dt)*K_plus_j;
        
        exit(1);
        blaze::DynamicMatrix<double,columnMajor> A_touse(to_use.size(), to_use.size());
        for (int i = 0; i < to_use.size(); i++) {
            for (int j = 0; j < to_use.size(); j++) {
                A_touse(i, j) = A(to_use[i], to_use[j]);
            }
        }
        lu(A_touse, L, U, P);
        int n = 0;
        while(n<nt){
            // PREDICTOR PHASE
            d1p = d + dt*v + ((dt*dt)/2) * (1-2*beta_b)*a;

            std::cout << " I am L245" <<std::endl;
            v1p = v + (1-gamma_b)*dt*a;

            std::cout << " I am L249" <<std::endl;
            //SOLUTION
            // b=fd1(to_use,n+1)-Kplus(to_use,to_use)*d1p(to_use);
            blaze::DynamicMatrix<double> K_plus_j_touse(to_use.size(), to_use.size());
            for (int i = 0; i < to_use.size(); i++) {
                for (int j = 0; j < to_use.size(); j++) {
                    K_plus_j_touse(i, j) =  K_plus_j(to_use[i], to_use[j]);
                }
            }

            blaze::DynamicVector<double> d1p_touse(to_use.size());
            for (int i = 0; i < to_use.size(); i++) {              
                     d1p_touse[i] = d1p[to_use[i]];
            }

            blaze::DynamicMatrix<double> fd1_touse(to_use.size(), to_use.size());
            for (int i = 0; i < to_use.size(); i++) {
                for (int j = 0; j < to_use.size(); j++) {
                    fd1_touse(i, j) =  fd1(to_use[i], to_use[j]);
                }
            }

            b = column(fd1_touse, n+1) - K_plus_j_touse * d1p_touse ;
            
            std::cout<< "u" <<U(0,0) << std::endl;
            std::cout<< "u" <<U(2,2) << std::endl;
            std::cout<< "u" <<U(2,3) << std::endl;
            std::cout<< "u" <<U(22,22) << std::endl;
            std::cout<< "u" <<U(95,95) << std::endl;
            std::cout<< "u" <<U.rows() << std::endl;
            std::cout<< "u" <<U.columns() << std::endl;

            //LU_decomposition
            z = inv(L) * b;

            a_int = inv(U) * z ;
            std::cout << " I am L273" <<std::endl;
            //a1 = zeros(length(to_cancel)+length(to_use),1);
            DynamicVector<double,columnVector > a1(to_cancel.size() + to_use.size());

            for (int o=0; o < to_use.size(); o++) {
                a1[to_use[o]] = a_int[o];
            }
            for (int o=0; o < to_cancel.size(); o++) {
                a1[to_cancel[o]] = 0;
            }
            std::cout << " I am L282" <<std::endl;
            // CORRECTOR PHASE
            d1 = d1p + beta_b * (dt*dt)*a1;
            v1 = v1p+(1-gamma_b)*dt*a1;
            std::cout << " I am L286" <<std::endl;
            /**
             * SUBSTITUTING
                U_d(:,n+1)=d1;
                U_v(:,n+1)=v1;
                U_a(:,n+1)=a1;
             */
            for( size_t i=0UL; i<sdof; ++i ) {     
                U_d(i,n+1) = d1[i];
                U_v(i,n+1) = v1[i];
                U_a(i,n+1) = a1[i];

            } 
            d = d1;
            v = v1;
            a = a1;

            for( size_t i=0UL; i<sdof; ++i ) {     
                WR(i,n+1) = d1[i];
            }
            n++;  
        }       

        i = i+1;

        if( i % 2 == 0){ //iteration number is even
            e = 1;
        }else{ //iteration number is odd 
            WR_STOR_2 = WR;
            blaze::DynamicMatrix<double> e_t(nt + 1, 1);
            for (int n1 = 0; n1 < nt + 1; ++n1) {
                auto tmp = blaze::column(WR_STOR_2, n1) - blaze::column(WR_STOR_1,n1);
                e_t(n1, 1) = blaze::l2Norm(tmp);
            }
            dd = abs(cc);
            //std::cout << "dd " << setprecision(16) << dd(0,8000) <<endl;
            e = max(dd);
            WR_STOR_1=WR;
        }

    }  
    
    std::ofstream fout0("2D_Jacobi_WR_plot_res.csv");
    fout0 << "dt,WR\n";
    for (std::size_t step = 0; step <= nt; step++){
        fout0 << dt*step<< ","
             <<WR(108,step) << "\n";
    }
    fout0.close();     

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

