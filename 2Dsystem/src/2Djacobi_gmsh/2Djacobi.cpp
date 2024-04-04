
#include <iostream>
#include <cmath>
#include <hpx/hpx_main.hpp>  // Need main source file

#include <blaze/Blaze.h>
#include <blaze/Math.h>
#include <blaze/math/Submatrix.h>
#include <blaze/Forward.h>
#include <blaze/math/Columns.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Subvector.h>
#include <blaze/math/Elements.h>

#include <yaml-cpp/yaml.h>

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
    // yaml start
    const std::string  inputFilename = input.yaml;

    YAML::Node config = YAML::LoadFile(inputFilename);

     // read file
    if (config["Mesh"]["File"])
        std::string filename = config["Mesh"]["File"].as<std::string>();
    else {
        std::cerr << "Error: Please specify the filename.\n";
        exit(1);
    }

    // read dimension
    if (config["Mesh"]["Dimension"])
        d_modelDeck_p->d_dim =config["Model"]["Dimension"].as<size_t>();
    else {
        std::cerr << "Error: Please specify the Dimension.\n";
        exit(1);
    }
    // yaml ends

    readGmsh("testGmsh.msh", gcoord, nodes, dim,element_type);
    
    
    //compute ELEMENT MATRICES and vectors and ASSEMBLY
    /*
    -----------------------------------------
     STIFFNESS MATRIX  && MASS MATRIX
    -----------------------------------------
    */
    geGlobal_k(nglx, ngly, iopt, emodule, poisson, nel, nnel, ndof);
    geGlobal_m(lumped, nel, nnel, ndof);

   //Initial acceleration
    
    a = blaze::solve(M, -K * d);


    // first col, initial a should be 0; 
    int iii = 1;
    for (int jj = 0; jj < gcoord.rows();++jj) {
        if (gcoord(jj, 0) == 0.0) {
            a[iii-1] = 0.0;
            a[iii] = 0.0;
        }
        iii += 2; // change here, due to the nodes are disorder
    }
    
    blaze::DynamicVector<double> a0 = a;

     /*
     -------------------------------------------------------------------------
     WAVEFORM RELAXATION SOLUTION (JACOBI)
     -------------------------------------------------------------------------
     */   

    // 1. Partitioning
    // for jacobi:... get K, M _plus and minus matrices
    DynamicMatrix<double> M_plus_j, K_plus_j, K_minus_j;
    M_plus_j = diagMatrix(M);
    K_plus_j = diagMatrix(K);
    K_minus_j = K_plus_j - K;

    // 2. Initialization
   
    blaze::DynamicVector<double> d_0, v_0, a_0; 
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
    DynamicMatrix <double> WR(sdof, nt+1, 0.0);
  
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

  
    DynamicMatrix<double,columnMajor> A_touse(rowsize, rowsize);
    DynamicMatrix <double> fd1;
    DynamicVector <double> dd(nt+1,0.0);
    DynamicVector<double> e_t(nt+1,0.0);

    double corrector_factor1 = beta_b * (dt*dt);
    double corrector_factor2 = (1-gamma_b)*dt;
    double corrector_factor3 = ((dt*dt)/2) * (1-2*beta_b);

    // Solution (time steps with Newmark)
    A = M + corrector_factor1*K_plus_j;  
    A_touse = rows(columns(A, to_use),to_use); 
    lu(A_touse, L, U, P);

    invert(L);
    invert(U);
    
    while((e > pow(10.0,-14))){
        //Initial conditions (reimposed at the beginning of every cycle)
        d = d_0;
        v = v0; 
        a = a0;

        //Force vector (over dofs and time)
        clear(fd1); 
        fd1 =  K_minus_j * WR;      
        
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
            // reconstruct the global acceleration
            reset(a1);
            
            elements(a1,to_use) = U*L * (column(rows(fd1,to_use), n+1) - rows(columns(K_plus_j, to_use),to_use) * elements(d1p, to_use));
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
            //cout << "n= " << n << endl;  
        }  
        
        i = i+1;

        if( i % 2 == 0){ //iteration number is even
            e = 1;
        }else{ //iteration number is odd 
            WR_STOR_2 = WR;
            auto diff = WR_STOR_2 - WR_STOR_1;
            for (size_t n1 = 0; n1 < nt + 1; ++n1) {
                //e_t(n1, 0) = blaze::l2Norm(blaze::column(WR_STOR_2, n1) - blaze::column(WR_STOR_1, n1));
                e_t[n1] = blaze::l2Norm(blaze::column(diff, n1));
            }
            
            dd = abs(e_t);
            e = max(dd);
            std::cout<<"i="<< i << std::endl;
            std::cout<<"e \n" << e<<std::endl;
            WR_STOR_1=WR;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "runtime " << elapsed.count() << " s\n";
    
    //for plotting

    // std::ofstream fout1("2D_Jaco_WR_newMesh_node1.csv");
    // fout1 << "dt,WR\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout1 << dt*step<< ","
    //          <<WR(2,step) << "\n";
    // }
    // fout1.close();

    // // std::ofstream fout2("2D_Jaco_(U_d)_newMesh_node1.csv");
    // // fout2 << "dt,U_dN\n";
    // // for (std::size_t step = 0; step <= nt; step++){
    // //     fout2 << dt*step<< ","
    // //         <<U_d(2,step) << "\n";
    // // }
    // // fout2.close();

    // std::ofstream fout2("2D_Jaco_WR_newMesh_node2.csv");
    // fout2 << "dt,WR\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout2 << dt*step<< ","
    //          <<WR(4,step) << "\n";
    // }
    // fout2.close();

    // std::ofstream fout3("2D_Jaco_WR_newMesh_node12.csv");
    // fout3 << "dt,WR\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout3 << dt*step<< ","
    //          <<WR(24,step) << "\n";
    // }
    // fout3.close();

    // std::ofstream fout4("2D_Jaco_WR_newMesh_node13.csv");
    // fout4 << "dt,WR\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout4 << dt*step<< ","
    //          <<WR(26,step) << "\n";
    // }
    // fout4.close();

    // // std::ofstream fout4("2D_Jaco_(U_d)_newMesh_node13.csv");
    // // fout4 << "dt,U_dN\n";
    // // for (std::size_t step = 0; step <= nt; step++){
    // //     fout4 << dt*step<< ","
    // //         <<U_d(27,step) << "\n";
    // // }
    // // fout4.close();


    // std::ofstream fout5("2D_Jaco_WR_newMesh_node14.csv");
    // fout5 << "dt,WR\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout5 << dt*step<< ","
    //          <<WR(28,step) << "\n";
    // }
    // fout5.close();

    // // std::ofstream fout2("2D_Jaco_(U_d)_newMesh_node14.csv");
    // // fout2 << "dt,U_dN\n";
    // // for (std::size_t step = 0; step <= nt; step++){
    // //     fout2 << dt*step<< ","
    // //         <<U_d(14,step) << "\n";
    // // }
    // // fout2.close();


    // std::ofstream fout6("2D_Jaco_WR_newMesh_node15.csv");
    // fout6 << "dt,WR\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout6 << dt*step<< ","
    //          <<WR(30,step) << "\n";
    // }
    // fout6.close();

    // // std::ofstream fout8("2D_Jaco_(U_d)_newMesh_node15.csv");
    // // fout8 << "dt,U_dN\n";
    // // for (std::size_t step = 0; step <= nt; step++){
    // //     fout8 << dt*step<< ","
    // //         <<U_d(15,step) << "\n";
    // // }
    // // fout8.close();

    return EXIT_SUCCESS;  
}

