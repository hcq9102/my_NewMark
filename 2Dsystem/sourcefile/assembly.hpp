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

//compute ELEMENT MATRICES and vectors and ASSEMBLY


// get element matrix
blaze::DynamicMatrix <double> geGlobal_k(const int ngx, const int ngy, const int iop, const int emodul, const double poisso, const int ne, const int nne, const int ndo){
    /*
    -----------------------------------------
     STIFFNESS MATRIX
    -------------------------------------------
    */
    feglqd2(ngx,ngy);                       // [point2,weight2]=feglqd2(nglx,ngly);  
    fematiso(iop,emodul,poisso);           //constitutive matrix matmtx=fematiso(iopt,emodule,poisson);

    double x = 0;
    double y = 0;
    double wtx =0.0;
    double wty =0.0;
    double detjacob=0.0;
    //int count = 0;

    for (int iel = 0; iel < ne; ++iel) { // loop over the total number of elements
        blaze::DynamicMatrix <double> kel = kel0;      //initialization of element (stiffness) matrix as all zeros
        for (int i = 0; i < nne; ++i) {
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
        indexk = feeldofk(nd,nne,ndo);                     
        //assemble element matrices in the global one
        K=feasmbl1_m(K,kel,indexk); //K=feasmbl1(K,kel,index);  
  
    }
    return K;
}

blaze::DynamicMatrix <double> geGlobal_m(const int lump, const int ne, const int nne, const int ndo){
    // -----------------------------------
    // ELEMENT (mass) matrix
    // -----------------------------------
   
    blaze::DynamicMatrix <double> mel = mel0;

    if (lump == 0) {
        mel = ((rho * elArea) / 36) * consistent_mel;
    } else if (lump == 1) {
        mel = ((rho * elArea) / 4) * blaze::IdentityMatrix<double>(8);
    }
   //loop for/over the total number of elements (considering one element at the time)
   for (int ielm = 0; ielm < ne; ++ielm) {
    
        for (int im = 0; im < nne; ++im){ //extract nodes for the (iel)-th element
            ndm[im] = nodes(ielm,im);
        }
 
    indexm = feeldofm(ndm,nnel,ndo);
    M=feasmbl1_m(M,mel,indexm); //assemble element matrices in the global one  
   }

   return M;
}

