#pragma once

#include <iostream>
#include <cmath>
#include <blaze/Blaze.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/Forward.h>
#include <blaze/math/Columns.h>
#include <blaze/math/Rows.h>
#include <blaze/Math.h>

#include <fstream>
#include <iomanip>

//-----------------------------------------
// INPUT DATA for parameters
//-----------------------------------------
const size_t nel=40;             //number of elements
const size_t nnel=4;             //number of nodes per element (rectangles=4)
const size_t ndof=2;             //number of degrees of freedom per node (2D=2)
const size_t nnode=54;           //total number of nodes in the system
const size_t sdof=nnode*ndof;    //total number of dofs in the system
const size_t edof=nnel*ndof;     //degrees of freedom per element
const int emodule=1;          //elastic modulus
const size_t E=emodule;          //elastic modulus different notation
double poisson=0.3;        //Poisson's ratio
const int nglx=2, ngly=2;     //2x2 Gauss-Legendre quadrature (sampling points in the two directions x,y)
const int nglxy=nglx*ngly;    //number of sampling points per element

const int iopt=1;             //plane stress analysis
//const size_t iopt=2;             //plane strain analysis

const size_t lumped=1;           //lumped=0 (consistent element mass matrix), =1 (lumped element mass matrix)集总元素质量矩阵

double elArea=0.125;       //area of each element (0.125 is a value specific to the rectangular elements considered in the problem)
const size_t rho=1;              //density of the material

double dt=0.01;            //time step
const size_t nt=5000;            //number of time steps

double beta_b=0.2500;         //Newmark predictor corrector parameter beta
double gamma_b=0.5000;        //Newmark predictor corrector parameter gamma         

/*-----------------------------------------
input data for NODAL COORDINATE values节点坐标值
gcoord(i,j) where i=node n and j = x or y (nodes written such that first
-----------------------------------------
node coordinates are in the first matrix row, x in the first column, y in
the second etc.)
*/
blaze::DynamicMatrix<double,blaze::columnMajor> gcoord{ {0.0, 0.0},{0.5,0.0},{1.0 ,0.0},{1.5, 0.0}, {2.0, 0.0}, {2.5, 0.0}, {3.0, 0.0}, {3.5, 0.0}, {4.0, 0.0},
                                                 {0.0, 0.25}, {0.5,0.25}, {1.0, 0.25}, {1.5,0.25}, {2.0, 0.25}, {2.5, 0.25 },{3.0, 0.25}, {3.5, 0.25}, {4.0, 0.25},
                                                 {0.0, 0.5 }, {0.5,0.5}, {1.0, 0.5}, {1.5,0.5}, {2.0, 0.5}, {2.5, 0.5},{3.0, 0.5}, {3.5, 0.5}, {4.0, 0.5},
                                                 {0.0, 0.75}, {0.5,0.75}, {1.0, 0.75}, {1.5,0.75}, {2.0, 0.75}, {2.5, 0.75},{3.0, 0.75}, {3.5, 0.75}, {4.0, 0.75},
                                                 {0.0, 1.0 }, {0.5,1.0}, {1.0, 1.0}, {1.5,1.0}, {2.0, 1.0}, {2.5, 1.0},{3.0, 1.0}, {3.5, 1.0}, {4.0, 1.0},
                                                 {0.0, 1.25}, {0.5,1.25}, {1.0,1.25}, {1.5,1.25}, {2.0, 1.25}, {2.5, 1.25 },{3.0, 1.25}, {3.5, 1.25}, {4.0, 1.25},
                                              };
/*
-----------------------------------------
 input data for NODAL CONNECTIVITY for each element
 nodes (i,j) where i= element n and j= connected nodes 
-----------------------------------------
 (basically in each row there are the element nodes in a counterclock
 order - 1st row, 1st element - 2nd row, 2nd element
*/

blaze::DynamicMatrix<double> nodes{ { 1,2,11,10 }, { 2,3,12,11 },{ 3,4,13,12 },{ 4,5,14,13 }, { 5,6,15,14 }, { 6,7,16,15 }, { 7,8,17,16 },{ 8,9,18,17 }, 
                        { 10,11,20,19 }, { 11,12,21,20 },{ 12,13,22,21 }, { 13,14,23,22 }, { 14,15,24,23 }, { 15,16,25,24 },{ 16,17,26,25 }, { 17,18,27,26 },
                        { 19,20,29,28 }, { 20,21,30,29 },{ 21,22,31,30 }, { 22,23,32,31 }, { 23,24,33,32 }, { 24,25,34,33 },{ 25,26,35,34 }, { 26,27,36,35 }, 
                        { 28,29,38,37 }, { 29,30,39,38 },{ 30,31,40,39 }, { 31,32,41,40 }, { 32,33,42,41 }, { 33,34,43,42 },{ 34,35,44,43 }, { 35,36,45,44 },
                        { 37,38,47,46 }, { 38,39,48,47 },{ 39,40,49,48 }, { 40,41,50,49 }, { 41,42,51,50 }, { 42,43,52,51 },{ 43,44,53,52 }, { 44,45,54,53 }
                    };

/*-----------------------------------------
 initialization of matrices and vectors
-----------------------------------------
*/

blaze::DynamicVector<double> ff(sdof);           //global system force vector
blaze::DynamicMatrix <double> K(sdof, sdof);     //global stiffness matrix
blaze::DynamicMatrix <double> M(sdof, sdof);     //global mass matrix   


blaze::DynamicVector<double> eldisp(edof);       //element displacement vector
//blaze::DynamicVector<double> d(sdof);            //system displacement vector
blaze::DynamicVector<double> v(sdof);            //system velocity vector
blaze::DynamicVector<double> a(sdof);            //system acceleration vector
blaze::DynamicMatrix <double> fd(sdof,nt+1);     //forces matrices (over dofs and time)

blaze::DynamicMatrix <double> kinmtx(3,edof);     //kinematic matrix
blaze::DynamicMatrix <double> matmtx(3,3);        //constitutive matrix
blaze::DynamicMatrix <double> kinmtx2(3,edof);

// initialize element stiffness/ mass matrix
//k
blaze::DynamicVector<double> nd(nnel); 
blaze::DynamicVector<double> xcoord(nnel);
blaze::DynamicVector<double> ycoord(nnel);    
blaze::DynamicVector<double> indexk(1,edof);           //index vector
blaze::ZeroMatrix<double> kel0(edof,edof);      //initialization of element (stiffness) matrix
//blaze::DynamicMatrix <double> kel(edof,edof);
// m
blaze::DynamicVector<double> ndm(nnel);
blaze::DynamicVector<double> indexm(1,edof); 
blaze::DynamicMatrix <double> mel(edof,edof);       //initialization of element (mass) matrix
blaze::ZeroMatrix<double> mel0(edof,edof);  
blaze::DynamicMatrix <double> consistent_mel{
               {4, 0, 2, 0, 1, 0, 2, 0},
               {0, 4, 0, 2, 0, 1, 0, 2},
               {2, 0, 4, 0, 2, 0, 1, 0},
               {0, 2, 0, 4, 0, 2, 0, 1},
               {1, 0, 2, 0, 4, 0, 2, 0},
               {0, 1, 0, 2, 0, 4, 0, 2},
               {2, 0, 1, 0, 2, 0, 4, 0},
               {0, 2, 0, 1, 0, 2, 0, 4}};

blaze::DynamicMatrix <double> unitmatrix{
               {1, 0, 0, 0, 0, 0, 0, 0},
               {0, 1, 0, 0, 0, 0, 0, 0},
               {0, 0, 1, 0, 0, 0, 0, 0},
               {0, 0, 0, 1, 0, 0, 0, 0},
               {0, 0, 0, 0, 1, 0, 0, 0},
               {0, 0, 0, 0, 0, 1, 0, 0},
               {0, 0, 0, 0, 0, 0, 1, 0},
               {0, 0, 0, 0, 0, 0, 0, 1}};

blaze::ZeroMatrix<double> fd1_0(sdof,nt+1); 
blaze::DynamicMatrix <double> jacobi2(2,2);
blaze::ZeroVector<double> a1_zero(sdof);


//-----------------------------------------
// preprocess 
//-----------------------------------------


//initialization
blaze::DynamicMatrix <double> pointx;
blaze::DynamicMatrix <double> weightx;

blaze::DynamicMatrix <double> pointy;
blaze::DynamicMatrix <double> weighty;


//sampling points & weights
int ngl = std::max(nglx,ngly);
// initialization
blaze::DynamicMatrix <double> point2(ngl,2); 
blaze::DynamicMatrix <double> weight2(ngl,2);

blaze::DynamicVector<double, blaze::rowVector> shape( 4UL );
blaze::DynamicVector<double> dhdx(nnel), dhdy(nnel), dhdr(nnel), dhds(nnel);





//---------------------------------------------------------------------
//Initial imposed displacement 'd'-we pull the points horizontally
blaze::DynamicVector<double> d {0.0, 0.0, 1.5, 0.0, 3.0, 0.0, 4.5, 0.0, 6.0, 0.0, 7.5, 0.0, 9.0, 0.0, 10.5, 0.0, 12.0, 0.0,
                                0.0, 0.0, 1.5, 0.0, 3.0, 0.0, 4.5, 0.0, 6.0, 0.0, 7.5, 0.0, 9.0, 0.0, 10.5, 0.0, 12.0, 0.0,
                                0.0, 0.0, 1.5, 0.0, 3.0, 0.0, 4.5, 0.0, 6.0, 0.0, 7.5, 0.0, 9.0, 0.0, 10.5, 0.0, 12.0, 0.0,
                                0.0, 0.0, 1.5, 0.0, 3.0, 0.0, 4.5, 0.0, 6.0, 0.0, 7.5, 0.0, 9.0, 0.0, 10.5, 0.0, 12.0, 0.0,
                                0.0, 0.0, 1.5, 0.0, 3.0, 0.0, 4.5, 0.0, 6.0, 0.0, 7.5, 0.0, 9.0, 0.0, 10.5, 0.0, 12.0, 0.0,
                                0.0, 0.0, 1.5, 0.0, 3.0, 0.0, 4.5, 0.0, 6.0, 0.0, 7.5, 0.0, 9.0, 0.0, 10.5, 0.0, 12.0, 0.0};

// blaze::DynamicVector<double> to_use{3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
//                                     21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
//                                     39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,
//                                     57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,
//                                     75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,
//                                     93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108};
std::vector<size_t> to_use{2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
                           20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
                           38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
                           56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,
                           74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,
                           92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107};

 //blaze::DynamicVector<double> to_cancel{1,2,19,20,37,38,55,56,73,74,91,92};
std::vector<size_t> to_cancel{0,1,18,19,36,37,54,55,72,73,90,91};

//std::vector<size_t>
//----------------------------------------------------------------------------------


// diagnoal matrice
blaze::DynamicMatrix<double> diagMatrix(blaze::DynamicMatrix<double> &Y){
    blaze::DynamicMatrix<double> D(Y.rows(), Y.columns());
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


//-----------------------------------------------------------------------
blaze::DynamicMatrix <double> feasmbl1_m(blaze::DynamicMatrix <double> &X,blaze::DynamicMatrix <double> &x,blaze::DynamicVector<double, blaze::columnVector> &idx_m){
    int edof = idx_m.size();
    for(int i =0; i<edof; i++){
        int ii = idx_m[i];
        for(int j = 0; j<edof; j++){
            int jj = idx_m[j];
            X(ii,jj) = X(ii,jj) + x(i,j);
        }
    }
    return X;
}


void federiv2(int n,blaze::DynamicVector<double> &r,blaze::DynamicVector<double> &s,blaze::DynamicMatrix <double> &z){
    for(int i  = 0; i< n; i++ ){
        dhdx[i] = z(0,0) * r[i] + z(0,1) * s[i];
        dhdy[i] = z(1,0) * r[i] + z(1,1) * s[i];
    }
} 

void fekine2D(int n,blaze::DynamicVector<double> &dx,blaze::DynamicVector<double> &dy ){
    for(int i  = 0; i< n; i++ ){
        int i1 = 2 *i;
        int i2 = i1 + 1;
        kinmtx2(0,i1) = dx[i];
        kinmtx2(1,i2) = dy[i];
        kinmtx2(2,i1) = dy[i];
        kinmtx2(2,i2) = dx[i];
    }
}

blaze::DynamicVector<double> feeldofk(blaze::DynamicVector<double> &nd_k, int x, int y){
    int edof = x * y;
    blaze::DynamicVector<double> resk(edof);
    int k = 0;

    for(int i = 0; i<x;i++){
        int start = (nd_k[i] -1) * y;
        for(int j =0; j<y;j++){
            resk[k] = start + j;
            k++;
        }
    }
    return resk;
    
}

blaze::DynamicVector<double> feeldofm(blaze::DynamicVector<double> &nd_m, int x, int y){
    int edof = x * y;
    blaze::DynamicVector<double> resm(edof);
    int k = 0;

    for(int i = 0; i<x;i++){
        int start = (nd_m[i] -1) * y;
        for(int j =0; j<y;j++){
            resm[k] = start + j;
            k++;
        }
    }
    return resm;
    
}

//find corresponding integration points and weights
void feglqd1x(int z){
    if(z == 1){
        pointx ={{0.0, 0.0}};
        weightx ={{2.0, 2.0}};
    } 
    else if(z == 2){
        pointx ={{-0.577350269189626, -0.577350269189626},{0.577350269189626,0.577350269189626}};
        weightx ={{1.0, 1.0},{1.0, 1.0}};
    }
}

void feglqd1y(int y){
    if(y == 1){
        pointy ={{0.0, 0.0}};
        weighty ={{2.0, 2.0}};
    } 
    else if(y == 2){
        pointy ={{-0.577350269189626, -0.577350269189626},{0.577350269189626,0.577350269189626}};
        weighty ={{1.0, 1.0},{1.0, 1.0}};
    }
}

void feglqd2(int x, int y){
      
    // find corresponding integration points and weights
    feglqd1x(x);
    feglqd1y(y);
    //quadrature for two-dimension
    for(int i = 0; i<x;i++){
        point2(i,0)= pointx(i,0);
        weight2(i,0) = weightx(i,0);
    }
    for(int j = 0; j<y;j++){
        point2(j,1)= pointy(j,0);
        weight2(j,1) = weighty(j,0);  // pointy(inty) = pointx(j,2) ?
    }
}

template<typename T>
void feisoq4(T x, T y){
    
    // shape functions
    shape[0] = 0.25*(1-x)*(1-y); //N(1)
    shape[1] =0.25*(1+x)*(1-y);  //N(2)
    shape[2] =0.25*(1+x)*(1+y);  //N(3)
    shape[3] =0.25*(1-x)*(1+y);  //N(4)
    
    // derivatives 
    dhdr[0] = -(1-y)/4;
    dhdr[1] = 0.25*(1-y);
    dhdr[2] = 0.25*(1+y);
    dhdr[3] = -0.25*(1+y);

    dhds[0]=-0.25*(1-x);
    dhds[1]=-0.25*(1+x);
    dhds[2]=0.25*(1+x);
    dhds[3]=0.25*(1-x);
}


blaze::DynamicMatrix <double> fejacobi2( int n, blaze::DynamicVector<double> &r,blaze::DynamicVector<double> &s, blaze::DynamicVector<double> &xcrd,blaze::DynamicVector<double> &ycrd ){
    blaze::DynamicMatrix <double> res{{0,0},{0,0}};
    for(int i = 0; i< n; i++){
        res(0,0)=res(0,0)+r[i]*xcrd[i];
        res(0,1)=res(0,1)+r[i]*ycrd[i];
        res(1,0)=res(1,0)+s[i]*xcrd[i];
        res(1,1)=res(1,1)+s[i]*ycrd[i];
    }
    return res;
}


void fematiso(int x, int y, double z){ 
    if(x == 1){
        blaze::StaticMatrix<double,3UL, 3UL> temp1= {{1,z,0},{z,1,0},{0,0,(1-z)/2}};
        matmtx = (y / (1-z*z))* temp1;
    }else if(x == 2){
        blaze::StaticMatrix<double,3UL, 3UL> temp2{{1-z,z,0},{z,1-z,0},{0,0,(1-2*z)/2}};
        matmtx = (y / ((1-2*z)*(1+z)))* temp2;
    }

}

// decompose
// diagnoal matrice
blaze::DynamicMatrix<double> lowerMatrix(blaze::DynamicMatrix<double>& Y){
    blaze::DynamicMatrix<double> X(Y);
    for(size_t i = 0UL; i <Y.rows(); i++){
        for(size_t j = i+1; j< Y.columns(); j++){
            X(i,j) = 0.0;
        }
    }
    return X;
}

void print_X_value(blaze::DynamicMatrix <double> &X){

    std::ofstream file("2DWR_matrix.csv"); // create an output file stream
    if (file.is_open())
    {
        // loop over the rows of the matrix and write them to the file
        for (size_t i = 0; i < X.rows(); ++i)
        {
            // loop over the columns of the matrix and write them to the file
            for (size_t j = 0; j < X.columns(); ++j)
            {
                // write the matrix element to the file
                file << std::fixed << X(i, j);

                // add a comma separator for all but the last column
                if (j < X.columns() - 1)
                    file << ",";
            }
            // add a newline character after each row
            file << "\n";
        }

        file.close(); // close the file
        std::cout << "Matrix output to file:2DWR matrix.txt\n";
    }
    else
    {
        std::cerr << "Error: could not open file for writing\n";
    }

}

void print_vector_value(blaze::DynamicVector<double> &X){

        // 打开输出文件流
        std::ofstream outfile("v1p.csv");

        // 将向量的每个元素输出到文件中，并在每个元素后面添加一个换行符
        for (size_t i = 0; i < X.size(); i++) {
            outfile << X[i] << "\n";
        }

        // 关闭输出文件流
        outfile.close();
        std::cout << "vector output to file:v1p.csv\n";
}
    

