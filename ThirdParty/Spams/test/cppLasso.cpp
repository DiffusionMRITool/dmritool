#include "linalg.h"
#include "decomp.h"

int main(int argc, char** argv) {
   int m = 10;
   int p = 20;
   /// external allocation for the matrix
   double* prD = new double[m*p];
   spams::Matrix<double> D(prD,m,p); 
   D.setAleat(); 
   D.normalize();

   /// Allocate a matrix of size m x p
   spams::Matrix<double> D2(m,p); 
   D2.setAleat();
   D2.normalize();

   int n = 100;
   spams::Matrix<double> X(m,n);
   X.setAleat();
   
   /// create empty sparse matrix
   spams::SpMatrix<double> spa;

   spams::lasso2(X,D,spa,10,0.15);  // first simple example
   D.print("D");
   X.print("X");
   spa.print("spa");

   /// extern allocation for the matrix D requires
   /// manually unallocating prD
   delete[](prD);
}
