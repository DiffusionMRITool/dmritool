/**
 *       @file  spamsLassoTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  01-02-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */  

#include "linalg.h"
#include "decomp.h"


int main(int argc, char** argv) 
{


  spams::Matrix<double> D(3,2);
  D(0,0) = 1;   D(0,1) = 1;
  D(1,0) = -1;  D(1,1) = 2;
  D(2,0) = 2;   D(2,1) = 1;
  D.print("D");

  double lambda = 1.0;

  spams::Matrix<double> X(3,1);
  X(0,0)=2, X(1,0)=2, X(2,0)=-3;
  spams::SpMatrix<double> alpha;

  // lasso
  spams::lasso2<double>(X,D,alpha,2,lambda);
  alpha.print("alpha");
  
  // lasso (Cholesky-based)
  spams::lasso<double>(X,D,alpha,2,lambda);
  alpha.print("alpha");
  
  // lasso (non-negativity)
  spams::constraint_type mode = spams::PENALTY;
  spams::lasso<double>(X,D,alpha,2,lambda,0.0, mode, true);
  alpha.print("alpha");

}
