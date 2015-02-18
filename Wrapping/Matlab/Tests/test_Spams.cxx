/**
 *       @file  test_Spams.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "02-22-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkSpamsWeightedLassoSolver.h"
#include "linalg.h"
#include "decomp.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{

  // vnl_matrix<double> A;
  // vnl_vector<double> b, w, x;
  
  // const double array_A [] =
  //   {
  //   1,  1, 
  //   -1, 2,
  //   2,  1
  //   };
  // const double array_b [] =
  //   {
  //   2, 2, 3
  //   };

  // A = vnl_matrix<double>(array_A, 3, 2);
  // utl::PrintVnlMatrix(A, "A");
  // b = vnl_vector<double>(array_b, 3);
  // utl::PrintVnlVector(b, "b");

  // w.set_size(2);
  // w.fill(5);
  // w[0]=2;
  
  // std::cout << "weighted lasso using spams" << std::endl << std::flush;
  // typedef itk::SpamsWeightedLassoSolver<double> SoverType;
  // SoverType::Pointer solver = SoverType::New();
  // solver->SetA(A);
  // solver->Setb(b);
  // solver->Setw(w);
  // solver->Solve();
  // solver->Print(std::cout<<"solver = ");
  // x = solver->Getx();
  // std::cout << "x = " << x << std::endl << std::flush;
  // std::cout << "cost = " << solver->EvaluateCostFunction() << std::endl << std::flush;
  
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
  // spams::lasso<double>(X,D,alpha,2,lambda);
  // alpha.print("alpha");
  
  // lasso (non-negativity)
  spams::constraint_type mode = spams::PENALTY;
  // spams::lasso<double>(X,D,alpha,2,lambda,0.0, mode, true);
  // alpha.print("alpha");
}
