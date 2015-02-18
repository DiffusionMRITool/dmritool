/**
 *       @file  itkSpamsWeightedLassoSolverTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "02-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */
#include "itkSpamsWeightedLassoSolver.h"

typedef utl::NDArray<double,2> MatrixType;
typedef utl::NDArray<double,1> VectorType;

/**
 * \brief  test itkSpamsWeightedLassoSolver
 */
int 
main (int argc, char const* argv[])
{
  MatrixType A;
  VectorType b, w, x;
  
  const double array_A [] =
    {
    1,  1, 
    -1, 2,
    2,  1
    };
  const double array_b [] =
    {
    2, 2, 3
    };

  A = MatrixType(array_A, 3, 2);
  utl::PrintUtlMatrix(A, "A");
  b = VectorType(array_b, 3);
  utl::PrintUtlVector(b, "b");

  w.ReSize(2);
  w.Fill(5);
  w[0]=2;
  
  std::cout << "weighted lasso using spams" << std::endl << std::flush;
  typedef itk::SpamsWeightedLassoSolver<double> SolverType;
  SolverType::Pointer solver = SolverType::New();
  solver->SetA(SolverType::MatrixPointer(new SolverType::MatrixType(A)));
  solver->Setb(SolverType::VectorPointer(new SolverType::VectorType(b)));
  solver->Setw(SolverType::VectorPointer(new SolverType::VectorType(w)));
  solver->Solve();
  solver->Print(std::cout<<"solver = ");
  x = solver->Getx();
  std::cout << "x = " << x << std::endl << std::flush;
  std::cout << "cost = " << solver->EvaluateCostFunction() << std::endl << std::flush;

  return 0;
}
