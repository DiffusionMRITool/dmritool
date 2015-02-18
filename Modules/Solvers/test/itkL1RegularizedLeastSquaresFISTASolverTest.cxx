/**
 *       @file  itkL1RegularizedLeastSquaresFISTASolverTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-8-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkL1RegularizedLeastSquaresFISTASolver.h"
#include "vnl/vnl_matrix.h"
#include "utl.h"


/**
 * \brief  test itkL1RegularizedLeastSquaresFISTASolver
 */
int 
main (int argc, char const* argv[])
{
  typedef utl::NDArray<double,2> MatrixType;
  typedef utl_shared_ptr<MatrixType> MatrixPointer;
  typedef utl::NDArray<double,1> VectorType;
  typedef utl_shared_ptr<VectorType> VectorPointer;

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


  typedef itk::L1RegularizedLeastSquaresFISTASolver<double> FISTASolverType;
  FISTASolverType::Pointer solver = FISTASolverType::New();

  VectorType x0(2);
  x0[0]=-500; x0[1]=500;
  
  solver->SetUseL2SolverForInitialization(false);
  solver->SetA(MatrixPointer(new MatrixType(A)));
  solver->Setb(VectorPointer(new VectorType(b)));
  solver->Setw(VectorPointer(new VectorType(w)));
  solver->SetDebug(true);
  solver->SetMinRelativeChangeOfCostFunction(0.001);
  solver->Solve(x0);
  solver->Print(std::cout<<"solver = ");
  x = solver->Getx();
  std::cout << "x = " << x << std::endl << std::flush;

  std::cout << "least square for initialization" << std::endl << std::flush;
  solver->SetUseL2SolverForInitialization(true);
  solver->SetA(MatrixPointer(new MatrixType(A)));
  solver->Setb(VectorPointer(new VectorType(b)));
  solver->Setw(VectorPointer(new VectorType(w)));
  solver->Solve();
  solver->Print(std::cout<<"solver = ");
  x = solver->Getx();
  std::cout << "x = " << x << std::endl << std::flush;
  std::cout << "cost = " << solver->EvaluateCostFunction() << std::endl << std::flush;

  
  return 0;
}
