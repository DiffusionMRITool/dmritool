/**
 *       @file  test_SHBasis.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "09-20-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */


#include "itkSHBasisGeneratorTestCLP.h"
#include "utl.h"



/**
 * \brief  test Spherical Harmonic basis 
 */
int 
main (int argc, char const* argv[])
{
  
  // GenerateCLP
  PARSE_ARGS;
  
  typedef float  TScalarType;

  typedef utl::NDArray<TScalarType,2> MatrixType;
  utl_shared_ptr<MatrixType> mat = utl::ReadGrad<TScalarType>(_dataOrientationsFile, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  std::cout << "mat = \n" << *mat << std::endl << std::flush;
  
  utl_shared_ptr<MatrixType> matT3 = utl::ReadGrad<TScalarType>(3, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  std::cout << "mat t3 = \n" << *matT3 << std::endl << std::flush;

  
  utl_shared_ptr<MatrixType> basisMatrix = utl::ComputeSHMatrix(4, *matT3, CARTESIAN_TO_SPHERICAL);

  if (_verbose)
    std::cout << "basisMatrix: \n" << *basisMatrix;

  return 0;
}
