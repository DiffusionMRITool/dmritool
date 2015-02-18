/**
 *       @file  itkDiffusionTensorTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-02-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkDiffusionTensor.h"
#include "utl.h"

#include "utlRotationMatrixFromVectors.h"



/**
 * \brief  test \ref Tensor3x3
 */
int 
main (int argc, char const* argv[])
{
  vnl_matrix<double> matrix(3,3,0.0);

  // for ( int i = 0; i < 3; i += 1 ) 
  //   {
  //   matrix(i,i)= i+1;
  //   }

  for ( int i = 0; i < 3; i += 1 ) 
    {
    for ( int j = 0; j < 3; j += 1 ) 
      {
      matrix(i,j)= 1.0/double(1+i+j);
      }
    }

  std::cout << "matrix = \n" << matrix << std::endl << std::flush;

  itk::DiffusionTensor<double> tensor;
  tensor.SetVnlMatrix(matrix);
  

  itk::DiffusionTensor<double> tensor_inv = tensor.Inv();
  
  std::cout << "tensor = \n" << tensor << std::endl << std::flush;
  std::cout << "tensor_inv = \n" << tensor_inv << std::endl << std::flush;


  vnl_matrix<double> mat_b(3,3,0.0);
  itk::DiffusionTensor<double> tensor_2;
  tensor_2.SetIdentity();

  utlPrintVar1(true, tensor.EuclideanDistance(tensor_2)); 
  utlPrintVar1(true, tensor.KLDistance(tensor_2)); 
  utlPrintVar1(true, tensor.GeodesicDistance(tensor_2)); 
  utlPrintVar1(true, tensor.LogEucDistance(tensor_2)); 

  typedef itk::Vector<double,3> VectorType;
  typedef itk::Matrix<double,3,3> MatrixType;

  VectorType vec1, vec2;
  vec1[0]=vec1[1]=vec1[2]=1.0;
  vec2[0]=vec2[1]=vec2[2]=-2.0;
  vec2[2]=3.0;

  MatrixType matRotation;
  utl::RotationMatrixFromVectors<VectorType, MatrixType>(vec1, vec2, matRotation);
  std::cout << "rotation = \n" << matRotation << std::endl << std::flush;



  return 0;
}
