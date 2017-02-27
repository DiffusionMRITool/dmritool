/**
 *       @file  itkSHCoefficientsRotationGTest.cxx
 *      @brief  
 *     Created  "02-26-2015
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "utlGTest.h"
#include "itkSHCoefficientsRotation.h"
#include "utl.h"
#include "utlRotationMatrixFromVectors.h"

typedef itk::SHCoefficientsRotation<double> SHRotationFilter;
typedef utl::NDArray<double,1> VectorType;
typedef utl::NDArray<double,2> MatrixType;
typedef utl_shared_ptr<MatrixType> MatrixPointer;

inline VectorType
__GenerateRandomSH(const int rank)
{
  int dim = utl::RankToDimSH(rank);
  VectorType vec(dim);
  for ( int i = 0; i < dim; ++i ) 
    vec[i] = utl::Random<double>(-10.0,10.0);
  return vec;
}

inline MatrixType
__GeterateRandomRotation()
{
  VectorType v1(3), v2(3);
  for ( int i = 0; i < 3; ++i ) 
    {
    v1[i] = utl::Random<double>(-1.0,1.0);
    v2[i] = utl::Random<double>(-1.0,1.0);
    }
  MatrixType mat(3,3);
  utl::RotationMatrixFromVectors(v1, v2, mat);
  return mat;
}

#define __Test_RandomRotate(order)                                                                 \
do                                                                                                 \
{                                                                                                  \
  SHRotationFilter::Pointer shRotate = SHRotationFilter::GetInstance();                            \
  shRotate->SetTessOrder(order);                                                                   \
  shRotate->SetMaxRank(10);                                                                        \
  shRotate->Initialize();                                                                          \
  VectorType shVec_2 = shRotate->GetRotatedSHCoefficients(shVec, rotateMatrix.GetTranspose());     \
  VectorType sfVec_2 = (*shMatrix) *shVec_2;                                                       \
                                                                                                   \
  VectorType sfdiff = sfVec - sfVec_2;                                                             \
  double diffNorm = sfdiff.GetTwoNorm();                                                           \
  double norm = sfVec.GetTwoNorm();                                                                \
  EXPECT_NEAR(diffNorm/norm, 0, 1e-6);                                                             \
  EXPECT_NEAR_VECTOR(sfVec, sfVec_2, shVec.Size(), mean*1e-4);                                     \
} while ( 0 )

TEST(itkSHCoefficientsRotation, RandomRotate)
{
  int rank = 8;
  MatrixPointer gradT4 = utl::ReadGrad<double>(4, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  MatrixType rotateMatrix = __GeterateRandomRotation();

  MatrixPointer shMatrix = utl::ComputeSHMatrix(rank, *gradT4, CARTESIAN_TO_SPHERICAL);
  MatrixPointer gradT4Rotate(new MatrixType(gradT4->Rows(), gradT4->Cols()));
  *gradT4Rotate = (rotateMatrix * gradT4->GetTranspose()).GetTranspose();
  MatrixPointer shMatrixRotate = utl::ComputeSHMatrix(rank, *gradT4Rotate, CARTESIAN_TO_SPHERICAL);

  VectorType shVec = __GenerateRandomSH(rank);
  VectorType sfVec = (*shMatrixRotate) * shVec;
  double mean = sfVec.GetElementAbsolute().GetMean();

  __Test_RandomRotate(3);
  __Test_RandomRotate(4);
  __Test_RandomRotate(5);

  // EXPECT_NEAR_VECTOR(sf, sfVec_2, shVec.Size(), 1e-8);

  // utl::PrintUtlVector(sfVec, "sfVec");
  // utl::PrintUtlVector(sfVec_2, "sfVec_2");
}

