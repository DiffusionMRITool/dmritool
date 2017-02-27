/**
 *       @file  utl4thOrderTensorGTest.cxx
 *      @brief  
 *     Created  "06-21-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "utlGTest.h"
#include "utlNDArray.h"
#include "utlVNL.h"
#include "utlDMRI.h"


typedef utl::NDArray<double,1> UtlVectorType;
typedef utl::NDArray<double,2> UtlMatrixType;
typedef utl::NDArray<double,4> Utl4thTensorType;

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const T val1, const T val2 )
{
  vnl_vector<double> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
}

template <class T>
UtlVectorType
__GenerateUtlVector(const int N, const T val1, const T val2)
{
  UtlVectorType v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
}

template <class T>
UtlMatrixType
__GenerateUtlMatrix(const int rows, const int cols, const T val1, const T val2)
{
  UtlMatrixType mat(rows, cols);
  for ( int i = 0; i < rows*cols; i += 1 ) 
    mat[i] = utl::Random<double>(val1,val2);
  return mat;
}

TEST(utl4thOrderTensor, Constructors)
{
  int Dim = 4;
  unsigned shape[4];
  shape[0]=2, shape[1]=3, shape[2]=4, shape[3]=5;
  
  int N = 1;
  for ( int i = 0; i < Dim; ++i ) 
    N *= shape[i];
  
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
    
  Utl4thTensorType tensor(vec0.data_block(), shape), tensor2(shape);
  // utl::PrintUtlNDArray(tensor, "tensor");
  EXPECT_NEAR_VECTOR(vec0, tensor, N, 1e-10);
  tensor2=utl::VnlVectorToStdVector(vec0);
  EXPECT_NEAR_VECTOR(vec0, tensor2, N, 1e-10);
}

TEST(utl4thOrderTensor, Symmetrize)
{
  int Dim = 4;
  unsigned shape[4];
  shape[0]=3, shape[1]=3, shape[2]=3, shape[3]=3;
  
  int N = 1;
  for ( int i = 0; i < Dim; ++i ) 
    N *= shape[i];
  
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);

    {
    Utl4thTensorType tensor(vec0.data_block(), shape);
    tensor.MajorSymmetrize();
    EXPECT_TRUE(tensor.IsMajorSymmetric());
    }

    {
    Utl4thTensorType tensor(vec0.data_block(), shape);
    tensor.MinorSymmetrize();
    EXPECT_TRUE(tensor.IsMinorSymmetric());
    }
}

TEST(utl4thOrderTensor, Product)
{
  unsigned shape[4];
  shape[0]=2, shape[1]=3, shape[2]=4, shape[3]=5;
  
  int N = 1;
  for ( int i = 0; i < 4; ++i ) 
    N *= shape[i];
  
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);

  UtlMatrixType mat1(vec0.data_block(), shape[0], shape[1]), mat2(vec0.data_block()+shape[0]*shape[1], shape[2], shape[3]), mat3, mat4;
  Utl4thTensorType tensor;
  // mat1.PrintWithIndex(std::cout<<"mat1");
  // mat2.PrintWithIndex(std::cout<<"mat2");

  utl::OuterProduct(mat1, mat2, tensor);

    {
    for ( int i = 0; i < shape[0]; ++i ) 
      for ( int j = 0; j < shape[1]; ++j ) 
        for ( int k = 0; k < shape[2]; ++k ) 
          for ( int l = 0; l < shape[3]; ++l ) 
            {
            EXPECT_NEAR(tensor(i,j,k,l), mat1(i,j)*mat2(k,l),1e-10) << i <<  ", "<< j << ", " << k << ", " << l;
            }

    utl::InnerProduct(tensor, mat2, mat4);
    mat4.Scale(1.0/mat2.GetSquaredTwoNorm());
    EXPECT_NEAR_UTLMATRIX(mat4,mat1, 1e-10);

    utl::InnerProduct(mat1, tensor, mat4);
    mat4.Scale(1.0/mat1.GetSquaredTwoNorm());
    EXPECT_NEAR_UTLMATRIX(mat4,mat2, 1e-10);

    EXPECT_NEAR(InnerProduct(tensor,tensor), mat2.GetSquaredTwoNorm()*mat1.GetSquaredTwoNorm(), 1e-10);
    EXPECT_NEAR(InnerProduct(tensor,tensor), tensor.GetSquaredTwoNorm(), 1e-10);
    }
}

TEST(utl4thOrderTensor, Convert2th4thOrderTensor)
{
  unsigned shape[4];
  shape[0]=shape[1]=shape[2]=shape[3]=3;

  UtlMatrixType mat(6,6), mat1, mat2, mat3;
  Utl4thTensorType tensor, tensor1;
  mat = __GenerateUtlMatrix<double>(6,6,-2.0,2.0);
//   mat = {
// 0.365901000000000,	0.503819000000000,	1.59671000000000 , 0.964604000000000 , -1.53631000000000 , -0.194382000000000, 
// 1.31128000000000 , -1.50478000000000 , 1.26309000000000  ,	0.816242000000000,	0.895654000000000,	0.795273000000000,
// -1.83091000000000,	-1.77465000000000,	-1.95480000000000,	-1.85136000000000,	1.44264000000000 , -1.80099000000000 ,
// 0.283206000000000,	-1.59728000000000,	1.75744000000000 , -0.544227000000000,	-0.96231100000000,	-1.81135000000000,
// -1.43618000000000,	-1.68439000000000,	0.196733000000000,	1.58803000000000 , 1.37835000000000	 , -1.99786000000000 ,
// 1.77427000000000 , -0.255747000000000,	0.505962000000000,	1.37098000000000 , -1.29114000000000 , 0.969649000000000 
//   };
  // utl::PrintUtlNDArray(mat, "mat");
  utl::Convert2To4Tensor(mat, tensor);

    {
    // convert back 
    utl::Convert4To2Tensor(tensor, mat1);
    EXPECT_NEAR_UTLMATRIX(mat, mat1, 1e-10);
    }

    {
    // keep inner product: tensor [*] tensor
    mat1 = __GenerateUtlMatrix<double>(6,6,-2.0,2.0);
    utl::Convert2To4Tensor(mat1, tensor1);
    EXPECT_NEAR(utl::InnerProduct(mat, mat1), utl::InnerProduct(tensor,tensor1), 1e-10);
    }
    
    {
    // keep inner product:  matrix [*] tensor [*] matrix
    mat1 = __GenerateUtlMatrix<double>(3,3,-2.0,2.0);
    mat2 = __GenerateUtlMatrix<double>(3,3,-2.0,2.0);
    mat1.Symmetrize();
    mat2.Symmetrize();

    UtlVectorType v1, v2;
    utl::Convert2To1Tensor(mat1, v1);
    utl::Convert2To1Tensor(mat2, v2);

      {
      // Outproduct(mat1, mat2) [*] tensor  == mat1 * tensor *mat2
      utl::OuterProduct(mat1, mat2, tensor1);
      utl::InnerProduct(tensor, mat2, mat3);
      EXPECT_NEAR(utl::InnerProduct(tensor, tensor1), utl::InnerProduct(mat1, mat3), 1e-10);
      utl::InnerProduct(mat1, tensor, mat3);
      EXPECT_NEAR(utl::InnerProduct(tensor, tensor1), utl::InnerProduct(mat3, mat2), 1e-10);
      }

      {
      // mat1 [*] tensor [*] mat3  == v1 * mat *v2
      utl::InnerProduct(mat1, tensor, mat3);
      EXPECT_NEAR(utl::InnerProduct(mat3, mat2), utl::InnerProduct(v1, mat*v2), 1e-10);
      }
    }
  
    {
    // keep inner product: tensor [*] matrix
    mat1 = __GenerateUtlMatrix<double>(3,3,-2.0,2.0);
    mat1.Symmetrize();

    UtlVectorType v1, v2, v3;
    utl::Convert2To1Tensor(mat1, v1);

    utl::InnerProduct(tensor, mat1, mat2);
    utl::InnerProduct(mat, v1, v2);
    utl::Convert2To1Tensor(mat2, v3);
    EXPECT_NEAR_UTLVECTOR(v2, v3, 1e-10);
    }

    {
    // keep eigen-decomposition.
    UtlVectorType valReal, valImg, valReal0, valImg0;
    // UtlMatrixType vecRealR0, vecImgR0, vecRealL0, vecImgL0;

    // mat.EigenDecompositionNonSymmetricMatrix(valReal0, valImg0, vecRealR0, vecImgR0);
    // mat.EigenDecompositionNonSymmetricMatrix(valReal0, valImg0, vecRealR0, vecImgR0, vecRealL0, vecImgL0);

    std::vector<UtlMatrixType> matRealR, matImgR, matRealL, matImgL;

    tensor.EigenDecompositionWithMinorSymmetry(valReal, valImg, matRealR, matImgR);
    // utl::PrintUtlVector(valReal, "vecReal");
    // utl::PrintUtlVector(valImg, "vecImg");
    // for ( int i = 0; i < valReal.Size(); ++i ) 
    //   {
    //   std::cout << "i = " << i << std::endl << std::flush;
    //   utl::PrintUtlMatrix(matRealR[i], "matRealR");
    //   utl::PrintUtlMatrix(matImgR[i], "matImgR");
    //   }
    for ( int i = 0; i < valReal.Size(); ++i ) 
      {
      utl::InnerProduct(tensor, matRealR[i], mat1);
      EXPECT_NEAR(utl::ToMatrix<double>(mat1-matRealR[i]%valReal[i] + matImgR[i]%valImg[i])->GetTwoNorm(), 0, 1e-8 );
      utl::InnerProduct(tensor, matImgR[i], mat1);
      EXPECT_NEAR(utl::ToMatrix<double>(mat1-matRealR[i]%valImg[i] - matImgR[i]%valReal[i])->GetTwoNorm(), 0, 1e-8 );
      }
    
    tensor.EigenDecompositionWithMinorSymmetry(valReal, valImg, matRealR, matImgR, matRealL, matImgL);
    typedef std::complex<double> VT;
    utl::NDArray<VT,2> matR, matTmp, matL, matLH;
    for ( int i = 0; i < valReal.Size(); ++i ) 
      {
      utl::InnerProduct(tensor, matRealR[i], mat1);
      EXPECT_NEAR(utl::ToMatrix<double>(mat1-matRealR[i]%valReal[i] + matImgR[i]%valImg[i])->GetTwoNorm(), 0, 1e-8 );
      utl::InnerProduct(tensor, matImgR[i], mat1);
      EXPECT_NEAR(utl::ToMatrix<double>(mat1-matRealR[i]%valImg[i] - matImgR[i]%valReal[i])->GetTwoNorm(), 0, 1e-8 );

      utl::InnerProduct(matRealL[i], tensor, mat1);
      EXPECT_NEAR(utl::ToMatrix<double>(mat1-matRealL[i]%valReal[i] - matImgL[i]%valImg[i])->GetTwoNorm(), 0, 1e-8 );
      utl::InnerProduct(matImgL[i], tensor, mat1);
      EXPECT_NEAR(utl::ToMatrix<double>(mat1+matRealL[i]%valImg[i] - matImgL[i]%valReal[i])->GetTwoNorm(), 0, 1e-8 );

      //  Inner(A, V) = \lambda V
      matR = utl::ComplexCombine(matRealR[i], matImgR[i]);
      utl::InnerProduct(utl::ComplexCombine(tensor,0.0), matR, matTmp);
      EXPECT_NEAR(utl::ToMatrix<VT>(matTmp-matR%VT(valReal[i], valImg[i]))->GetTwoNorm(), 0, 1e-8 );

      //  Inner(Conj(U), A) = \lambda Conj(U)
      matL = utl::ComplexCombine(matRealL[i], matImgL[i]);
      matLH = utl::Conj(matL);
      utl::InnerProduct(matLH, utl::ComplexCombine(tensor,0.0), matTmp);
      EXPECT_NEAR(utl::ToMatrix<VT>(matTmp-matLH%VT(valReal[i], valImg[i]))->GetTwoNorm(), 0, 1e-8 );
      }

    // A V = V D
    {
    utl::NDArray<VT,2> matR(6,6), matRInv, matRe, D(6,6,VT(0.0));
    for ( int i = 0; i < 6; ++i ) 
      {
      utl::NDArray<double,1> vecR1, vecR2;
      utl::Convert2To1Tensor(matRealR[i], vecR1);
      utl::Convert2To1Tensor(matImgR[i], vecR2);
      for ( int j = 0; j < 6; ++j ) 
        matR(j,i) = VT( vecR1[j], vecR2[j]);
      D(i,i) = VT(valReal[i], valImg[i]);
      }
    matRInv = matR.InverseMatrix();
    matRe = matR * D * matRInv;
    EXPECT_NEAR_MATRIX_COMPLEX(mat, matRe, 6, 6, 1e-10);
    }
    
    // U^H  A  = D U^H
    {
    utl::NDArray<VT,2> matL(6,6), matLInv, matRe, D(6,6,VT(0.0));
    for ( int i = 0; i < 6; ++i ) 
      {
      utl::NDArray<double,1> vecL1, vecL2;
      utl::Convert2To1Tensor(matRealL[i], vecL1);
      utl::Convert2To1Tensor(matImgL[i], vecL2);
      for ( int j = 0; j < 6; ++j ) 
        matL(j,i) = VT( vecL1[j], vecL2[j]);
      D(i,i) = VT(valReal[i], valImg[i]);
      }
    matL = matL.GetConjugateTranspose();
    matLInv = matL.InverseMatrix();
    matRe = matLInv * D * matL;
    EXPECT_NEAR_MATRIX_COMPLEX(mat, matRe, 6, 6, 1e-10);
    }

    }
}
